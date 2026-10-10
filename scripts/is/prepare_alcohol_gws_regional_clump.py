#!/usr/bin/env python3
"""Build full Japanese exposure GWAS regional clump inputs for 1000G EAS PLINK2.

Reads every original 7.68M GWAS row, filters ALL p<1e-4 within each of
seven ±500kb windows (includes all GWS within windows), exact GRCh37
allele pair to 1000G EAS pvar, allows REF/ALT reversed, excludes non-SNV,
palindromic and allele-discordant rows. Original numeric p=0 is stored
unmodified in audit, PLINK-safe p=1e-300. No fake novel independent IVs.
"""
import argparse,csv,gzip,json,math
from collections import Counter,defaultdict
from pathlib import Path
from extract_alcohol_eas_region_reference_for_clump import ROOT as REFROOT,SENTINELS
SRC=Path("/srv/is-analysis/data/is/gwas_exposure/japanese_alcohol_2024/1_Alcohol_intake_Unstratified.tsv.gz")
OUT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_gws_regional_clumping")
def normalized(chr_,pos,ref,alt):
    return (chr_.removeprefix("chr"),int(pos),ref.upper(),alt.upper())
def same_pair(a,b,c,d):
    return (a,b)==(c,d) or (a,b)==(d,c)
def valid_snv(*alleles):
    return all(len(x)==1 and x in "ACGT" for x in alleles)
def is_palindrome(ref,alt):
    return {ref,alt} in ({"A","T"},{"C","G"})
def make_reference(chrom,refroot):
    pvar=refroot/f"chr{chrom}.EAS504.regions.pvar"
    source=json.loads((refroot/f"chr{chrom}.region_source_audit.json").read_text())
    if source["EAS_samples"]!=504 or source["plink_sample_count"]!=504 or not source["validated"]:
        raise ValueError("Raw EAS sample and region source provenance invalid")
    mapping=defaultdict(list)
    with pvar.open() as f:
        for row in csv.DictReader((line for line in f if not line.startswith("##")),delimiter="\t"):
            chromosome=row.get("#CHROM") or row.get("CHROM")
            if chromosome!=chrom:raise ValueError("Mixed reference chromosomes")
            ref=row["REF"].upper();alt=row["ALT"].upper()
            if not valid_snv(ref,alt):continue
            mapping[(chrom,int(row["POS"]))].append((row["ID"],ref,alt))
    if len(mapping)<100:raise ValueError("Reference PLINK SNPs unexpectedly missing")
    return mapping,source
def make(root,refroot,source,source_min_rows=7000000):
    refs={c:make_reference(c,refroot)[0] for c in SENTINELS}
    regions={c:[(p-500000,p+500000) for p in positions]
             for c,positions in SENTINELS.items()}
    accepted=defaultdict(list);stats={c:Counter() for c in SENTINELS}
    with gzip.open(source,"rt") as f:
        reader=csv.DictReader(f,delimiter="\t")
        assert reader.fieldnames==["SNP","CHR","POS","EA","NEA","EAF","BETA","SE","P","HetP","N"]
        total=0
        for x in reader:
            total+=1
            c=x["CHR"].removeprefix("chr")
            if c not in refs:continue
            pos=int(x["POS"])
            if not any(lo<=pos<=hi for lo,hi in regions[c]):continue
            stats[c]["all_rows_in_500kb_windows"]+=1
            p=float(x["P"])
            if not math.isfinite(p) or p<0 or p>1:raise ValueError("Nonfinite source p-value")
            if p>=1e-4:continue
            stats[c]["source_p_lt_1e4"]+=1
            if p<5e-8:stats[c]["source_p_lt_5e8"]+=1
            parts=x["SNP"].split("_")
            if len(parts)!=4 or parts[0].removeprefix("chr")!=c or parts[1]!=str(pos):
                stats[c]["id_chr_pos_mismatch"]+=1;continue
            ref,alt=parts[2].upper(),parts[3].upper()
            ea=x["EA"].upper();nea=x["NEA"].upper()
            if not valid_snv(ref,alt,ea,nea) or not same_pair(ref,alt,ea,nea):
                stats[c]["non_snv_or_ea_pair_discordance"]+=1;continue
            if is_palindrome(ref,alt):
                stats[c]["palindrome_removed"]+=1;continue
            matches=[(ident,ref0,alt0) for ident,ref0,alt0 in refs[c].get((c,pos),[])
                     if same_pair(ref,alt,ref0,alt0)]
            if len(matches)!=1:
                stats[c]["not_single_exact_biallelic_reference_match"]+=1;continue
            rid,rref,ralt=matches[0]
            beta=float(x["BETA"]);se=float(x["SE"]);eaf=float(x["EAF"])
            if not all(map(math.isfinite,(beta,se,eaf))) or se<=0 or not(0<eaf<1):
                stats[c]["bad_beta_se_eaf"]+=1;continue
            # Study ALT alignment and reference ALT alignment are separately stored.
            src_alt_beta=beta if ea==alt else -beta
            ref_alt_beta=src_alt_beta if alt==ralt else -src_alt_beta
            ref_alt_eaf=(eaf if ea==ralt else 1-eaf)
            accepted[c].append({"ID":rid,"P":max(p,1e-300),
                "P_original":p,"P_clamped_from_zero":int(p==0),
                "chr":c,"pos":pos,"reference_ALT":ralt,"GWAS_ALT":alt,
                "reference_ALT_beta":ref_alt_beta,"reference_ALT_eaf_Japanese":ref_alt_eaf,
                "source_neffect":int(float(x["N"])),
                "ref_alt_swapped_source":int(rref!=ref),
                "source_allele_A1":ea,"source_p":p})
            stats[c]["exact_reference_allele_matches"]+=1
        if total<source_min_rows:raise ValueError("Unexpected truncated Japanese source")
    root.mkdir(parents=True,exist_ok=True)
    all_associations=[]
    for c,items in accepted.items():
        dedup={}
        for r in items:
            if r["ID"] not in dedup or r["P"]<dedup[r["ID"]]["P"]:
                dedup[r["ID"]]=r
        items=list(dedup.values())
        stats[c]["clump_input_rows_after_dedup"]=len(items)
        stats[c]["clump_input_GWS"]=sum(x["P_original"]<5e-8 for x in items)
        if not any(x["P_original"]<5e-8 for x in items):
            raise ValueError("No matched Japanese GWS in "+c)
        for a in SENTINELS[c]:
            candidates=[x for x in items if x["pos"]==a]
            if len(candidates)!=1:raise ValueError(f"Missing original sentinel chr{c}:{a}")
        clump=root/f"chr{c}.alcohol_regional.clump_input.tsv"
        with clump.open("w",newline="") as f:
            w=csv.DictWriter(f,delimiter="\t",fieldnames=["ID","P"])
            w.writeheader();w.writerows({k:x[k] for k in ("ID","P")} for x in items)
        all_associations.extend(items)
    detailed=root/"IS_ALCOHOL_SOURCE_REF_MATCHED_REGIONAL_ASSOCIATIONS.tsv"
    with detailed.open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(all_associations[0]))
        w.writeheader();w.writerows(all_associations)
    report={"source_original_total_rows":total,
        "chrom_source_stats":{c:dict(v) for c,v in stats.items()},
        "original_source":"Koyanagi2024_Zenodo_10038152_MD5_VERIFIED",
        "reference":"1000G_phase3_v5b_EAS504_GRCh37",
        "regions":{c:regions[c] for c in regions},
        "p1_index_threshold":5e-8,"p2_secondary_threshold":1e-4,
        "p_zero_clamped_for_PLINK":sum(x["P_clamped_from_zero"] for x in all_associations),
        "regional_clumping_only":True,
        "genome_wide_clumping_claimed":False,
        "horizontal_pleiotropy_excluded":False,
        "sample_independence_verified":False}
    (root/"IS_ALCOHOL_REGIONAL_GWAS_REFERENCE_COVERAGE.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True)
    return report
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,default=SRC)
    p.add_argument("--reference",type=Path,default=REFROOT)
    p.add_argument("--out",type=Path,default=OUT)
    a=p.parse_args();make(a.out,a.reference,a.source)
