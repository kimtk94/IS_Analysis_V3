#!/usr/bin/env python3
"""Audit *why* four EAS AIS G0022 SuSiE CS SNPs are not in GTEx QTL.

Uses official UCSC hg19ToHg38 chain, reciprocal hg38ToHg19 position check,
full 720-kb GRCh38 nominal-QTL extracts for each of four tissues.
Does not equate absent source assay with a negative molecular association.
"""
import argparse,csv,gzip,json,math
from collections import Counter,defaultdict
from pathlib import Path
from audit_g0022_gtex_v8_dapg import LiftOver,complement,variant_liftover
from prepare_g0022_eqtl_catalogue_multitissue import ROOT,RAW,OUT,SOURCES,source_file
from prepare_g0022_eqtl_catalogue_aorta_coloc import parse_nominal
CHAIN_FORWARD=Path("/srv/is-analysis/data/is/reference/hg19ToHg38.over.chain.gz")
CHAIN_REVERSE=Path("/srv/is-analysis/data/reference/hg38ToHg19.over.chain")
def orient_alleles(snp,strand):
    _,pos,ref,alt=snp.split(":")
    if strand=="-":ref,alt=complement(ref),complement(alt)
    return ref,alt
def audit(root,raw,out,forward,reverse):
    if LiftOver is None:raise RuntimeError("pyliftover not installed for forward chain")
    model=json.loads((root/"G0022_FULL_AIS_SUSIE_SUMMARY.json").read_text())
    cs={s for lst in model["credible_set_variant_ids"].values() for s in lst}
    if len(cs)!=4:raise ValueError("Unexpected full-locus credible-set source")
    with (root/"variants.tsv").open() as f:
        variants={x["variant_id"]:x for x in csv.DictReader(f,delimiter="\t") if x["variant_id"] in cs}
    if set(variants)!=cs:raise ValueError("All four CS variants must have EAS allele QC")
    lo=LiftOver(str(forward));back=LiftOver(str(reverse))
    lifted={}
    for snp in cs:
        contig,pos,ref,alt=snp.split(":")
        if contig!="12":raise ValueError("Unexpected chromosome")
        result=lo.convert_coordinate("chr12",int(pos)-1)
        if len(result)!=1:raise ValueError("Ambiguous forward liftOver: "+snp)
        chrom,new0,strand,_=result[0]
        if chrom!="chr12" or strand not in ("+","-"):
            raise ValueError("Unexpected mapped contig/strand")
        rr=back.convert_coordinate(chrom,int(new0))
        # Reference mapping can be ambiguous for parts of genome. Require reciprocal.
        if len(rr)!=1 or rr[0][0]!="chr12" or int(rr[0][1])!=int(pos)-1:
            raise ValueError("Unstable reciprocal liftOver: "+snp)
        target_ref,target_alt=orient_alleles(snp,strand)
        lifted[snp]=(int(new0)+1,target_ref,target_alt,strand)
    rows=[];dataset_summary=[]
    for study in SOURCES:
        qtlfile=source_file(raw,study)
        if not qtlfile.is_file():raise FileNotFoundError(qtlfile)
        wanted={x[0] for x in lifted.values()}
        index=defaultdict(list);count=0
        with gzip.open(qtlfile,"rt") as f:
            for line in f:
                if not line.strip():continue
                a=parse_nominal(line);count+=1
                if int(a["position"]) in wanted:index[int(a["position"])].append(a)
        tissue,sample_size=SOURCES[study]
        row_count=0
        for snp,(position,ref,alt,strand) in sorted(lifted.items(),key=lambda x:x[1][0]):
            hits=index.get(position,[])
            exact=[h for h in hits if h["ref"]==ref and h["alt"]==alt]
            swaps=[h for h in hits if h["ref"]==alt and h["alt"]==ref]
            same_pair=[h for h in hits if set((h["ref"],h["alt"]))=={ref,alt}]
            query=variants[snp]
            state=("NOT_TESTED_VARIANT_POSITION_ABSENT" if not hits else
                "VARIANT_POSITION_PRESENT_DIFFERENT_ALLELES" if not same_pair else
                "EXACT_REF_ALT_QTL_VARIANT_ASSAYED" if exact else
                "SWAPPED_QTL_REF_ALT_REQUIRES_GENOME_RECHECK")
            if exact:row_count+=1
            rows.append({
              "dataset_id":study,"tissue":tissue,
              "gwas_CS_variant_GRCh37":snp,"QTL_variant_position_GRCh38":position,
              "expected_ref_GRCh38":ref,"expected_alt_GRCh38":alt,
              "liftover_strand":strand,
              "gwas_alt_eaf_EAS":query["gwas_alt_eaf"],
              "gwas_eas_reference_alt_af":query["reference_alt_af"],
              "qtl_records_at_position_any_gene":len(hits),
              "qtl_records_exact_ref_alt":len(exact),
              "qtl_records_swapped_ref_alt":len(swaps),
              "source_alt_eaf_if_exact":";".join(sorted(set(h["ac"]+"/"+h["an"] for h in exact))),
              "source_status":state,
              "interpretation":"ASSAY_COVERAGE_AUDIT_NOT_MOLECULAR_EFFECT_PROOF"})
        dataset_summary.append({"dataset_id":study,"tissue":tissue,
             "regional_association_rows":count,"credible_set_positions_assayed_exactly":row_count})
    out.mkdir(parents=True,exist_ok=True)
    with (out/"G0022_FULL_AIS_FOUR_CS_QTL_SOURCE_COVERAGE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
    if len(rows)!=len(SOURCES)*4:raise RuntimeError("Incomplete source audit")
    summary={"source":"OFFICIAL_UCSC_HG19_TO_HG38_PLUS_GTEX_EQTL_CATALOGUE",
       "dataset_coverage":dataset_summary,"credible_set_snps":len(cs),
       "dataset_variant_tests":len(rows),"by_coverage_status":dict(Counter(x["source_status"] for x in rows)),
       "status":"SOURCE_ASSAY_COVERAGE_AUDITED_NO_BIOLOGICAL_NEGATIVE",
       "note":"Absence from tested source SNP set prevents colocalization but does not establish absence of eQTL"}
    (out/"G0022_FULL_AIS_FOUR_CS_QTL_SOURCE_COVERAGE_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    for r in rows:
        print(r["dataset_id"],r["gwas_CS_variant_GRCh37"],
              r["QTL_variant_position_GRCh38"],r["source_status"])
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--raw",type=Path,default=RAW)
    p.add_argument("--out",type=Path,default=OUT)
    p.add_argument("--forward",type=Path,default=CHAIN_FORWARD)
    p.add_argument("--reverse",type=Path,default=CHAIN_REVERSE)
    a=p.parse_args()
    audit(a.root,a.raw,a.out,a.forward,a.reverse)
