#!/usr/bin/env python3
"""Verify all 8,956 GRCh38 coloc sites against native BBJ IS GWAS and hg19->hg38.

READ-ONLY regarding input files. The native Japanese Sakaue/Kanai 2020
autosomal ZIP member documents AF_Allele2 / BETA. Exact original v
identities, effect allele, beta, EAF, imputation quality and reference
coordinate conversions are compared without assuming non-palindromic status.

Native GWAS effect-allele provenance does NOT verify GTEx allele-count
preprocessing, donor ancestry match, GTEx genotype LD, or causal colocalization.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math,zipfile,sys
from collections import Counter,defaultdict
from pathlib import Path

REQUIRED=("locus","dataset_key","gene_base","variant_id","match_key",
          "ref","alt","gwas_beta","gwas_eaf","harmonization")
def fnum(x):
    n=float(x)
    if not math.isfinite(n):raise ValueError("Non-finite numeric")
    return n
def pass_value(expected,obs,tol=1e-11):
    return abs(expected-obs) <= tol*max(1,abs(expected),abs(obs))
def sha_file(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
    return h.hexdigest()
def collect(recovered_dir,master):
    keys={(r["locus"],r["dataset_key"],r["gene_base"]):int(r["nsnps"])
          for r in master}
    if len(keys)!=len(master) or len(master)!=646:
        raise ValueError("646 unique assays required")
    per_site={};count=0;conflicts=[]
    for i,(locus,ds,gene) in enumerate(sorted(keys)):
        name=f"{locus}__{ds}__{gene}.tsv"
        path=recovered_dir/name
        if not path.is_file():raise FileNotFoundError(path)
        n=0
        with path.open(newline="") as f:
            rd=csv.DictReader(f,delimiter="\t")
            if not set(REQUIRED).issubset(rd.fieldnames or []):
                raise ValueError("Required source columns missing in "+name)
            for row in rd:
                n+=1;count+=1
                if (row["locus"],row["dataset_key"],row["gene_base"])!=(locus,ds,gene):
                    raise ValueError("Assay association key differs")
                if row["harmonization"]!="EXACT_REF_ALT_GRCH38":
                    raise ValueError("Unexpected harmonization")
                native_key=row["variant_id"]
                target=row["match_key"]
                nseg=native_key.split(":");tseg=target.split(":")
                if len(nseg)!=4 or len(tseg)!=4 or tseg[2:]!=nseg[2:]:
                    raise ValueError("Source genome-build/REF:ALT mismatch")
                if row["ref"]!=tseg[2] or row["alt"]!=tseg[3]:
                    raise ValueError("GTEx reference-alleles mismatched")
                if len(tseg[2])!=1 or len(tseg[3])!=1:
                    raise ValueError("Explicitly exclude unsupported indels in audit")
                meta={"locus":locus,"variant_id_GRCh37":native_key,
                    "match_key_GRCh38":target,"gwas_beta_cached":fnum(row["gwas_beta"]),
                    "gwas_eaf_cached":fnum(row["gwas_eaf"]),
                    "ref_GRCh38":row["ref"],"alt_GRCh38":row["alt"],"count_assay_rows":0,
                    "per_assay_tissues":set(),"cache_consistency_status":"PASS"}
                key=(locus,target)
                if key not in per_site:per_site[key]=meta
                p=per_site[key]
                if (p["variant_id_GRCh37"]!=native_key
                    or not pass_value(p["gwas_beta_cached"],meta["gwas_beta_cached"])
                    or not pass_value(p["gwas_eaf_cached"],meta["gwas_eaf_cached"])):
                    if len(conflicts)<30:conflicts.append({"key":key,"new_variant":native_key})
                    p["cache_consistency_status"]="SOURCE_REUSED_GWAS_INCONSISTENT"
                p["count_assay_rows"]+=1
                p["per_assay_tissues"].add(ds)
        if n!=keys[(locus,ds,gene)]:
            raise ValueError(f"SNP count {n} disagrees with archived master {name}")
        if i%100==99:print("COLLECT",i+1,"/646",flush=True)
    if len(per_site)!=8956:
        raise ValueError(f"Expected 8956 unique locus SNPs, got {len(per_site)}")
    if any(x["cache_consistency_status"]!="PASS" for x in per_site.values()):
        raise ValueError(f"Cached GWAS metrics inconsistent across assays: {conflicts}")
    native_set={x["variant_id_GRCh37"] for x in per_site.values()}
    print("COLLECT_COMPLETE","source_rows",count,"unique_locus_variants",len(per_site),
          "unique_native_variants",len(native_set),flush=True)
    return per_site,native_set,count
def stream_native(native_zip,native_set):
    found={};duplicates=set();records=0
    with zipfile.ZipFile(native_zip) as archive:
        member=[x for x in archive.namelist() if x.endswith(".auto.txt.gz")]
        if len(member)!=1:raise ValueError("Expected exactly one original BBJ autosomal GWAS")
        with archive.open(member[0]) as bytes_stream:
            with gzip.open(bytes_stream,"rt",newline="") as text:
                rd=csv.DictReader(text,delimiter="\t")
                fields={"v","CHR","POS","Allele1","Allele2","AF_Allele2","BETA"}
                if not fields.issubset(rd.fieldnames or []):
                    raise ValueError("Native allele coding not as documented")
                for row in rd:
                    records+=1
                    vid=row["v"]
                    if vid in native_set:
                        if vid in found:duplicates.add(vid)
                        else:found[vid]={
                            "source_allele1":row["Allele1"].upper(),
                            "source_allele2":row["Allele2"].upper(),
                            "source_beta":fnum(row["BETA"]),
                            "source_AF_Allele2":fnum(row["AF_Allele2"]),
                            "source_v":vid,
                            "source_position":row["POS"],
                            "source_chromosome":row["CHR"]
                        }
                    if records%4000000==0:
                        print("NATIVE_STREAM",records,"matched",len(found),flush=True)
    print("NATIVE_STREAM_COMPLETE","rows",records,"found",len(found),"duplicates",len(duplicates),flush=True)
    return found,duplicates,records
def audit(sites,native,duplicates,liftover):
    rows=[]
    for k,r in sorted(sites.items()):
        parts=r["variant_id_GRCh37"].split(":")
        chrom,coord,ref,alt=parts
        source=native.get(r["variant_id_GRCh37"])
        dest=r["match_key_GRCh38"].split(":")
        notes=[];status="PASS_NATIVE_EFFECT_BETA_AF_AND_CHAIN"
        allele_side=""
        lifted=liftover.convert_coordinate("chr"+chrom,int(coord)-1)
        if len(lifted)!=1:
            notes.append("CHAIN_NON_UNIQUE")
        else:
            new_chrom,pos,strand,*_ = lifted[0]
            if new_chrom.removeprefix("chr")!=dest[0] or int(pos)+1!=int(dest[1]) or strand!="+":
                notes.append("CHAIN_BUILD_STRAND_OR_POS_MISMATCH")
        if source is None:
            notes.append("NATIVE_VARIANT_NOT_FOUND")
        elif r["variant_id_GRCh37"] in duplicates:
            notes.append("NATIVE_DUPLICATE_ID")
        else:
            if (source["source_chromosome"]!=chrom or
                source["source_position"]!=coord or
                source["source_allele1"]!=ref or
                source["source_allele2"]!=alt):
                notes.append("NATIVE_ALLELE_REFALT_MISMATCH")
            allele_side=("ALT" if source["source_allele2"]==alt
                          else "REF" if source["source_allele2"]==ref else "OTHER")
            if not pass_value(source["source_beta"],r["gwas_beta_cached"]):
                notes.append("CACHED_GWAS_BETA_NOT_NATIVE")
            if not pass_value(source["source_AF_Allele2"],r["gwas_eaf_cached"]):
                notes.append("CACHED_GWAS_EAF_NOT_NATIVE")
        if notes:status="REVIEW"
        rows.append({
            "locus":r["locus"],"variant_id_GRCh37":r["variant_id_GRCh37"],
            "variant_id_GRCh38":r["match_key_GRCh38"],
            "gwas_cached_effect_beta":r["gwas_beta_cached"],
            "gwas_cached_EAF":r["gwas_eaf_cached"],
            "native_effect_allele_Allele2":source["source_allele2"] if source else "",
            "native_other_allele_Allele1":source["source_allele1"] if source else "",
            "native_effect_vs_GRCh37_alt":allele_side,
            "native_beta":source["source_beta"] if source else "",
            "native_AF_Allele2":source["source_AF_Allele2"] if source else "",
            "lift_target":(str(lifted[0][0])+":"+str(int(lifted[0][1])+1)) if len(lifted)==1 else "",
            "lift_strand":lifted[0][2] if len(lifted)==1 else "",
            "n_gene_tissue_rows_reusing_SNP":r["count_assay_rows"],
            "n_GTEx_tissues_containing_SNP":len(r["per_assay_tissues"]),
            "status":status,"issues":";".join(notes),
            "GTEx_effect_allele_not_independently_checked":True,
            "scientific_status":"BBJ_EFFECT_DIRECTION_SOURCE_AUDIT_NOT_CAUSAL_GENE",
        })
    return rows
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",required=True,type=Path)
    p.add_argument("--master",required=True,type=Path)
    p.add_argument("--native-zip",type=Path,default=Path("/srv/is-analysis/data/is/east_asia/japan/bbj/hum0197.v3.BBJ.IS.v1.zip"))
    p.add_argument("--chain",type=Path,default=Path("/srv/is-analysis/data/is/reference/hg19ToHg38.over.chain.gz"))
    p.add_argument("--vendor-pyliftover",type=Path,default=Path("/srv/is-analysis/data/is/reference/vendor_pyliftover"))
    p.add_argument("--out-dir",required=True,type=Path)
    a=p.parse_args()
    for x in (a.native_zip,a.master,a.chain):
        if not x.is_file():raise FileNotFoundError(x)
    if a.out_dir.exists() and any(a.out_dir.iterdir()):
        raise FileExistsError("Refuse overwrite")
    with a.master.open(newline="") as f:master=list(csv.DictReader(f,delimiter="\t"))
    sites,keys,n_source_rows=collect(a.source_dir,master)
    native,dups,n_native=stream_native(a.native_zip,keys)
    sys.path.insert(0,str(a.vendor_pyliftover))
    from pyliftover import LiftOver
    lo=LiftOver(str(a.chain))
    rows=audit(sites,native,dups,lo)
    if len(rows)!=8956:raise ValueError("Full 8956 SNP source audit scope changed")
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"IS_8956_BBJ_NATIVE_EFFECT_ALLELE_CHAIN_QC.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
    counts=Counter(x["status"] for x in rows)
    tags=Counter(tag for x in rows for tag in x["issues"].split(";") if tag)
    side=Counter(x["native_effect_vs_GRCh37_alt"] for x in rows)
    report={
     "status":"NATIVE_BBJ_GWAS_EFFECT_ALLELE_AND_CHAIN_SNP_AUDIT",
     "source_GWAS_zip":str(a.native_zip),
     "source_GWAS_size_bytes":a.native_zip.stat().st_size,
     "native_autosomal_GWAS_rows_streamed":n_native,
     "source_coloc_gene_tissue_files":646,
     "source_assay_SNP_rows":n_source_rows,
     "distinct_GRCh37_native_ids":len(keys),
     "distinct_locus_GRCh38_SNP_ids":len(sites),
     "native_ids_found":len(native),
     "status_counts":dict(counts),
     "issue_counts":dict(tags),
     "native_effect_allele_vs_GRCh37_ALT":dict(side),
     "source_reference_file_sha256":{"master":sha_file(a.master),"chain":sha_file(a.chain)},
     "expected_native_effect_allele":"Allele2",
     "beta_and_EAF_compared_to_original_GWAS":True,
     "GRCh37_to_38_coordinate_and_plus_strand_checked":True,
     "GRCh38_actual_reference_bases_independently_checked":False,
     "GTEx_regression_effect_allele_source_provenance_confirmed_here":False,
     "population_frequency_discrepancy_is_strand_flip_proof":False,
     "causal_locus_conclusion":"NOT_ESTABLISHED",
     "no_originals_modified":True
    }
    (a.out_dir/"IS_8956_BBJ_NATIVE_PROVENANCE_SUMMARY.json").write_text(json.dumps(report,indent=2)+"\n")
    print("IS_ALL_8956_NATIVE_GWAS_CHAIN_AUDIT_COMPLETE",json.dumps({"rows":len(rows),"status":dict(counts),"allele_direction":dict(side),"issue_counts":dict(tags)}),flush=True)

if __name__=="__main__":main()
