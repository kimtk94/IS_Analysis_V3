#!/usr/bin/env python3
"""Independent GRCh37->GRCh38 coordinate and REF/ALT audit for TPMI 433.21.
Uses local UCSC hg19ToHg38 chain + public UCSC hg38 sequence, compares
Ensembl GRCh38 rsID coordinates. No TPMI GWAS download or source changes.
"""
import argparse,csv,json,sys,time,urllib.request
from collections import defaultdict,Counter
from pathlib import Path

UCSC="https://api.genome.ucsc.edu/getData/sequence"
ENSEMBL="https://rest.ensembl.org/variation/human/"
def get_json(url):
    req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"IS_Analysis_V3_GRCh38_Reference_QC/1.0"})
    with urllib.request.urlopen(req,timeout=40) as response:
        return json.load(response)

def intervals_for_targets(targets,lo):
    out=[]
    for r in targets:
        chrom="chr"+r["chrom"]
        mapped=lo.convert_coordinate(chrom,int(r["pos"])-1)
        if not mapped or len(mapped)!=1:
            raise ValueError("Liftover not unique: "+r["variant_id"])
        newchrom,p0,strand,*_=mapped[0]
        if newchrom!=chrom or strand!="+" or int(p0)!=p0:
            raise ValueError("Liftover different chrom/strand/coordinate: "+r["variant_id"])
        out.append({**r,"chr38":newchrom,"pos38":int(p0)+1})
    return out

def fetch_reference_bases(records):
    grouped=defaultdict(list)
    for r in records:grouped[r["chr38"]].append(r["pos38"])
    bases={}
    for chrom,locs in grouped.items():
        # UCSC API is zero-based half-open; sequence[0] is start.
        start=min(locs)-1
        end=max(locs)
        url=f"{UCSC}?genome=hg38;chrom={chrom};start={start};end={end}"
        resp=get_json(url)
        dna=resp.get("dna","").upper()
        if (resp.get("genome")!="hg38" or resp.get("chrom")!=chrom or
            int(resp.get("start",-1))!=start or int(resp.get("end",-1))!=end or
            len(dna)!=(end-start)):
            raise ValueError("Unexpected UCSC sequence metadata")
        for pos in locs:
            base=dna[pos-1-start]
            if base not in "ACGT":raise ValueError("Non-ACGT reference")
            bases[(chrom,pos)]=base
    return bases

def check_ensembl_rs(r,variation):
    rsid=r["allele_verified_rsids"]
    if variation.get("name")!=rsid:return "RSID_RESPONSE_MISMATCH"
    matching=[m for m in variation.get("mappings",[])
              if m.get("assembly_name")=="GRCh38" and
              m.get("seq_region_name")==r["chrom"] and
              m.get("start")==r["pos38"] and
              m.get("end")==r["pos38"] and m.get("strand")==1]
    if len(matching)!=1:return "ENSEMBL_MAPPING_MISMATCH_OR_DUPLICATE"
    alleles=set(str(matching[0].get("allele_string","")).upper().split("/"))
    if not {r["ref"].upper(),r["alt"].upper()}.issubset(alleles):
        return "ENSEMBL_ALLELE_SET_DISCORDANCE"
    return "ENSEMBL_RSID_POSITION_ALLELES_PASS"

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--chain",type=Path,default=Path("/srv/is-analysis/data/is/reference/hg19ToHg38.over.chain.gz"))
    p.add_argument("--vendor-root",type=Path,default=Path("/srv/is-analysis/data/is/reference/vendor_pyliftover"))
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    with a.input.open(newline="") as f:targets=list(csv.DictReader(f,delimiter="\t"))
    if len(targets)!=11 or len({x["variant_id"] for x in targets})!=11:
        raise ValueError("Expected 11 distinct GRCh37 CS variants")
    if any(x["allele_mapping_status"]!="VERIFIED_ALLELE_SET" for x in targets):
        raise ValueError("Upstream GRCh37 VEP audit incomplete")
    if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError(a.out_dir)
    sys.path.insert(0,str(a.vendor_root))
    from pyliftover import LiftOver
    lo=LiftOver(str(a.chain))
    records=intervals_for_targets(targets,lo)
    ref_bases=fetch_reference_bases(records)
    out=[]
    for r in records:
        ref38=ref_bases[(r["chr38"],r["pos38"])]
        source={r["ref"].upper(),r["alt"].upper()}
        status="CHAIN_UCSC_REF_MISMATCH"
        validated=""
        rs_status="NOT_QUERIED"
        if ref38 in source:
            alt38=(source-{ref38}).pop()
            validated=f"{r['chr38'][3:]}:{r['pos38']}:{ref38}:{alt38}"
            status="CHAIN_UCSC_GRCH38_REF_ALT_PASS"
            try:
                v=get_json(ENSEMBL+r["allele_verified_rsids"]+"?content-type=application/json")
                rs_status=check_ensembl_rs(r,v)
            except Exception as e:
                rs_status="ENSEMBL_API_UNAVAILABLE_"+type(e).__name__
        out.append({"gene":r["locus"],"rsid":r["allele_verified_rsids"],
                    "GRCh37":r["variant_id"],"GRCh38":validated,
                    "chain_position38":r["pos38"],"chain_strand":"+",
                    "ucsc_ref38":ref38,"reference_status":status,
                    "ensembl_rsid_audit":rs_status})
        print("MAP",r["allele_verified_rsids"],validated,status,rs_status,flush=True)
        time.sleep(0.1)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    dest=a.out_dir/"TPMI_11_SNP_GRCH38_CHAIN_UCSC_ENS_AUDIT.tsv"
    with dest.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter="\t")
        w.writeheader();w.writerows(out)
    manifest={"status":"REFERENCE_CROSS_BUILD_VALIDATION_ONLY",
        "n":len(out),
        "chain_and_ucsc_pass":sum(x["reference_status"]=="CHAIN_UCSC_GRCH38_REF_ALT_PASS" for x in out),
        "ensembl_rsid_pass":sum(x["ensembl_rsid_audit"]=="ENSEMBL_RSID_POSITION_ALLELES_PASS" for x in out),
        "unresolved":sum(x["reference_status"]!="CHAIN_UCSC_GRCH38_REF_ALT_PASS" or x["ensembl_rsid_audit"]!="ENSEMBL_RSID_POSITION_ALLELES_PASS" for x in out),
        "TPMI_GWAS_loaded":False,
        "TPMI_phenotype":"433.21",
        "ALDH2_finemapping":"BLOCKED",
        "ADH1B_finemapping":"EXPLORATORY"}
    (a.out_dir/"TPMI_GRCH38_CROSSBUILD_AUDIT_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("TPMI_CROSSBUILD_AUDIT_COMPLETE",json.dumps(manifest),flush=True)

if __name__=="__main__":main()
