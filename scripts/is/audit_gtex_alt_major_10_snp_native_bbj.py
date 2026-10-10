#!/usr/bin/env python3
"""Follow up 10 ALT-major GTEx sites absent from prior BBJ/GTEx coloc cache.

With already saved QTD subset and original BBJ autosomal ZIP, inverse-map
GRCh38 QTL SNP coordinates via LOCAL hg19ToHg38 UCSC chain file, then check
whether BBJ native original GWAS contains the GRCh37 SNP and allele pair.

A missing SNP in coloc is NOT proof of a pipeline bug: there may be GWAS
absence or original analysis exclusions. Palindromic SNPs require special care.
No new databases, no input overwrite, no automatic allele flip.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,zipfile
from pathlib import Path
from collections import defaultdict,Counter

def read(path):
    with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""):h.update(b)
    return h.hexdigest()
def inverse_blocks(chain,chrom):
    out=[]
    current=None
    with gzip.open(chain,"rt") as f:
        for raw in f:
            line=raw.strip()
            if not line:
                current=None
                continue
            if line.startswith("chain "):
                a=line.split()
                if len(a)<13:raise ValueError("Bad chain format")
                score=int(a[1]);tchr=a[2];tstrand=a[4];tpos=int(a[5])
                qchr=a[7];qstrand=a[9];qpos=int(a[10])
                current={
                    "accept":tchr==chrom and qchr==chrom and
                             tstrand=="+" and qstrand=="+",
                    "score":score,"t":tpos,"q":qpos,"id":a[12],
                }
                continue
            if current is None:continue
            a=line.split()
            if len(a) not in (1,3):raise ValueError("Bad chain block format")
            sz=int(a[0])
            if current["accept"]:
                out.append((current["q"],current["q"]+sz,current["t"],
                            current["score"],current["id"]))
            if len(a)==3:
                current["t"]+=sz+int(a[1])
                current["q"]+=sz+int(a[2])
            else:
                current=None
    return out

def invert_one(blocks,pos0):
    candidates=[(t+(pos0-q0),score,chainid)
                for q0,qend,t,score,chainid in blocks if q0<=pos0<qend]
    if not candidates:return {"status":"CHAIN_NO_MAPPING","pos37":None,"candidate_count":0}
    top=max(score for _,score,_ in candidates)
    best={p for p,score,_ in candidates if score==top}
    if len(best)!=1:return {"status":"CHAIN_AMBIGUOUS_TOP_SCORE",
                            "pos37":None,"candidate_count":len(candidates)}
    return {"status":"TOP_CHAIN_SELECTED_UNIQUE_POS","pos37":next(iter(best))+1,
            "candidate_count":len(candidates),"top_score":top}

def scan_native(zip_path,target_pos):
    observed=defaultdict(list);n=0
    with zipfile.ZipFile(zip_path) as z:
        members=[x for x in z.namelist() if x.endswith(".auto.txt.gz")]
        if len(members)!=1:raise ValueError("Autosomal native GWAS member nonunique")
        with z.open(members[0]) as binfile:
            with gzip.open(binfile,"rt",newline="") as text:
                reader=csv.DictReader(text,delimiter="\t")
                for r in reader:
                    n+=1
                    if r["CHR"]=="12" and r["POS"] in target_pos:
                        observed[r["POS"]].append({
                         "v":r["v"],"allele1":r["Allele1"],
                         "allele2":r["Allele2"],"AF_Allele2":r["AF_Allele2"],
                         "beta":r["BETA"],"p":r["p.value"]
                        })
                    if n%5000000==0:
                        print("BBJ_NATIVE_SCAN",n,flush=True)
    return observed,n

def audit(qtd_rows,blocks,native_zip):
    sites=defaultdict(list)
    for r in qtd_rows:
        sites[r["site_GRCh38"]].append(r)
    if len(sites)!=10:
        raise ValueError(f"Expected 10 distinct ALT-major sites, found {len(sites)}")
    base=[]
    for site,rs in sorted(sites.items()):
        chrom,pos,ref,alt=site.split(":")
        if chrom!="12":raise ValueError("Expected L003 chr12")
        inv=invert_one(blocks,int(pos)-1)
        base.append({"site":site,"ref":ref,"alt":alt,
                     "QTL_gene_tissue_records":len(rs),
                     "raw_QTD_ALT_AF_min":min(float(r["QTD_ALT_ac_an"]) for r in rs),
                     "raw_QTD_ALT_AF_max":max(float(r["QTD_ALT_ac_an"]) for r in rs),
                     "palindromic":(ref,alt) in {("A","T"),("T","A"),("C","G"),("G","C")},
                     "pos37":inv.get("pos37"),"chain_status":inv["status"],
                     "chain_candidate_count":inv["candidate_count"]})
    target={str(x["pos37"]) for x in base if x["pos37"] is not None}
    native,nrows=scan_native(native_zip,target)
    for x in base:
        if x["pos37"] is None:
            x.update({"native_status":"UNMAPPED","source_gwas_allele":"UNKNOWN",
                     "native_allele2_AF":"","native_p":"","native_variant_candidates":""})
            continue
        matches=native.get(str(x["pos37"]),[])
        exact=[r for r in matches if
               r["allele1"]==x["ref"] and r["allele2"]==x["alt"]]
        reverse=[r for r in matches if
               r["allele1"]==x["alt"] and r["allele2"]==x["ref"]]
        if len(exact)==1:
            x.update({"native_status":"BBJ_ORIGINAL_SAME_REF_ALT_PAIR",
                "source_gwas_allele":"ALT","native_allele2_AF":exact[0]["AF_Allele2"],
                "native_p":exact[0]["p"]})
        elif len(reverse)==1 and not exact:
            x.update({"native_status":"BBJ_ORIGINAL_ALLELES_REVERSED_REVIEW",
                "source_gwas_allele":"NEEDS_HARMONIZATION",
                "native_allele2_AF":reverse[0]["AF_Allele2"],
                "native_p":reverse[0]["p"]})
        elif not matches:
            x.update({"native_status":"NOT_REPORTED_AT_GRCH37_POSITION",
                "source_gwas_allele":"NOT_REPORTED","native_allele2_AF":"",
                "native_p":""})
        else:
            x.update({"native_status":"GWAS_ALLELE_PAIR_UNRESOLVED",
                "source_gwas_allele":"UNKNOWN","native_allele2_AF":"",
                "native_p":""})
        x["native_variant_candidates"]=";".join(r["v"] for r in matches)
        x["scientific_interpretation"]="DESCRIPTIVE_COVERAGE_CHECK_NOT_INFERRED_ALLELE_ERROR"
    return base,nrows

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--qtd-alt-major-rows",required=True,type=Path)
    parser.add_argument("--chain",required=True,type=Path)
    parser.add_argument("--native-bbj-zip",required=True,type=Path)
    parser.add_argument("--out-dir",required=True,type=Path)
    a=parser.parse_args()
    if a.out_dir.exists() and any(a.out_dir.iterdir()):
        raise FileExistsError("Do not overwrite prior audit")
    rows=read(a.qtd_alt_major_rows)
    blocks=inverse_blocks(a.chain,"chr12")
    if not blocks:raise ValueError("No chr12 chain segments")
    audit_rows,whole_gwas_n=audit(rows,blocks,a.native_bbj_zip)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"IS_10_GTEX_ALT_MAJOR_SNP_NATIVE_BBJ_AVAILABILITY.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(audit_rows[0]),delimiter="\t")
        w.writeheader();w.writerows(audit_rows)
    summary={
       "status":"SOURCE_AVAILABILITY_CHECK_ONLY_NOT_PIPELINE_CAUSATION",
       "QTD_ALT_major_variant_sites":len(audit_rows),
       "source_BBJ_autosomal_rows_scanned":whole_gwas_n,
       "native_gwas_status":dict(Counter(x["native_status"] for x in audit_rows)),
       "chain_status":dict(Counter(x["chain_status"] for x in audit_rows)),
       "input_SHA256":{
        "original_QTD_ALT_major_row_ledger":digest(a.qtd_alt_major_rows),
        "UCSC_chain":digest(a.chain),
       },
       "cause_of_absence_in_646_original_coloc":"PENDING_EVIDENCE",
       "BBJ_GWAS_effect_allele_ALT_confirmed_when_same_refalt":True,
       "GTEx_regression_effect_allele_independently_attested":False,
       "no_allele_flipping_performed":True,"no_causal_gene_promotion":True
    }
    (a.out_dir/"IS_10_GTEX_BBJ_COVERAGE_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("IS_QTD_10_SNP_BBJ_NATIVE_COVERAGE_COMPLETE",json.dumps(summary["native_gwas_status"]),flush=True)
    for r in audit_rows:print("SNP",r["site"],r["chain_status"],r["native_status"],r.get("native_variant_candidates",""),flush=True)
if __name__=="__main__":
    main()
