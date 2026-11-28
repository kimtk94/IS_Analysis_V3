#!/usr/bin/env python3
"""Re-evaluate the complete archived 36-assay/20-gene chr10 family against published results.
Four historical anchor genes and CALHM2/NEURL1 context are distinct from the fixed assay universe.
Exact ID or frozen symbol/explicit aliases only; no causal ranking or statistical promotion.
"""
import argparse,collections,csv,hashlib,json,re
from pathlib import Path
from openpyxl import load_workbook
ROOT=Path(__file__).resolve().parents[2]
EVID=ROOT/"docs/is/evidence/verified_p3a_20261124"
COVERAGE_SHA="6fc64e2068b44e45a1bf16e1162e20062e4913cca9fc1d3b2e9b702cd118e344"
MODEL_SHA="023750fde4a2a57c1604793dfa2b78ba65a29dbe06ff849e17382a559906eb2a"
PUBLISHER_SHA="7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4"
ANCHORS={"FGF5","ALDH2","SH3PXD2A","COL4A2"}
def require(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def load_verified(p,expected):
    require(sha(p)==expected,"Pinned evidence digest mismatch: "+p.name)
    return json.loads(p.read_text())
def candidates():
    cov=load_verified(EVID/"IS_P3A_36x2_original_cis_allele_coverage_HOLD.json",COVERAGE_SHA)
    mod=load_verified(EVID/"IS_P3A_36_original_QTL_SE_MAF_AN_residual_proxy_N_model_HOLD.json",MODEL_SHA)
    model={}
    for m in mod["groups"]:
        k=(m["gene_id"],m["QTL_dataset"])
        require(k not in model and not m["formal_coloc_abf_eligible"] and not m["formal_coloc_susie_eligible"],"Archived gate drift")
        model[k]=m
    pairs=collections.defaultdict(dict)
    for r in cov["rows"]:
        k=(r["gene_id"],r["QTL_dataset"]);a=r["ancestry"]
        require(k in model and a in {"EAS","EUR"} and a not in pairs[k],"Natural key error")
        require(r["gene_symbol"]==model[k]["gene_symbol"] and not r["formal_coloc_abf_eligible"] and not r["formal_coloc_susie_eligible"],"Scientific gate drift")
        pairs[k][a]=r
    require(len(model)==len(pairs)==36 and len(cov["rows"])==72 and set(model)==set(pairs),"Fixed assay denominator")
    out={}
    for (gid,ds),v in sorted(pairs.items()):
        require(set(v)=={"EAS","EUR"},"Missing ancestry")
        g=v["EAS"]["gene_symbol"]
        require(g==v["EUR"]["gene_symbol"],"Gene symbol drift")
        d=out.setdefault(gid,{"gene_id":gid,"archived_symbol":g,"assays":[]})
        require(d["archived_symbol"]==g,"Gene mapping drift")
        d["assays"].append({"dataset":ds,"EAS_strict_fraction":v["EAS"]["fraction_QTL_strict_nonpal_with_GWAS_match"],"EUR_strict_fraction":v["EUR"]["fraction_QTL_strict_nonpal_with_GWAS_match"],"EAS_all_cis_fraction":v["EAS"]["fraction_QTL_ALL_original_cis_with_nonpal_GWAS_match"],"EUR_all_cis_fraction":v["EUR"]["fraction_QTL_ALL_original_cis_with_nonpal_GWAS_match"],"EAS_matched_rows":v["EAS"]["matched_nonpal_SNP_assay_records"],"EUR_matched_rows":v["EUR"]["matched_nonpal_SNP_assay_records"],"regression_N_certified":False,"sdY_certified":False,"status":"HOLD"})
    require(len(out)==20 and len({v["archived_symbol"] for v in out.values()})==20,"Fixed gene denominator")
    return out
def tokens(v):
    if not isinstance(v,str):return []
    # Preserve punctuation in historical RP11 symbols. Never match substrings.
    return [s.strip("()[]") for s in re.split(r"[\s,;|/]+",v.strip()) if s.strip("()[]")]
def published(path,genes):
    require(path.stat().st_size==756910 and sha(path)==PUBLISHER_SHA,"Pinned GIGASTROKE publisher workbook differs")
    symbols={v["archived_symbol"]:gid for gid,v in genes.items()}
    # NEURL is the historical symbol for the exact ENSG00000107954 ID; no guessed renaming.
    require(symbols.get("NEURL")=="ENSG00000107954","Historical NEURL mapping drift")
    symbols["NEURL1"]="ENSG00000107954"
    # Known symbol aliases of these stable IDs; match separately and report alias provenance.
    alias={"BORCS7":"ENSG00000166275","ATP5MD":"ENSG00000173915"}
    for a,gid in alias.items():require(gid in genes,"Alias ID not in frozen universe")
    symbols.update(alias)
    for g in ANCHORS:symbols[g]="anchor:"+g
    selected={};matches=collections.defaultdict(list)
    wb=load_workbook(path,read_only=True,data_only=True,keep_links=False)
    try:
        for s in wb.worksheets:
            m=re.fullmatch(r"Suppl[.]?\s*Table[.]?\s*(6|7|27|30|38)(?:[a-zA-Z])?",s.title,re.I)
            if not m:continue
            require((s.max_row or 0)<=100000 and (s.max_column or 0)<=300,"Workbook size bound")
            n=int(m[1]);head=[];count=0
            for i,row in enumerate(s.iter_rows(values_only=True),1):
                count=i
                vals=["" if v is None else str(v) for v in row]
                if i<=8:head.append(vals[:50])
                found=collections.defaultdict(list)
                for j,v in enumerate(row,1):
                    for t in tokens(v):
                        stable=t.split(".")[0]
                        if stable in genes:found[stable].append({"column":j,"token":t,"mode":"stable_ensembl_id"})
                        elif t in symbols:found[symbols[t]].append({"column":j,"token":t,"mode":"explicit_alias" if t in alias or t=="NEURL1" else "exact_symbol"})
                for k,fields in found.items():
                    matches[k].append({"sheet":s.title,"table":n,"row":i,"matched_fields":fields,"source_values":vals[:60]})
            selected[s.title]={"table":n,"rows_scanned":count,"first_rows":head}
    finally:wb.close()
    require({6,7,27,30,38}<=set(v["table"] for v in selected.values()),"Expected publication tables not present")
    return selected,matches
def main():
    p=argparse.ArgumentParser();p.add_argument("--workbook",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    g=candidates();s,m=published(a.workbook,g)
    rows=[]
    for gid,v in sorted(g.items(),key=lambda x:x[1]["archived_symbol"]):
        q=v["assays"];mch=m.get(gid,[])
        rows.append({**v,"gene":v["archived_symbol"],"EAS_80pct_assays":sum(x["EAS_strict_fraction"]>=.8 for x in q),"EUR_80pct_assays":sum(x["EUR_strict_fraction"]>=.8 for x in q),"both_80pct_assays":sum(x["EAS_strict_fraction"]>=.8 and x["EUR_strict_fraction"]>=.8 for x in q),"published_exact_token_rows_by_table":dict(collections.Counter(str(x["table"]) for x in mch)),"published_rows":mch,"scientific_status":"HOLD_NO_NEW_COLOC","novelty_status":"REVIEW_PUBLISHED_ROW_CONTEXT" if mch else "NO_MATCH_IN_SELECTED_TABLES_NOT_NOVELTY"})
    res={"status":"COMPLETE_ARCHIVED_FAMILY_PUBLICATION_AUDIT","archived_gene_count":len(rows),"archived_assays":sum(len(r["assays"]) for r in rows),"ancestry_cells":72,"independent_archived_GWAS_loci":1,"EAS_80pct_assays":sum(r["EAS_80pct_assays"] for r in rows),"EUR_80pct_assays":sum(r["EUR_80pct_assays"] for r in rows),"both_80pct_assays":sum(r["both_80pct_assays"] for r in rows),"formal_coloc_eligible":0,"genes":rows,"anchors_not_in_36_assay_family":[{"gene":s,"published_rows":m.get("anchor:"+s,[])} for s in sorted(ANCHORS)],"publisher_sheets":s,"source_sha256":{"coverage":COVERAGE_SHA,"model":MODEL_SHA,"gigastroke_workbook":PUBLISHER_SHA},"limitations":["The 20 genes/36 groups belong to one chr10 region; the broad positional discovery universe of 2225 genes is not audited here. Four historical anchors are separate comparator contexts.","Exact token occurrence does not prove a reported positive association; read table header/row/phenotype and study. Aliases must be reviewed against stable IDs. Missing matches are not evidence of novelty.","80% is a prior technical screen, not a universal scientific threshold. Coverage differs by ancestry and is based on strict nonpal SNPs; full-cis coverage is lower.","No independent validation, no new H4, no causal ranking, no new MR. All 36 groups retain scientific HOLD."]}
    a.out.mkdir(parents=True,exist_ok=True);(a.out/"candidate-audit.json").write_text(json.dumps(res,indent=2,allow_nan=False)+"\n")
    with (a.out/"candidate-summary.tsv").open("w",newline="") as f:
        w=csv.writer(f,delimiter="\t");w.writerow(["gene_id","archived_symbol","assays","EAS_80pct","EUR_80pct","both_80pct","table_6_rows","table_7_rows","table_27_rows","table_30_rows","table_38_rows","status"])
        for d in rows:w.writerow([d["gene_id"],d["archived_symbol"],len(d["assays"]),d["EAS_80pct_assays"],d["EUR_80pct_assays"],d["both_80pct_assays"],*[d["published_exact_token_rows_by_table"].get(str(i),0) for i in (6,7,27,30,38)],d["scientific_status"]])
    print(json.dumps({**{k:v for k,v in res.items() if k not in ("genes","anchors_not_in_36_assay_family","publisher_sheets")},"genes":[{k:v for k,v in r.items() if k not in ("assays","published_rows")} for r in rows],"anchor_counts":{x["gene"]:dict(collections.Counter(str(y["table"]) for y in x["published_rows"])) for x in res["anchors_not_in_36_assay_family"]},"sheet_headers":{k:v["first_rows"][:4] for k,v in s.items()}},indent=2,allow_nan=False))
if __name__=="__main__":main()
