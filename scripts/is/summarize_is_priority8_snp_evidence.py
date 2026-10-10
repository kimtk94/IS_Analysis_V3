#!/usr/bin/env python3
"""Priority-eight direct SNP coloc audit; never modifies original inputs."""
import csv, json, hashlib, math
from pathlib import Path
import argparse
GENES=("FGF5","CALHM2","NEURL1","C4orf22","INA","SH3PXD2A","COL4A2","COL4A1")
PRIORS=(1e-6,3e-6,1e-5,3e-5,1e-4)
def get(path):
    with path.open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))
def run(replay,input_qc,maf,master,output):
    if output.exists():raise FileExistsError("NO_OVERWRITE")
    a,b,c,h=map(get,(replay,input_qc,maf,master))
    if tuple(map(len,(a,b,c,h)))!=(120,8,24,646):raise ValueError("DENOMINATOR_MISMATCH")
    audit=[]
    for gene in GENES:
        x=[r for r in a if r["gene"]==gene]
        y=[r for r in b if r["gene"]==gene]
        z=[r for r in c if r["gene"]==gene]
        if len(x)!=15 or len(y)!=1 or len(z)!=3:raise ValueError("MISSING_GENE_"+gene)
        q=y[0]
        if {(r["qtl_n_mode"],float(r["p12"])) for r in x}!={(n,p) for n in ("FIRST","MIN","MAX") for p in PRIORS}:
            raise ValueError("INVALID_CONDITIONS")
        f=lambda p:next(r for r in x if r["qtl_n_mode"]=="FIRST" and float(r["p12"])==p)
        base=f(1e-5)
        minimum=next(r for r in x if r["qtl_n_mode"]=="MIN" and float(r["p12"])==1e-5)
        orig=[r for r in h if r["locus"]==q["locus"] and r["dataset_key"]==q["tissue"] and r["gene_symbol"]==gene]
        if len(orig)!=1 or int(orig[0]["nsnps"])!=int(q["n_snps"]):raise ValueError("ORIGINAL_MISMATCH")
        if abs(float(orig[0]["qtl_n"])-float(q["qtl_n_first"]))>1e-8:raise ValueError("QTL_N_MISMATCH")
        for k in range(5):
            if abs(float(orig[0][f"PP.H{k}"])-float(base[f"PP_H{k}"]))>1e-8:raise ValueError("HYPOTHESIS_NOT_REPLAYED")
        for r in x:
            p=[float(r[f"PP_H{k}"]) for k in range(5)]
            if not all(math.isfinite(v) and 0<=v<=1 for v in p) or abs(sum(p)-1)>1e-8:
                raise ValueError("H0_H4_INVALID")
        base_maf=next(r for r in z if r["scenario"]=="UNFILTERED_HISTORICAL")
        strict_maf=next(r for r in z if r["maf_filter_max_diff"] not in ("NA","") and float(r["maf_filter_max_diff"])==0.1)
        if abs(float(base_maf["PP_H4"])-float(base["PP_H4"]))>1e-8:raise ValueError("MAF_REFERENCE_MISMATCH")
        audit.append({
          "gene":gene,"locus":q["locus"],"dataset_key":q["tissue"],"nsnps":int(q["n_snps"]),
          "qtl_n_first":int(float(q["qtl_n_first"])),
          "qtl_n_min":int(float(q["qtl_n_min"])),
          "qtl_n_max":int(float(q["qtl_n_max"])),
          "pp_h3_baseline":float(base["PP_H3"]),
          "pp_h4_baseline":float(base["PP_H4"]),
          "pp_h4_p12_1e6":float(f(1e-6)["PP_H4"]),
          "pp_h4_p12_1e4":float(f(1e-4)["PP_H4"]),
          "pp_h4_min_qtl_n":float(minimum["PP_H4"]),
          "excluded_maf_delta_gt_0_1":int(strict_maf["n_excluded"]),
          "lead_eqtl_kept_maf_delta_le_0_1":strict_maf["lead_eqtl_retained"]=="TRUE",
          "posthoc_filtered_pp_h4":float(strict_maf["PP_H4"]),
          "status":"ORIGINAL_FIVE_HYPOTHESES_SNP_REPLAY_PASS_NO_CAUSAL_CLAIM"
        })
    source_inputs={}
    for row in b:
        file=Path(row["source_file"])
        if not file.is_file() or file.stat().st_size<1000:
            raise ValueError("RAW_SOURCE_FILE_LOST_"+row["gene"])
        source_inputs[row["gene"]]={
            "source_file":str(file),
            "bytes":file.stat().st_size,
            "sha256":hashlib.sha256(file.read_bytes()).hexdigest()
        }
    if len(source_inputs)!=8:
        raise ValueError("EIGHT_SOURCE_HASHES_REQUIRED")
    source={"replay":replay,"source_input_qc":input_qc,"maf_qc":maf,"original_646_master":master}
    obj={"schema":"IS_PRIORITY8_SNP_ABF_20261010_V1",
      "status":"8_SELECTED_GENE_TISSUE_PAIRS_SNP_LEVEL_REPLAYED",
      "original_legacy_tests":646,"verified_replays":8,"not_snp_replayed":638,
      "runs":120,"source_p1":1e-4,"source_p2":1e-4,"source_p12":1e-5,
      "case_control_n":174686,"case_count":22664,
      "LD_ancestry_verified":False,"causal_gene_mechanism_established":False,
      "note":"MAF difference exclusion is posthoc and does not establish harmonization. Min/Max scalar sample-N sensitivity is not per-SNP missingness modeling.",
      "hashes":{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in source.items()},
      "per_gene_source_SNP_file_provenance":source_inputs,
      "genes":audit}
    output.mkdir(parents=True)
    (output/"IS_PRIORITY8_SOURCE_VERIFIED.json").write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n")
    with (output/"IS_PRIORITY8_SOURCE_VERIFIED.tsv").open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(audit[0]),delimiter="\t")
        w.writeheader();w.writerows(audit)
    lines=["# Direct SNP ABF replay: priority eight source-matched gene–tissue pairs","",
       "8/646 original selected historical GTEx gene–tissue tests numerically reproduced with original five H0–H4 hypotheses.",
       "120 runs = 8 pairs x 5 p12 priors x 3 scalar QTL N assumptions. Other 638/646 remain unreplayed.",
       "No novel gene-level causal mechanism is proven. QTL/GWAS ancestry difference, LD, tissue selection and multi-signal issues remain.",
       "","| Gene | Original H4 | p12=1e-6 | p12=1e-4 | MIN QTL N H4 | SNPs | Removed by MAF difference 0.1 | Lead eQTL kept? |",
       "|---|---:|---:|---:|---:|---:|---:|---|"]
    for v in audit:
        lines.append(f'| {v["gene"]} | {v["pp_h4_baseline"]:.3f} | {v["pp_h4_p12_1e6"]:.3f} | {v["pp_h4_p12_1e4"]:.3f} | {v["pp_h4_min_qtl_n"]:.3f} | {v["nsnps"]} | {v["excluded_maf_delta_gt_0_1"]} | {"yes" if v["lead_eqtl_kept_maf_delta_le_0_1"] else "**NO**"} |')
    lines +=["","MAF-based exclusion is diagnostic only and cannot be used to increase apparent evidence. The full SNP inputs remain separate from Git.","",
       "Source provenance and input SHA256 checksums are in IS_PRIORITY8_SOURCE_VERIFIED.json."]
    (output/"IS_PRIORITY8_RESULT_REPORT.md").write_text("\n".join(lines)+"\n")
    print("PRIORITY8_SOURCE_VERIFIED 8/646 PASS; 638 SNP inputs not independently replayed")
    return obj
if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("replay","input_qc","maf","master","output"):
        p.add_argument("--"+name,required=True,type=Path)
    x=p.parse_args();run(x.replay,x.input_qc,x.maf,x.master,x.output)
