#!/usr/bin/env python3
"""Reproduce MD5-gated EUR AIS + existing EAS/Japanese gene-universe extension.
No GWAS or canonical results are modified. Output is candidate evidence, not causal claims.
"""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v2_ancestry")
BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
DATA=Path("/srv/is-analysis/data/is/processed/gigastroke/broad_v1/GCST90104540_AIS_EUR_GRCh37.allele_pairs.tsv.gz")
CODE=Path(__file__).resolve().parent
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()
def main(args):
    qc=Path(str(args.eur)+".qc.json")
    if not args.eur.exists() or not qc.exists():
        raise FileNotFoundError("Real source-verified EUR normalized input and QC metadata required")
    metadata=json.loads(qc.read_text())
    if metadata.get("source_md5_verified")!="1bfe4ae8a5042fb24cb0a562b05b0d2f":
        raise ValueError("Source EUR GWAS MD5 verification failed")
    if metadata.get("counts",{}).get("valid",0)<=0:raise ValueError("Empty normalized GWAS")
    cmds=[
      ["build_ancestry_expanded_discovery.py","--eur",str(args.eur),
       "--base",str(args.base),"--out",str(args.out)],
      ["map_all_genes_gencode_v19.py","--root",str(args.out)],
      ["compare_anchor_gene_coverage.py","--root",str(args.out)],
      ["build_ancestry_gene_evidence_v3.py","--out",str(args.out),"--base",str(args.base)]
    ]
    args.out.mkdir(parents=True,exist_ok=True)
    anchor=args.base/"LEGACY_ANCHOR_CANDIDATES.tsv"
    copied=args.out/"LEGACY_ANCHOR_CANDIDATES.tsv"
    if copied.exists() and sha(copied)!=sha(anchor):
        raise ValueError("Existing anchor file differs from canonical baseline")
    if not copied.exists():
        copied.write_bytes(anchor.read_bytes())
    for cmd in cmds:
        p=subprocess.run([sys.executable,str(CODE/cmd[0])]+cmd[1:],
                         capture_output=True,text=True)
        if p.returncode!=0:
            raise RuntimeError("Step failed "+cmd[0]+"\n"+p.stderr+"\n"+p.stdout[-1000:])
    g=json.loads((args.out/"IS_ALL_GENE_UNIVERSE_SUMMARY.json").read_text())
    ev=json.loads((args.out/"IS_ANCESTRY_GENE_EVIDENCE_V3_SUMMARY.json").read_text())
    disc=json.loads((args.out/"IS_ANCESTRY_EXPANSION_DISCOVERY_SUMMARY.json").read_text())
    baseline=json.loads((args.base/"IS_ALL_GENE_UNIVERSE_SUMMARY.json").read_text())
    coverage=json.loads((args.out/"IS_CANDIDATE_COVERAGE_SUMMARY.json").read_text())
    if disc["original_eas_japanese_components"]!=30 or disc["combined_rows"]!=g["regions"]:
        raise ValueError("Unexpected source region dimensions")
    if ev["eas_gene_region_pairs"]!=baseline["region_gene_associations"]:
        raise ValueError("EAS baseline gene-region association regression")
    if ev["gene_region_pairs"]!=g["region_gene_associations"]:
        raise ValueError("Lost ancestry gene-region association")
    if ev["eur_region_pairs_inappropriately_reusing_coloc"]!=0:
        raise ValueError("Cross-population coloc propagation")
    if coverage["legacy_anchor_count"]!=9:
        raise ValueError("Failed to preserve original nine anchors")
    files=["IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv","IS_EUR_AIS_PROVISIONAL_REGIONS.tsv",
           "IS_ALL_GENE_WINDOW_UNIVERSE.tsv","IS_LEGACY_ANCHOR_COVERAGE.tsv",
           "IS_ANCESTRY_GENE_EVIDENCE_V3.tsv"]
    result={"status":"PASS","eur_source_verified":metadata["source_md5_verified"],
      "total_regions":g["regions"],"eur_clusters":disc["eur_distance_clusters"],
      "eur_gws_clusters":disc["eur_gws_clusters"],"unique_gene_ids":g["unique_gene_ids"],
      "gene_region_pairs":ev["gene_region_pairs"],"legacy_anchor_count":9,
      "eur_eas_coloc_leaks":0,"sha256":{f:sha(args.out/f) for f in files},
      "scientific_caveat":"Provisional regions, positional genes only; no independent locus, causal gene or EUR LD evidence."}
    (args.out/"IS_ANCESTRY_EXPANSION_REPRODUCIBILITY.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--base",type=Path,default=BASE)
    ap.add_argument("--eur",type=Path,default=DATA)
    ap.add_argument("--out",type=Path,default=ROOT)
    main(ap.parse_args())
