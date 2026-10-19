#!/usr/bin/env python3
"""Generate immutable one-cell full original-646 variant-N sdY experiment.

Clones the prior VERIFIED QTL scalar-N notebook, swaps only R analysis and
output contracts, and asserts that all 646 gene-tissue input sources jointly
contain precisely the already native-BBJ/chain-verified 8,309 GWAS variants.
"""
import argparse,ast,hashlib,json,subprocess,tempfile
from pathlib import Path

def build(source_notebook,source_r,target):
    if target.exists():raise FileExistsError("NOTEBOOK_NO_OVERWRITE")
    note=json.loads(source_notebook.read_text())
    if len(note["cells"])!=1 or note["cells"][0]["cell_type"]!="code":
        raise ValueError("SINGLE_CELL_NOTEBOOK_EXPECTED")
    py="".join(note["cells"][0]["source"])
    r=source_r.read_text()
    if r.count("QTL_VARIANT_N_SDY_STATUS")!=1:
        raise ValueError("SDY_SOURCE_SIGNATURE")
    with tempfile.TemporaryDirectory() as t:
        f=Path(t)/"check.R";f.write_text(r)
        q=subprocess.run(["Rscript","--vanilla","-e",
          "invisible(parse(file=commandArgs(trailingOnly=TRUE)[1]));cat('SDY_SOURCE_PARSE_OK\\n')",
          str(f)],text=True,capture_output=True,timeout=20)
        if q.returncode:raise RuntimeError(q.stderr)
    lines=py.splitlines(keepends=True)
    matches=[i for i,line in enumerate(lines) if line.startswith("R_CODE = ")]
    if len(matches)!=1:raise ValueError("EMBEDDED_R_CODE_MISSING")
    lines[matches[0]]="R_CODE = "+repr(r)+"\n"
    py="".join(lines)
    replacements={
      'is_phase10_11_646_qtl_n_stress':'is_phase10_11_646_variant_n_sdy',
      'IS_PHASE10_11_646_QTL_N_SCALAR_SENSITIVITY_V1':'IS_PHASE10_11_646_VARIANT_N_SDY_CALIBRATION_V1',
      'IS_PHASE10_11_646_QTL_N_R_LOG.txt':'IS_PHASE10_11_646_VARIANT_N_SDY_R_LOG.txt',
      'IS_PHASE10_11_646_qtl_n.R':'IS_PHASE10_11_646_variant_n_sdy.R',
      'QTL_N_SENSITIVITY_R_PARSE_OK':'VARIANT_N_SDY_R_PARSE_OK',
      'IS_646_QTL_N_SCALAR_SENSITIVITY.tsv':'IS_VARIANT_N_SDY_CALIBRATION_STATUS_646.tsv',
      'IS_646_QTL_N_SCALAR_SUMMARY.tsv':'IS_VARIANT_N_SDY_CALIBRATION_SUMMARY.json',
      'IS_646_QTL_N_SOURCE_SHA256.tsv':'IS_646_VARIANT_N_SDY_SOURCE_SHA256.tsv',
      'IS_PHASE10_11_646_QTL_N_SCALAR_COLAB_V1':'IS_PHASE10_11_646_VARIANT_N_SDY_COLAB_V1',
      'FULL_646_SCALAR_N_STRESS_COMPLETE':'FULL_646_VARIANT_N_SDY_CALIBRATION_COMPLETE',
      'N_MIN_MEDIAN_FIRST_SCALAR_STRESS_NOT_PER_SNP_N_MODEL_OR_CAUSAL_GENE':'VARIANT_N_ASSISTED_COMMON_SDY_ONLY_NOT_PER_SNP_SE_CORRECTION',
      'IS_646_N_SCALAR_COLAB_MANIFEST.json':'IS_646_VARIANT_N_SDY_COLAB_MANIFEST.json',
      'IS_646_N_SCALAR_FULL_RUN_LOG.txt':'IS_646_VARIANT_N_SDY_FULL_RUN_LOG.txt',
      'IS_PHASE10_11_646_QTL_N_STRESS=RUN_ENDED_AND_STATUS_ACCOUNTED':'IS_PHASE10_11_646_VARIANT_N_SDY=RUN_ENDED_AND_STATUS_ACCOUNTED',
      'N_SCALAR_646_SUMMARY_FILE=':'VARIANT_N_SDY_646_SUMMARY_FILE=',
    }
    for old,new in replacements.items():
        if py.count(old)!=1:raise ValueError(f"NOTEBOOK_TEMPLATE_SIGNATURE_COUNT_{py.count(old)}:{old}")
        py=py.replace(old,new)
    # Disclose exact statistical model and replace misleading scenario list.
    s='"sensitivity_scenarios":["ORIGINAL_FIRST_EQ_MAX","MEDIAN_SNP_N","MIN_SNP_N"],'
    if py.count(s)!=1:raise ValueError("OLD_SCENARIO_FIELD_MISSING")
    py=py.replace(s,'"variance_calibration_scenarios":["HISTORICAL_COMMON_N_SDY","PER_SNP_N_ASSISTED_COMMON_SDY"],')
    s='"n_assumption_status":"THREE_SCALAR_N_STRESS_NOT_SNP_SPECIFIC_MODEL",'
    if py.count(s)!=1:raise ValueError("OLD_MODEL_STATUS_MISSING")
    py=py.replace(s,'"n_assumption_status":"PER_VARIANT_QTL_N_ONLY_IN_COMMON_SDY_ESTIMATE_NO_VARIANT_SE_RECALC",')
    marker='if shutil.which("Rscript") is None:'
    if py.count(marker)!=1:raise ValueError("SETUP_ANCHOR_MISSING")
    catalog_code="""# Exact independent native BBJ chain source matching across ALL 646 QTL inputs.
EXPECTED_BBJ_GWAS_SNP_CATALOG_SHA256 = "d1e237d3c0ad8281adf3d9fe78f69091505eeeab6d152ae0bc54edd5e964853f"
SNP_CATALOG=WORK/"IS_FOUR_BBJ_LOCI_NATIVE_VERIFIED_GWAS_CATALOG.tsv"
unique_gwas={}
for fn in input_files:
    with fn.open(encoding="utf8",newline="") as h:
        reader=csv.DictReader(h,delimiter="\\t")
        for z in reader:
            k=(z["locus"],z["variant_id"])
            v=(z["match_key"],z["gwas_beta"],z["gwas_eaf"])
            if k in unique_gwas and unique_gwas[k]!=v:
                raise ValueError(f"GWAS_SOURCE_EFFECT_DISAGREES_ACROSS_GTEX_ASSAYS: {k}")
            unique_gwas[k]=v
with SNP_CATALOG.open("w",encoding="utf8",newline="") as h:
    writer=csv.writer(h,delimiter="\\t",lineterminator="\\n")
    writer.writerow(["locus","variant_id","match_key","gwas_beta","gwas_eaf"])
    for (locus,variant),v in sorted(unique_gwas.items()):
        writer.writerow([locus,variant,*v])
print("GWAS_UNIQUE_VARIANTS_ALL_646=",len(unique_gwas))
print("GWAS_CATALOG_SHA256=",sha(SNP_CATALOG))
assert len(unique_gwas)==8309, "UNEXPECTED_UNIQUE_GWAS_VARIANTS"
assert sha(SNP_CATALOG)==EXPECTED_BBJ_GWAS_SNP_CATALOG_SHA256, "NATIVE_BBJ_EFFECT_ALLELE_CATALOG_NOT_SAME_AS_ALL_646"
print("FULL_646_GWAS_NATIVE_ALLELE_CATALOG_VERIFIED_BY_EXACT_SET_AND_EFFECTS")
"""
    py=py.replace(marker,catalog_code+marker)
    s='manifest={'
    if py.count(s)!=1:raise ValueError("MANIFEST_MARKER")
    py=py.replace(s,'shutil.copyfile(SNP_CATALOG,OUTPUT/"IS_FOUR_BBJ_LOCI_NATIVE_VERIFIED_GWAS_CATALOG.tsv")\n'+s)
    s='    "g1_script_sha256":sha(G1_SCRIPT),'
    if py.count(s)!=1:raise ValueError("SOURCE_G1_SIGNATURE")
    py=py.replace(s,s+'\n    "native_BBJ_8309_catalog_sha256":sha(SNP_CATALOG),\n    "native_BBJ_8309_unique_variants_covered_by_all_646":len(unique_gwas),')
    py=py.replace("One-cell Phase10/11 original 646-assay QTL scalar N sensitivity (min/median/first).",
        "One-cell 646 historical BBJ x GTEx variant-specific QTL-N assisted common-sdY calibration.")
    ast.parse(py)
    note["cells"][0]["source"]=py.splitlines(keepends=True)
    note["cells"][0]["execution_count"]=None
    note["cells"][0]["outputs"]=[]
    note["metadata"]["colab"]["name"]=target.name
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(note,indent=1,ensure_ascii=False)+"\n")
    assert "".join(json.loads(target.read_text())["cells"][0]["source"])==py
    print("ONE_CELL_COLAB_CREATED",target.name,
          "R_SHA256",hashlib.sha256(r.encode()).hexdigest(),
          "NOTEBOOK_SHA256",hashlib.sha256(target.read_bytes()).hexdigest(),
          "CODE_LINE_COUNT",len(py.splitlines()))
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--template",required=True,type=Path)
    p.add_argument("--source-r",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args();build(a.template,a.source_r,a.out)
