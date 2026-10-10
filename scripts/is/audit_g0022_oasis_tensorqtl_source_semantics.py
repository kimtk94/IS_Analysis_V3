#!/usr/bin/env python3
"""Bounded original OASIS raw-cis source vs Japan Omics Bokeh crosswalk.

One LINC01409 B_Activated-L2 first-TAR-file source; no G0022 full nominal data.
Confirm web field lineage, not disease causal inference or actual PLINK A1/ALT.
"""
import argparse, csv, gzip, hashlib, io, json, math, struct, urllib.request, zlib
from pathlib import Path
from audit_g0022_oasis_browser_qtl import parse
from audit_g0022_egead1054_nominal_archive import URL, EXPECTED_ZIP_SIZE, RANGE_BYTES, REQUIRED_COLUMNS

BROWSER_URL="https://japan-omics.jp/gene/OASIS-B_Activated-L2?input_value=ENSG00000237491"
AUTHOR_QTL_URL="https://github.com/REdahiro/OASIS_project/blob/main/scripts/cis-eQTL_mapping/tensorQTL_cis_nominal.py"
AUTHOR_COLOC_URL="https://github.com/REdahiro/OASIS_project/blob/main/scripts/colocalization/coloc.R"
TENSORQTL_DOCS="https://github.com/broadinstitute/tensorqtl/blob/master/docs/outputs.md"

def fetch_first_bytes():
    with urllib.request.urlopen(urllib.request.Request(URL,method="HEAD"),timeout=20) as r:
        if int(r.headers.get("Content-Length","0"))!=EXPECTED_ZIP_SIZE or "bytes" not in r.headers.get("Accept-Ranges",""):
            raise ValueError("Source metadata changed")
    req=urllib.request.Request(URL,headers={"Range":f"bytes=0-{RANGE_BYTES-1}"})
    with urllib.request.urlopen(req,timeout=25) as r:
        if r.status!=206 or r.headers.get("Content-Range")!=f"bytes 0-{RANGE_BYTES-1}/{EXPECTED_ZIP_SIZE}":
            raise ValueError("HTTP byte range not honored")
        prefix=r.read(RANGE_BYTES+1)
    if len(prefix)!=RANGE_BYTES:
        raise ValueError("Unexpected remote bytes")
    return prefix

def source_rows_from_prefix(prefix, text_limit=110000):
    if prefix[:4]!=b"PK\x03\x04" or len(prefix)!=RANGE_BYTES:
        raise ValueError("Invalid bounded ZIP prefix")
    n,x=struct.unpack_from("<HH",prefix,26)
    if prefix[30:30+n]!=b"eQTL_summary_statistics.tar":
        raise ValueError("Wrong ZIP member")
    z=prefix[30+n+x:]
    partial_tar=zlib.decompressobj(-15).decompress(z,300000)
    if partial_tar[:100].split(b"\0")[0]!=b"B_Activated_PC15_MAF0.05.cis_nominal.txt.gz":
        raise ValueError("Wrong TAR first gene-cell member")
    inflated=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(partial_tar[512:],text_limit)
    text=inflated.decode("utf8")
    raw=list(csv.DictReader(text.splitlines()[:-1],delimiter="\t"))
    if not raw or list(raw[0])!=list(REQUIRED_COLUMNS):
        raise ValueError("Original columns changed")
    if any(r["phenotype_id"]!="ENSG00000237491" or r["gene"]!="LINC01409" for r in raw):
        raise ValueError("Partial raw gene changed")
    if len(raw)<500:
        raise ValueError("Too few original rows")
    return raw

def compare_source_browser(raw, web, tolerance=0.000501):
    by_key={x["variant_id_hg38"]:x for x in web}
    if len(by_key)!=len(web):
        raise ValueError("Browser has duplicate SNP identities")
    joined=[]
    for a in raw:
        key=a["variant_id"].replace("_",":")
        if key not in by_key:
            raise ValueError("Original tested variant absent in browser: "+key)
        b=by_key[key]
        diffs={
            "p":abs(float(a["pval_nominal"])-float(b["pval_nominal"])),
            "beta":abs(float(a["slope"])-float(b["effect_size"])),
            "se":abs(float(a["slope_se"])-float(b["effect_size_SE"])),
            "af":abs(float(a["af"])-float(b["maf"])),
        }
        if any(not math.isfinite(x) or x>tolerance for x in diffs.values()):
            raise ValueError("Raw Bokeh rounding mismatch "+key+" "+str(diffs))
        joined.append({
            "variant_id":key,
            "source_p":a["pval_nominal"],"browser_p":b["pval_nominal"],
            "source_slope":a["slope"],"browser_effect_size":b["effect_size"],
            "source_slope_se":a["slope_se"],"browser_effect_size_SE":b["effect_size_SE"],
            "source_af":a["af"],"browser_labeled_maf":b["maf"],
            "max_abs_rounding_difference":max(diffs.values()),
        })
    if not any(float(r["source_af"])>.5 and float(r["browser_labeled_maf"])>.5 for r in joined):
        raise ValueError("High AF evidence missing")
    return joined

def execute(outdir,prefix,page):
    raw=source_rows_from_prefix(prefix)
    web=parse(page.decode("utf8","replace"),"LINC01409")
    rows=compare_source_browser(raw,web)
    outdir.mkdir(exist_ok=True,parents=True)
    detail=outdir/"OASIS_ORIGINAL_TO_BROWSER_CIS_NOMINAL_NUMERIC_CROSSWALK.tsv"
    with detail.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
    audit={
        "audit":"OASIS_EGEAD1054_ORIGINAL_FIRST_CELL_RAW_BROWSER_CROSSWALK_V1",
        "origin":"E-GEAD-1054 B_Activated_PC15_MAF0.05.cis_nominal.txt.gz",
        "origin_gene":"LINC01409",
        "original_raw_rows_compared":len(raw),
        "browser_total_variant_rows":len(web),
        "row_by_variant_match":len(rows),
        "raw_AF_greater_than_half":sum(float(r["source_af"])>.5 for r in rows),
        "browser_labeled_MAF_greater_than_half":sum(float(r["browser_labeled_maf"])>.5 for r in rows),
        "max_abs_column_difference":max(float(r["max_abs_rounding_difference"]) for r in rows),
        "browser_field_p":"RAW_pval_nominal_ROUNDED",
        "browser_field_effect_size":"RAW_tensorQTL_slope_ROUNDED",
        "browser_field_effect_size_SE":"RAW_tensorQTL_slope_se_ROUNDED",
        "browser_field_maf":"RAW_tensorQTL_af_ROUNDED_NOT_MINOR_FREQUENCY",
        "original_prefix_sha256":hashlib.sha256(prefix).hexdigest(),
        "browser_html_sha256":hashlib.sha256(page).hexdigest(),
        "original_5cell_G0022_source_data_directly_checked":False,
        "source_genotype_PLINK_A1_equals_variant_ALT_externally_verified":False,
        "G0022_full_cis_coloc_ready":False,
        "source_author_mapping_code":AUTHOR_QTL_URL,
        "source_author_coloc_script":AUTHOR_COLOC_URL,
        "tensorQTL_documentation":TENSORQTL_DOCS,
        "evidence_gate":"STRONG_ORIGINAL_BROWSER_COLUMN_LINEAGE_ON_CHR1_CELL_ONLY"
    }
    (outdir/"OASIS_ORIGINAL_BROWSER_SEMANTICS_AUDIT.json").write_text(json.dumps(audit,indent=2)+"\n")
    return audit

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir",type=Path,required=True)
    ap.add_argument("--online",action="store_true")
    ap.add_argument("--source-prefix",type=Path)
    ap.add_argument("--browser-html",type=Path)
    a=ap.parse_args()
    a.outdir.mkdir(exist_ok=True,parents=True)
    if a.online:
        prefix_path=a.outdir/"DDBJ_ORIGINAL_ZIP_FIRST_262144_BYTES.bin"
        if prefix_path.exists():
            prefix=prefix_path.read_bytes()
        else:
            prefix=fetch_first_bytes()
            prefix_path.write_bytes(prefix)
        page_path=a.outdir/"OASIS_LINC01409_ACTIVATED_B_L2_ORIGINAL.html"
        if page_path.exists():
            page=page_path.read_bytes()
        else:
            req=urllib.request.Request(BROWSER_URL,headers={"User-Agent":"Mozilla/5.0"})
            with urllib.request.urlopen(req,timeout=20) as resp:
                page=resp.read(1200001)
            if len(page)>1200000:
                raise ValueError("Browser page too large")
            page_path.write_bytes(page)
    else:
        if not a.source_prefix or not a.browser_html:
            ap.error("Offline requires --source-prefix and --browser-html")
        prefix=a.source_prefix.read_bytes()
        page=a.browser_html.read_bytes()
    print(json.dumps(execute(a.outdir,prefix,page),indent=2))

if __name__=="__main__":
    main()
