#!/usr/bin/env python3
"""Read-only public cis-eQTL actual-access gates (no full summary downloads).

Stage 0 only: MAEEA Zenodo attached-file inventory, NBDC Japanese full archive
HEAD, eQTL Catalogue EUR tissue/cell all-association FTP HEAD, GTEx allpairs
requester-pays diagnostic, and JCTF filtered-pairs readme.
Results never assert colocalization readiness based on HTTP or title alone.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import urllib.error
from urllib.parse import urlparse
from urllib.request import Request,urlopen

URLS={
    "MAEEA_ZENODO":"https://zenodo.org/api/records/21296030",
    "JAPAN_HEALTHY_IMMUNE_EQTL_ARCHIVE":"https://humandbs-production.ddbj.nig.ac.jp/files/hum0099/hum0099.v1.eqtl.v1.zip",
    "JAPAN_JCTF_EQTL_README":"https://humandbs-production.ddbj.nig.ac.jp/files/hum0343/readme_wang_qs_et_al_2024_japan_covid19_taskforce.txt",
    "EUR_BLUEPRINT_MONOCYTE_ALL_ASSOC":"https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD000021/QTD000021.all.tsv.gz",
    "GTEX_V8_WHOLE_BLOOD_ALLPAIRS":"https://storage.googleapis.com/gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_eQTL_all_associations/Whole_Blood.allpairs.txt.gz",
}
MAX_META=350000


def query_head(url):
    if urlparse(url).scheme!="https":
        raise ValueError("Only HTTPS sources")
    req=Request(url,method="HEAD",headers={
        "User-Agent":"IS_Analysis_V3_QTL_PublicAccessGate/1.0"
    })
    status={"url":url,"access":"SOURCE_UNAVAILABLE","full_variant_eQTL_QC":"HOLD",
            "usable_for_coloc":"NOT_YET_VERIFIED"}
    try:
        with urlopen(req,timeout=35) as r:
            size=r.headers.get("Content-Length","")
            status["http_status"]=r.status
            status["archive_bytes_advertised"]=int(size) if size.isdigit() else None
            status["final_hostname"]=urlparse(r.geturl()).hostname
            status["access"]="REMOTE_OBJECT_EXISTS_NOT_CONTENT_QC" if r.status==200 else "HTTP_REVIEW"
    except urllib.error.HTTPError as e:
        status["http_status"]=e.code
        status["access"]="POSSIBLY_REQUESTER_PAYS_OR_RESTRICTED" if e.code in (400,401,403) else "HTTP_ERROR"
    except (OSError,TimeoutError,ValueError) as e:
        status["error_short"]=f"{type(e).__name__}: {str(e)[:100]}"
    return status


def limited_get(url,limit=MAX_META):
    req=Request(url,headers={"User-Agent":"IS_Analysis_V3_QTL_PublicAccessGate/1.0"})
    with urlopen(req,timeout=35) as resp:
        b=resp.read(limit+1)
        if len(b)>limit:
            raise ValueError("Metadata exceeds bounded cap")
        return b.decode("utf-8",errors="replace")


def maeea():
    url=URLS["MAEEA_ZENODO"]
    result={"url":url,"access":"SOURCE_UNAVAILABLE","all_cis_variant_level_data":"HOLD",
            "EAS_eQTL_coloc_ready":False}
    try:
        d=json.loads(limited_get(url))
        result["record_id"]=d.get("id")
        files=d.get("files",[])
        result["file_count"]=len(files)
        result["advertised_total_bytes"]=sum(int(x.get("size",0)) for x in files)
        result["file_names"]=[x.get("key") for x in files]
        result["access"]="RECORD_METADATA_ONLY_NO_FILES" if not files else "PUBLIC_FILE_MANIFEST_NOT_DATA_QC"
    except (OSError,ValueError,TimeoutError,urllib.error.URLError,
            json.JSONDecodeError) as e:
        result["error_short"]=f"{type(e).__name__}: {str(e)[:140]}"
    return result


def jctf_readme():
    result={"url":URLS["JAPAN_JCTF_EQTL_README"],"status":"SOURCE_UNAVAILABLE",
            "eQTL_input_full_cis_distribution":"HOLD"}
    try:
        content=limited_get(result["url"],8000).lower()
        has_se="slope_se" in content
        has_hg19="variant_id_hg19" in content
        filtered="p<0.05 or pip>0.001" in content.replace(" ","")
        # Original README has 'p<0.05 or PIP>0.001', normalized strip spaces
        condensed="".join(content.split())
        filtered=("p<0.05orpip>0.001" in condensed)
        result.update({"readme_has_se":has_se,"readme_has_hg19":has_hg19,
                       "filtered_pairs_only":filtered})
        result["status"]="PUBLISHED_SIGNIFICANCE_FILTERED_QTL" if has_se and has_hg19 and filtered else "READ_ME_UNEXPECTED_REVIEW"
    except (OSError,ValueError,TimeoutError,urllib.error.URLError) as e:
        result["error_short"]=f"{type(e).__name__}: {str(e)[:140]}"
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",default="artifacts/is_eQTL_actual_access_gates")
    a=p.parse_args()
    out=Path(a.out_dir)
    out.mkdir(parents=True,exist_ok=True)
    tests={
        "MAEEA_CBEA_EAS":maeea(),
        "JAPAN_HEALTHY_NHA000051":query_head(URLS["JAPAN_HEALTHY_IMMUNE_EQTL_ARCHIVE"]),
        "JCTF_NHA000193":jctf_readme(),
        "EUR_BLUEPRINT_MONOCYTE":query_head(URLS["EUR_BLUEPRINT_MONOCYTE_ALL_ASSOC"]),
        "GTEX_V8_ALLPAIRS":query_head(URLS["GTEX_V8_WHOLE_BLOOD_ALLPAIRS"]),
    }
    report={
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "scope":"PUBLIC_SOURCE_METADATA_AND_HEAD_ONLY",
        "full_eqtl_variant_level_downloads":False,
        "any_coloc_completed":False,
        "all_sources_still_require_allele_effect_se_variant_n_build_QC":True,
        "results":tests,
    }
    (out/"source_access.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    lines=["# IS eQTL actual public access gates", "",
           "A HTTP 200/HEAD is **not** complete cis-eQTL/statistical eligibility.",
           "",
           "| Source | Access status | Nominal all cis effects |",
           "|---|---|---|"]
    for k,v in tests.items():
        lines.append(f"| {k} | {v.get('access',v.get('status'))} | HOLD |")
    lines.extend(["","MAEEA 2026 is ancestry-matched but original attached data not yet shown accessible.",
                  "NHA000051 (n=105) is public Japanese eQTL, but ZIP source availability does NOT establish full nominal tested SNP–gene coverage.",
                  "JCTF eQTL (n=1,019) is filtered by p-value/PIP: do not assume standard coloc valid.",
                  "GTEx all-associations can require requester-pays project; failure of anonymous HEAD is not biological absence.",
                  "Need actual full cis-region rows with matching variant alleles, effect/SE/MAF and relevant ancestry before GO."])
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md,encoding="utf-8")
    print(md)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as f:
            f.write(md)


if __name__=="__main__":
    main()
