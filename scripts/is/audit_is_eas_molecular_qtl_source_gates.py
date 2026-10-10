#!/usr/bin/env python3
"""Scientific fail-closed ledger for East-Asian G0022 molecular QTL sources.

Public source descriptions establish release selection policy, but not actual
SNP/site availability or matched GWAS–QTL signal-level colocalization.
--probe-remote fetches only fixed-size HTTP byte ranges, NEVER archive bodies.
"""
import argparse
import json
import struct
import urllib.request
from datetime import date
from pathlib import Path

DDBJ_ARCHIVE = ("https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/"
                "E-GEAD-000/E-GEAD-420/E-GEAD-420.processed.zip")
DDBJ_OFFICIAL_BYTES = 43507433191

SOURCES = [
    dict(dataset="NHA000172", name="JCTF_v2_Japanese_eQTL_sQTL",
         ancestry="Japanese_COVID19", assay="cis_eQTL_sQTL", n=465,
         public_release=True, release_selection="cis_p_less_than_0.05_only",
         includes_non_significant_tests=False,
         publication_ready=False, access="unrestricted_summary_749MB",
         ref="https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000172",
         reason="NBDC states variant–gene cis results were selected for p<0.05. A filtered source is not the entire tested cis universe."),
    dict(dataset="NHA000193", name="JCTF_v3_Japanese_eQTL_pQTL",
         ancestry="Japanese_COVID19", assay="cis_eQTL_pQTL", n=1405,
         public_release=True, release_selection="p_less_than_0.05_or_PIP_greater_than_0.001",
         includes_non_significant_tests=False,
         publication_ready=False, access="unrestricted_summary_1.7GB",
         ref="https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA000193",
         reason="Rows are P/PIP thresholded; rs671 association evidence does not supply full locus coloc inputs."),
    dict(dataset="E-GEAD-420", name="ImmuNexUT_Japanese_nominal_eQTL",
         ancestry="Japanese_immune_mixed_healthy_immune_disease", assay="cis_eQTL_nominal", n=416,
         public_release=True, release_selection="nominal_all_pairs_including_non_significant",
         includes_non_significant_tests=True,
         publication_ready=False, access="unrestricted_43.5GB_nested_zip_of_tar_gz",
         ref="https://www.immunexut.org/faqs",
         reason="Official source advertises full nominal eQTL, but genotype-aligned rs671 coverage, complete gene/cell-type rows, source LD and effect allele are not yet checked. ZIP holds one nested tar.gz, not individually seekable gene files."),
    dict(dataset="ZENODO_21296030", name="MAEEA_Asian_whole_blood_cis_eQTL",
         ancestry="East_Asian_Chinese_Japanese", assay="cis_eQTL", n=2024,
         public_release=False, release_selection="all_SNP_gene_pairs_described_access_restricted",
         includes_non_significant_tests=None,
         publication_ready=False, access="restricted_in_prior_server_metadata_live_unverified",
         ref="https://zenodo.org/records/21296030",
         reason="Complete summary is described, but no authorized downloadable file verified. Recheck access on future dates."),
]

def gate(source):
    if source["includes_non_significant_tests"] is False:
        return "BLOCKED_THRESHOLD_SELECTED"
    if source["public_release"] is False:
        return "BLOCKED_ACCESS"
    if source["includes_non_significant_tests"] is True:
        return "SOURCE_PROMISING_GENOTYPE_COVERAGE_UNKNOWN"
    return "SOURCE_COMPLETENESS_UNKNOWN"


def assess(sources):
    results = []
    for source in sources:
        item = dict(source)
        item["gate"] = gate(source)
        item["rs671_and_four_cs_source_coverage_verified"] = False
        item["complete_g0022_gene_by_cell_type_statistics_verified"] = False
        item["ancestry_and_effect_alleles_harmonized"] = False
        item["qtl_independent_signal_and_signed_ld_verified"] = False
        item["valid_gwas_qtl_colocalizations"] = 0
        item["can_publish_causal_mediator"] = False
        if item["publication_ready"] or item["can_publish_causal_mediator"]:
            raise ValueError("Unverified QTL evidence cannot be publication promoted.")
        results.append(item)
    return results


def archive_index_probe(url=DDBJ_ARCHIVE, max_bytes=262144, timeout=12):
    """HEAD plus 2 bounded Range GETs; never download entire 43.5GB ZIP."""
    def request(method, headers=None):
        req = urllib.request.Request(url, method=method, headers={
            "User-Agent": "G0022-scientific-archive-index-audit/1.0",
            **(headers or {})})
        return urllib.request.urlopen(req, timeout=timeout)
    with request("HEAD") as r:
        size = int(r.headers["Content-Length"])
        accept_ranges = r.headers.get("Accept-Ranges", "")
    if size < 1000 or size > 100_000_000_000:
        raise ValueError("Unexpected archive size; refusing request")
    if "bytes" not in accept_ranges.lower():
        raise ValueError("Server does not advertise supported byte ranges")
    def bounded(start, stop):
        expected = stop - start + 1
        if expected > max_bytes or expected <= 0:
            raise ValueError("Range request exceeds maximum allowed bytes")
        with request("GET", {"Range": f"bytes={start}-{stop}"}) as r:
            if r.status != 206 or not r.headers.get("Content-Range", "").startswith(f"bytes {start}-{stop}/"):
                raise ValueError("Server did not honor Range request (refusing full download)")
            buf = r.read(expected + 1)
        if len(buf) != expected:
            raise ValueError("Ranged response length mismatch")
        return buf
    head = bounded(0, 4095)
    tail_start = size - min(size, max_bytes)
    tail = bounded(tail_start, size - 1)
    if not head.startswith(b"PK\x03\x04"):
        raise ValueError("ZIP local file header missing")
    local_method = struct.unpack_from("<H", head, 8)[0]
    fn_len = struct.unpack_from("<H", head, 26)[0]
    local_name = head[30:30 + fn_len].decode("utf8")
    central_idx = tail.rfind(b"PK\x01\x02")
    end_idx = tail.rfind(b"PK\x05\x06")
    if central_idx < 0 or end_idx < 0:
        raise ValueError("ZIP central directory not located in bounded tail")
    name_len, extra_len, comment_len = struct.unpack_from("<HHH", tail, central_idx + 28)
    central_name = tail[central_idx + 46:central_idx + 46 + name_len].decode("utf8")
    if local_name != central_name:
        raise ValueError("Header and central ZIP member names mismatch")
    return {
        "source": url, "http_content_length": size,
        "accept_ranges": accept_ranges,
        "first_local_member": local_name,
        "last_central_member": central_name,
        "local_zip_compression_method": local_method,
        "zip64_end_record_detected": b"PK\x06\x06" in tail,
        "only_first_and_last_members_verified": True,
        "targeted_gene_or_variant_seekable": False if local_name.endswith(".tar.gz") else None,
        "bytes_fetched_body_upper_bound": len(head) + len(tail),
        "disposition": "NESTED_ARCHIVE_NO_RANDOM_ACCESS_TO_RS671_BY_ZIP_INDEX" if local_name.endswith(".tar.gz") else "ARCHIVE_LAYOUT_UNDETERMINED",
        "warning": "Never infer that the actual source contains rs671 or every tested variant until local archive evidence is checked.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--probe-remote", action="store_true",
                        help="Read only 4KB first, <=256KB last of DDBJ ZIP; not full download.")
    args = parser.parse_args()
    if not args.outdir.exists():
        args.outdir.mkdir(parents=True)
    output = {"audit": "IS_G0022_EAS_MOLECULAR_QTL_ACCESS_AND_SOURCE_COMPLETENESS_V1",
              "date": str(date.today()), "candidate_locus": "G0022 chr12 ALDH2 rs671",
              "source_count": len(SOURCES), "sources": assess(SOURCES),
              "valid_new_colocs": 0, "causal_target_confirmed": False}
    if args.probe_remote:
        try:
            output["imm unexut_archive_range_probe".replace(" ", "")] = archive_index_probe()
        except Exception as e:
            output["immunexut_archive_range_probe_error"] = f"{type(e).__name__}: {e}"
    path = args.outdir / "IS_G0022_EAS_QTL_SOURCE_GATES.json"
    path.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
