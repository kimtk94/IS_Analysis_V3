#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import zipfile
from collections import defaultdict
from pathlib import Path

SUPP_FILENAME = "NIHMS2102371-supplement-Supplementary_Tables.xlsx"
SUPP_URLS = [
    f"https://pmc.ncbi.nlm.nih.gov/articles/instance/12435990/bin/{SUPP_FILENAME}",
    f"https://pmc.ncbi.nlm.nih.gov/articles/PMC12435990/bin/{SUPP_FILENAME}",
    f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12435990/bin/{SUPP_FILENAME}",
]
# Springer/Nature supplementary media URLs are stable once the article's
# internal media-object number is known, but the MOESM index is not exposed in
# all machine-readable views. Try a bounded set and validate workbook contents.
SUPP_URLS += [
    (
        "https://static-content.springer.com/esm/"
        "art%3A10.1038%2Fs41591-025-03872-8/MediaObjects/"
        f"41591_2025_3872_MOESM{i}_ESM.xlsx"
    )
    for i in range(1, 13)
]

# Candidate-level evidence reproduced directly from Table 1 / main text of:
# Hirohama et al., Nature Medicine 2025, doi:10.1038/s41591-025-03872-8.
# These rows allow Stage 3A to remain informative even when publisher
# supplementary-file endpoints reject scripted downloads.
HIROHAMA_MAIN_TEXT = {
    "ACP1": {
        "kidney_pqtl_present": 1,
        "egfr_gwas_snp": "rs79154857",
        "kidney_pqtl_snp": "rs62114548",
        "kidney_pqtl_alt": "G",
        "kidney_pqtl_beta": 0.986,
        "kidney_pqtl_p": 8.5e-79,
        "kidney_egfr_pph4": 0.999,
        "kidney_smr_p": 5.8e-14,
        "kidney_heidi_p": 0.167,
        "main_text_note": "Kidney pQTL identified; ACP1 was noted as not previously seen in kidney eQTL analyses.",
    },
    "GSTA1": {
        "kidney_pqtl_present": 1,
        "egfr_gwas_snp": "rs6423287",
        "kidney_pqtl_snp": "rs9382146",
        "kidney_pqtl_alt": "A",
        "kidney_pqtl_beta": 0.249,
        "kidney_pqtl_p": 3.0e-22,
        "kidney_egfr_pph4": 0.999,
        "kidney_smr_p": 1.5e-11,
        "kidney_heidi_p": 0.067,
        "main_text_note": "Kidney protein prioritized by both colocalization and SMR/HEIDI for eGFR.",
    },
    "INHBC": {
        "kidney_pqtl_present": 1,
        "egfr_gwas_snp": "rs7964492",
        "kidney_pqtl_snp": "rs7971133",
        "kidney_pqtl_alt": "T",
        "kidney_pqtl_beta": -0.625,
        "kidney_pqtl_p": 1.2e-25,
        "kidney_egfr_pph4": 0.875,
        "kidney_smr_p": 3.9e-16,
        "kidney_heidi_p": 0.086,
        "main_text_note": "Kidney pQTL identified; INHBC was noted as not previously seen in kidney eQTL analyses.",
    },
}
EQTL_META_URL = "https://figshare.com/ndownloader/files/33957947"
EQTL_TUBULE_URL = "https://figshare.com/ndownloader/files/38295906"
EQTL_GLOM_URL = "https://figshare.com/ndownloader/files/38295879"
SUSZTAK_AGREEMENT_URL = "https://susztaklab.com/agree.php"


def run(cmd):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run([str(x) for x in cmd], check=True)


def download(url: str, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file() and dest.stat().st_size > 0:
        return
    part = dest.with_suffix(dest.suffix + ".part")
    run(["curl", "-L", "--fail", "--retry", "5", "--retry-delay", "2", "-o", part, url])
    if not part.is_file() or part.stat().st_size == 0:
        raise RuntimeError(f"download failed or empty: {url}")
    part.replace(dest)


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def fnum(x):
    try:
        if x is None or str(x).strip() in {"", "NA", "NaN", "nan"}:
            return None
        return float(x)
    except Exception:
        return None


def classify_stage2(row):
    status = (row.get("comparison_status") or "").strip()
    h4 = fnum(row.get("max_PP.H4"))
    abf_h4 = fnum(row.get("stage2b_ABF_PP.H4"))
    if status == "ABF_and_SuSiE_shared_signal":
        return "core_shared_signal"
    if status == "SuSiE_nonconverged_or_failed":
        return "susie_unresolved"
    if status == "ABF_not_confirmed_by_SuSiE":
        if h4 is not None and h4 >= 0.5:
            return "abf_supported_susie_weaker"
        return "abf_susie_discordant"
    if (abf_h4 is not None and abf_h4 < 0.2) and (h4 is None or h4 < 0.2):
        return "nonshared_or_unresolved"
    return "unresolved"


def write_tsv(path: Path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as fh:
        if not fields:
            return
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def is_valid_xlsx(path: Path):
    if not path.is_file() or path.stat().st_size < 100_000:
        return False
    if not zipfile.is_zipfile(path):
        return False
    try:
        with zipfile.ZipFile(path) as zf:
            names = set(zf.namelist())
            return "[Content_Types].xml" in names and any(
                n.startswith("xl/worksheets/") for n in names
            )
    except zipfile.BadZipFile:
        return False


def workbook_looks_like_supp_tables(path: Path):
    if not is_valid_xlsx(path):
        return False
    try:
        import openpyxl
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        names = list(wb.sheetnames)
        wb.close()
    except Exception:
        return False
    # The published file contains Supplementary Tables 1-30 in separate tabs.
    return len(names) >= 20


def download_supplementary_workbook(dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)

    if workbook_looks_like_supp_tables(dest):
        return dest, "downloaded_or_cached", ""

    dest.unlink(missing_ok=True)
    errors = []
    for url in SUPP_URLS:
        part = dest.with_suffix(dest.suffix + ".part")
        part.unlink(missing_ok=True)
        try:
            run([
                "curl", "-L", "--fail", "--retry", "2", "--retry-delay", "1",
                "-A", "Mozilla/5.0 CKD-Stage3A",
                "-o", part, url
            ])
            if workbook_looks_like_supp_tables(part):
                part.replace(dest)
                return dest, "downloaded", url
            size = part.stat().st_size if part.exists() else 0
            errors.append(f"{url} -> not Supplementary Tables workbook ({size} bytes)")
        except subprocess.CalledProcessError as exc:
            errors.append(f"{url} -> curl exit {exc.returncode}")
        finally:
            part.unlink(missing_ok=True)

    # This is intentionally non-fatal. Publisher/PMC anti-bot behavior can
    # block scripted supplementary downloads while the article itself remains
    # publicly readable. Main-text Table 1 evidence plus independent kidney
    # eQTL datasets are still valid Stage 3A inputs.
    return None, "unavailable_nonfatal", " | ".join(errors)

def sheet_table_number(title: str):
    m = re.search(r"(?:table|supp(?:lementary)?\s*table)?\s*(\d{1,2})",
                  title, flags=re.I)
    return int(m.group(1)) if m else None


def evidence_type_from_table(n):
    if n == 3:
        return "kidney_cis_pqtl_significant"
    if n == 13:
        return "kidney_cis_eqtl_significant_same_study"
    if n in {16, 17, 18}:
        return "kidney_qtl_egfr_colocalization"
    return "supplementary_candidate_hit"


def normalize_cell(v):
    if v is None:
        return ""
    return str(v).strip()


def exact_gene_in_cell(cell: str, genes):
    x = cell.strip().upper()
    if x in genes:
        return x
    toks = [t for t in re.split(r"[\s,;/|]+", x) if t]
    hits = [t for t in toks if t in genes]
    return hits[0] if hits else None


def extract_candidate_supplement_hits(workbook: Path, genes, output_path: Path):
    try:
        import openpyxl
    except ImportError as exc:
        raise RuntimeError("openpyxl is required") from exc

    wb = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    hits = []
    gene_set = set(genes)

    for ws in wb.worksheets:
        rows = ws.iter_rows(values_only=True)
        pre = []
        for _ in range(12):
            try:
                pre.append(next(rows))
            except StopIteration:
                break
        if not pre:
            continue

        header_i = max(
            range(len(pre)),
            key=lambda i: sum(v not in (None, "") for v in pre[i])
        )
        header = [
            normalize_cell(v) or f"col_{j+1}"
            for j, v in enumerate(pre[header_i])
        ]
        table_no = sheet_table_number(ws.title)
        evtype = evidence_type_from_table(table_no)

        all_rows = pre[header_i + 1:]
        all_rows.extend(list(rows))

        for row_index, row in enumerate(all_rows, start=header_i + 2):
            matched = None
            matched_col = None
            for j, v in enumerate(row):
                gene = exact_gene_in_cell(normalize_cell(v), gene_set)
                if gene:
                    matched = gene
                    matched_col = header[j] if j < len(header) else f"col_{j+1}"
                    break
            if not matched:
                continue

            obj = {}
            for j, v in enumerate(row):
                key = header[j] if j < len(header) else f"col_{j+1}"
                if key in obj:
                    key = f"{key}_{j+1}"
                obj[key] = normalize_cell(v)

            hits.append({
                "gene_symbol": matched,
                "sheet": ws.title,
                "table_number": "" if table_no is None else table_no,
                "evidence_type": evtype,
                "matched_column": matched_col,
                "row_index": row_index,
                "row_json": json.dumps(
                    obj, ensure_ascii=False, separators=(",", ":")
                ),
            })

    write_tsv(
        output_path,
        hits,
        [
            "gene_symbol", "sheet", "table_number", "evidence_type",
            "matched_column", "row_index", "row_json"
        ],
    )
    return hits


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(
            path, "rt", encoding="utf-8", errors="replace", newline=""
        )
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def split_line(line: str, delimiter):
    if delimiter == "\t":
        return line.rstrip("\r\n").split("\t")
    return re.split(r"\s+", line.strip())


def filter_eqtl_file(path: Path, dataset: str, sample_n: int,
                     genes, out_path: Path):
    gene_set = set(g.upper() for g in genes)
    matched = []
    counts = defaultdict(int)

    with open_text(path) as fh:
        first = fh.readline()
        if not first:
            raise RuntimeError(f"empty eQTL file: {path}")
        delimiter = "\t" if "\t" in first else None
        header = split_line(first, delimiter)
        lower = [h.lower() for h in header]
        gene_cols = [
            i for i, h in enumerate(lower)
            if "gene" in h or h in {"phenotype_id", "phenotype"}
        ]

        for line_no, line in enumerate(fh, start=2):
            if not line.strip():
                continue
            vals = split_line(line, delimiter)
            if len(vals) < len(header):
                continue

            hit_gene = None
            search_cols = gene_cols if gene_cols else range(len(vals))
            for i in search_cols:
                if i >= len(vals):
                    continue
                gene = exact_gene_in_cell(vals[i], gene_set)
                if gene:
                    hit_gene = gene
                    break

            if not hit_gene:
                continue

            row = {
                header[i]: vals[i] if i < len(vals) else ""
                for i in range(len(header))
            }
            matched.append({
                "dataset": dataset,
                "sample_n": sample_n,
                "gene_symbol": hit_gene,
                "source_line": line_no,
                "row_json": json.dumps(
                    row, ensure_ascii=False, separators=(",", ":")
                ),
            })
            counts[hit_gene] += 1

    write_tsv(
        out_path,
        matched,
        ["dataset", "sample_n", "gene_symbol", "source_line", "row_json"],
    )
    return counts


def aggregate(stage2_rows, supplement_hits, eqtl_counts, output: Path):
    supp = defaultdict(lambda: defaultdict(int))
    for hit in supplement_hits:
        supp[hit["gene_symbol"]][hit["evidence_type"]] += 1

    rows = []
    for row in stage2_rows:
        gene = row["gene_symbol"].upper()
        stage2_class = classify_stage2(row)

        if stage2_class == "core_shared_signal":
            next_action = "kidney_tissue_replication_then_koges"
        elif stage2_class in {
            "abf_susie_discordant", "abf_supported_susie_weaker"
        }:
            next_action = "resolve_signal_architecture_then_kidney_tissue"
        elif stage2_class == "susie_unresolved":
            next_action = "ld_diagnostic_or_conditional_coloc"
        else:
            next_action = "supporting_or_negative_control_evidence"

        mt = HIROHAMA_MAIN_TEXT.get(gene, {})
        rows.append({
            "gene_symbol": gene,
            "stage2_class": stage2_class,
            "stage2_comparison_status": row.get("comparison_status", ""),
            "stage2_abf_h4": row.get("stage2b_ABF_PP.H4", ""),
            "stage2_susie_max_h4": row.get("max_PP.H4", ""),
            "hirohama_maintext_kidney_pqtl": mt.get("kidney_pqtl_present", 0),
            "hirohama_maintext_egfr_gwas_snp": mt.get("egfr_gwas_snp", ""),
            "hirohama_maintext_kidney_pqtl_snp": mt.get("kidney_pqtl_snp", ""),
            "hirohama_maintext_kidney_pqtl_beta": mt.get("kidney_pqtl_beta", ""),
            "hirohama_maintext_kidney_pqtl_p": mt.get("kidney_pqtl_p", ""),
            "hirohama_maintext_egfr_pph4": mt.get("kidney_egfr_pph4", ""),
            "hirohama_maintext_smr_p": mt.get("kidney_smr_p", ""),
            "hirohama_maintext_heidi_p": mt.get("kidney_heidi_p", ""),
            "hirohama_maintext_note": mt.get("main_text_note", ""),
            "kidney_pqtl_significant_supp_hits":
                supp[gene]["kidney_cis_pqtl_significant"],
            "kidney_eqtl_same_study_supp_hits":
                supp[gene]["kidney_cis_eqtl_significant_same_study"],
            "kidney_qtl_egfr_coloc_supp_hits":
                supp[gene]["kidney_qtl_egfr_colocalization"],
            "kidney_eqtl_meta686_hits":
                eqtl_counts["meta686"].get(gene, 0),
            "kidney_eqtl_tubule356_hits":
                eqtl_counts["tubule356"].get(gene, 0),
            "kidney_eqtl_glomerulus303_hits":
                eqtl_counts["glomerulus303"].get(gene, 0),
            "next_action": next_action,
        })

    write_tsv(output, rows)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage2c-default", type=Path, required=True)
    ap.add_argument("--raw-root", type=Path, required=True)
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument("--accept-susztak-terms", action="store_true")
    args = ap.parse_args()

    if not args.accept_susztak_terms:
        raise SystemExit(
            "Susztak Kidney Biobank terms were not accepted. Read "
            + SUSZTAK_AGREEMENT_URL
            + " and rerun with --accept-susztak-terms."
        )

    stage2_rows = list(read_tsv(args.stage2c_default))
    if not stage2_rows:
        raise SystemExit("Stage2C default table is empty")
    genes = [row["gene_symbol"].upper() for row in stage2_rows]

    args.raw_root.mkdir(parents=True, exist_ok=True)
    args.output_root.mkdir(parents=True, exist_ok=True)

    # Remove the obsolete short Europe PMC response from the first
    # implementation if it is still present.
    (args.raw_root / "hirohama2025_europepmc_supplementary.zip").unlink(
        missing_ok=True
    )
    workbook, supplement_status, supplement_detail = (
        download_supplementary_workbook(args.raw_root / SUPP_FILENAME)
    )

    eqtl_sources = {
        "meta686": (
            EQTL_META_URL,
            args.raw_root / "Kidney_eQTL_Meta_S686_Significant.q0.01.txt.gz",
            686,
        ),
        "tubule356": (
            EQTL_TUBULE_URL,
            args.raw_root / "Kidney_eQTL.TubsigeQTLsFormated.txt.gz",
            356,
        ),
        "glomerulus303": (
            EQTL_GLOM_URL,
            args.raw_root / "Kidney_eQTL.GlomsigeQTLsFormated.txt.gz",
            303,
        ),
    }

    for _, (url, path, _) in eqtl_sources.items():
        download(url, path)

    supplement_hits_path = (
        args.output_root / "STAGE3A_HIROHAMA_SUPPLEMENT_HITS.tsv"
    )
    if workbook is not None:
        supplement_hits = extract_candidate_supplement_hits(
            workbook, genes, supplement_hits_path
        )
    else:
        supplement_hits = []
        write_tsv(
            supplement_hits_path,
            [],
            [
                "gene_symbol", "sheet", "table_number", "evidence_type",
                "matched_column", "row_index", "row_json"
            ],
        )
        print(
            "WARNING: Hirohama Supplementary Tables workbook unavailable; "
            "continuing with main-text Table 1 evidence and kidney eQTL datasets.",
            flush=True,
        )

    eqtl_counts = {}
    for name, (_, path, sample_n) in eqtl_sources.items():
        eqtl_counts[name] = filter_eqtl_file(
            path,
            name,
            sample_n,
            genes,
            args.output_root / f"STAGE3A_EQTL_{name.upper()}_HITS.tsv",
        )

    evidence = aggregate(
        stage2_rows,
        supplement_hits,
        eqtl_counts,
        args.output_root / "STAGE3A_KIDNEY_EVIDENCE.tsv",
    )

    provenance = {
        "stage": "CKD Stage 3A kidney-specific evidence",
        "candidate_source": str(args.stage2c_default),
        "candidate_genes": genes,
        "hirohama2025": {
            "pmcid": "PMC12435990",
            "supplementary_filename": SUPP_FILENAME,
            "supplementary_urls": SUPP_URLS,
            "supplement_status": supplement_status,
            "supplement_detail": supplement_detail,
            "workbook": "" if workbook is None else str(workbook),
            "main_text_candidate_evidence": HIROHAMA_MAIN_TEXT,
            "kidney_pqtl_sample_n": 337,
            "same_study_eqtl_sample_n": 315,
        },
        "kidney_eqtl": {
            "meta686": EQTL_META_URL,
            "tubule356": EQTL_TUBULE_URL,
            "glomerulus303": EQTL_GLOM_URL,
        },
        "license_note": {
            "agreement": SUSZTAK_AGREEMENT_URL,
            "raw_data_policy":
                "Raw QTL files stay outside GitHub and are not synced by this runner.",
            "publication_note":
                "Review the current Susztak Kidney Biobank user agreement before publication or redistribution.",
        },
        "interpretation_note":
            "eQTL source files contain significant SNP-gene pairs only; zero hits means no significant hit in that published table, not proof of absent kidney expression or absent cis regulation.",
    }

    (args.output_root / "STAGE3A_PROVENANCE.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )

    print(json.dumps({"genes": genes, "evidence": evidence}, indent=2))
    print(f"CKD_STAGE3A_PASS output={args.output_root}")


if __name__ == "__main__":
    main()
