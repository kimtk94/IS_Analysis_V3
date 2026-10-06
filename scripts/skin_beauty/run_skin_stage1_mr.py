#!/usr/bin/env python3
"""
SKIN BEAUTY MR — Stage 1
UKB-PPP ST16 cis-independent pQTL -> Korean facial wrinkle GWAS

Primary outcomes:
  ebi-a-GCST90094903  Facial wrinkles (under eye)
  ebi-a-GCST90094904  Facial wrinkles (crow's feet)

Design:
- Uses the existing MR-ready ST16 PRIMARY_NON_MHC canonical exposure file when available.
- Primary instrument QC follows the validated ST16 registry: cis, conditional P < 5e-8,
  F > 10, non-MHC.
- Queries exact rsIDs or GRCh37 chr:position keys from OpenGWAS /associations with proxies=0.
- Uses OpenGWAS-safe batches with N(id)*N(variant)<=64, retry/backoff and persistent cache.
- Cross-ancestry harmonisation is deliberately conservative:
  palindromic A/T and C/G SNPs are excluded from primary MR.
- Single-SNP proteins: Wald ratio.
- >=2 SNPs: IVW multiplicative random-effects primary estimate.
- >=3 SNPs: MR-Egger slope/intercept diagnostic.
- Also reports weighted-median point estimate.
- Writes FDR-adjusted protein-level results.

No external Python packages are required.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import random
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any, Optional

API_BASE = "https://api.opengwas.io/api"

OUTCOMES = {
    "ebi-a-GCST90094903": "facial_wrinkles_under_eye",
    "ebi-a-GCST90094904": "facial_wrinkles_crows_feet",
}

ALIASES = {
    "protein": ["protein_id", "aptamer_id", "aptamer", "protein", "target_id", "assay_id"],
    "gene": ["gene_symbol", "gene", "symbol", "target"],
    "uniprot": ["uniprot", "uniprot_id", "uniprotid"],
    "variant": ["rsid", "snp", "variant", "variant_id", "rs_id"],
    "snp_source": ["snp"],
    "variant_id_source": ["variant_id_hg19", "variant_id"],
    "chr": ["chr_hg19", "chrom_hg19", "chromosome", "chr", "chrom"],
    "pos": ["pos_hg19", "position", "pos", "bp"],
    "ea": ["effect_allele", "effect_allele_exposure", "alt_allele", "ea", "a1", "tested_allele"],
    "oa": ["other_allele", "other_allele_exposure", "ref_allele", "nea", "oa", "a2", "non_effect_allele"],
    "eaf": ["eaf", "eaf_exposure", "alt_freq", "effect_allele_frequency", "a1freq", "af"],
    "beta": ["beta_exposure", "beta_cond", "beta", "effect", "b"],
    "se": ["se_exposure", "se_cond", "se", "stderr", "standard_error"],
    "p": ["p_exposure", "p_value", "pval", "p", "pvalue"],
    "mlogp": ["minus_log10p_exposure", "minus_log10p_cond", "minus_log10p", "mlogp", "neglog10p"],
    "f": ["f_stat", "fstat", "f_statistic", "f"],
    "f_gt_10": ["f_gt_10", "f_gt10"],
    "sig": ["gwas_significant", "genomewide_significant", "significant", "p_lt_1_7e11", "p_lt_17e11"],
    "mhc": ["in_extended_mhc_hg19", "mhc", "is_mhc", "in_mhc"],
    "cis_trans": ["cis_trans"],
    "primary_non_mhc": ["primary_non_mhc"],
}

TRUE_SET = {"1", "true", "t", "yes", "y", "pass"}
FALSE_SET = {"0", "false", "f", "no", "n", "fail", ""}
VALID_BASES = {"A", "C", "G", "T"}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def norm(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(x).strip().lower()).strip("_")


def safe_float(x: Any) -> Optional[float]:
    if x is None:
        return None
    try:
        v = float(str(x).strip())
        if not math.isfinite(v):
            return None
        return v
    except Exception:
        return None


def as_bool(x: Any) -> Optional[bool]:
    if x is None:
        return None
    s = str(x).strip().lower()
    if s in TRUE_SET:
        return True
    if s in FALSE_SET:
        return False
    return None


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="replace", newline="")
    return path.open(mode, encoding="utf-8", errors="replace", newline="")


def detect_delimiter(first: str) -> str:
    return "\t" if first.count("\t") >= first.count(",") else ","


def detect_columns(fieldnames: list[str]) -> dict[str, Optional[str]]:
    by_norm = {norm(c): c for c in fieldnames}
    out = {}
    for role, candidates in ALIASES.items():
        out[role] = next((by_norm[a] for a in candidates if a in by_norm), None)
    return out


def normalize_chromosome(x: str) -> str:
    s = re.sub(r"^chr", "", str(x).strip(), flags=re.I).upper()
    if s == "M":
        s = "MT"
    if re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT)", s):
        return s
    return ""


def parse_chrpos_from_source(*values: str) -> tuple[Optional[str], Optional[int], str]:
    """Recover GRCh37 chr/pos from canonical source identifiers.

    Accepted examples include `X:15692297_CGT_C` and
    `X:15692297:CGT:C:imp:v1`. Only chromosome and integer position are
    recovered here; alleles continue to come from the canonical exposure
    columns and are validated during harmonisation.
    """
    pat = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT|M)[:_](\d+)(?=[:_]|$)", re.I)
    for label, value in values:
        text = str(value or "").strip()
        if not text:
            continue
        m = pat.match(text)
        if not m:
            continue
        chrom = normalize_chromosome(m.group(1))
        try:
            pos = int(m.group(2))
        except Exception:
            continue
        if chrom and pos > 0:
            return chrom, pos, label
    return None, None, ""


def choose_exposure(default_path: str, root: Path) -> Path:
    p = Path(default_path)
    if p.exists():
        return p
    fallbacks = [
        Path("/srv/is-analysis/results/metabolic_resilience/stage1_pqtl/instruments/UKBPPP_CIS_INSTRUMENTS_PRIMARY_NON_MHC.tsv.gz"),
        Path("/srv/is-analysis/data/ckd/analysis_ready/UKBPPP_ST16_cis_independent.tsv.gz"),
        root / "data/skin_beauty/stage0_gwas/UKBPPP_CIS_INSTRUMENTS.tsv.gz",
    ]
    for x in fallbacks:
        if x.exists():
            return x
    raise FileNotFoundError(
        "ST16 exposure not found. Expected: "
        "/srv/is-analysis/data/ckd/analysis_ready/UKBPPP_ST16_cis_independent.tsv.gz"
    )


def read_exposure(path: Path, p_threshold: float) -> tuple[list[dict], dict]:
    with open_text(path) as fh:
        first = fh.readline()
        if not first:
            raise RuntimeError(f"Empty exposure file: {path}")
        delim = detect_delimiter(first)
        fieldnames = first.rstrip("\r\n").split(delim)
        cmap = detect_columns(fieldnames)

        required = ["variant", "ea", "oa", "beta", "se"]
        missing = [x for x in required if cmap.get(x) is None]
        if missing:
            raise RuntimeError(
                f"Exposure is missing required canonical fields {missing}. "
                f"Header={fieldnames}"
            )

        reader = csv.DictReader(fh, fieldnames=fieldnames, delimiter=delim)
        rows = []
        no_query_key_rows = []
        counters = defaultdict(int)

        for r in reader:
            counters["raw_rows"] += 1

            # Raw ST16 is a published conditional cis-pQTL table.  When the
            # source exposes cis_trans explicitly, retain cis signals only.
            if cmap.get("cis_trans"):
                cis_trans = str(r.get(cmap["cis_trans"], "")).strip().lower()
                if cis_trans and cis_trans != "cis":
                    counters["drop_non_cis"] += 1
                    continue

            # The preferred canonical file is already PRIMARY_NON_MHC.
            # If this flag is present, enforce it rather than re-deriving it.
            if cmap.get("primary_non_mhc"):
                pn = as_bool(r.get(cmap["primary_non_mhc"]))
                if pn is False:
                    counters["drop_not_primary_non_mhc"] += 1
                    continue

            raw_variant = str(r.get(cmap["variant"], "")).strip()
            snp_source = str(r.get(cmap["snp_source"], "")).strip() if cmap.get("snp_source") else ""
            variant_id_source = str(r.get(cmap["variant_id_source"], "")).strip() if cmap.get("variant_id_source") else ""
            valid_rsid = bool(re.fullmatch(r"rs\d+", raw_variant, flags=re.I))
            rsid = raw_variant.lower() if valid_rsid else ""

            chr_v = normalize_chromosome(r.get(cmap["chr"], "") if cmap.get("chr") else "")
            pos_v = safe_float(r.get(cmap["pos"])) if cmap.get("pos") else None
            pos_i = int(pos_v) if pos_v is not None and pos_v > 0 else None
            query_key_recovery_source = "canonical_chr_pos" if chr_v and pos_i is not None else ""

            # Some ST16 X-chromosome indels have a blank chr_hg19 field while
            # `snp`/`variant_id_hg19` still carry a complete GRCh37 X:position
            # identifier. Recover only the coordinate here; canonical alleles
            # remain authoritative for harmonisation.
            if not chr_v or pos_i is None:
                rec_chr, rec_pos, rec_source = parse_chrpos_from_source(
                    ("snp_source", snp_source),
                    ("variant_id_hg19_source", variant_id_source),
                    ("rsid_source", raw_variant),
                )
                if not chr_v and rec_chr:
                    chr_v = rec_chr
                    counters["recovered_missing_chr_from_source"] += 1
                if pos_i is None and rec_pos is not None:
                    pos_i = rec_pos
                    counters["recovered_missing_pos_from_source"] += 1
                if chr_v and pos_i is not None and rec_source:
                    query_key_recovery_source = rec_source
                    counters[f"query_key_recovered_from_{rec_source}"] += 1

            if valid_rsid:
                query_variant = rsid
                query_key_type = "rsid"
                query_key_recovery_source = "rsid"
            elif chr_v and pos_i is not None:
                query_variant = f"{chr_v}:{pos_i}"
                query_key_type = "chrpos"
                counters["query_key_chrpos_fallback"] += 1
            else:
                counters["drop_no_query_key"] += 1
                no_query_key_rows.append({
                    "row_number_data": counters["raw_rows"],
                    "instrument_id": str(r.get("instrument_id", "")).strip(),
                    "protein_id": str(r.get(cmap["protein"], "")).strip() if cmap.get("protein") else "",
                    "gene_symbol": str(r.get(cmap["gene"], "")).strip() if cmap.get("gene") else "",
                    "rsid_source": raw_variant,
                    "snp_source": snp_source,
                    "variant_id_hg19_source": variant_id_source,
                    "chr_hg19_source": chr_v,
                    "pos_hg19_source": str(r.get(cmap["pos"], "")).strip() if cmap.get("pos") else "",
                    "effect_allele": str(r.get(cmap["ea"], "")).strip() if cmap.get("ea") else "",
                    "other_allele": str(r.get(cmap["oa"], "")).strip() if cmap.get("oa") else "",
                    "eaf": str(r.get(cmap["eaf"], "")).strip() if cmap.get("eaf") else "",
                    "beta_exposure": str(r.get(cmap["beta"], "")).strip() if cmap.get("beta") else "",
                    "se_exposure": str(r.get(cmap["se"], "")).strip() if cmap.get("se") else "",
                    "p_exposure": str(r.get(cmap["p"], "")).strip() if cmap.get("p") else "",
                    "minus_log10p_exposure": str(r.get(cmap["mlogp"], "")).strip() if cmap.get("mlogp") else "",
                })
                continue

            bx = safe_float(r.get(cmap["beta"]))
            sex = safe_float(r.get(cmap["se"]))
            if bx is None or sex is None or sex <= 0 or bx == 0:
                counters["drop_bad_beta_se"] += 1
                continue

            fstat = safe_float(r.get(cmap["f"])) if cmap.get("f") else None
            if fstat is None:
                fstat = (bx / sex) ** 2
            if fstat <= 10:
                counters["drop_f_le_10"] += 1
                continue

            if cmap.get("f_gt_10"):
                fflag = as_bool(r.get(cmap["f_gt_10"]))
                if fflag is False:
                    counters["drop_f_flag"] += 1
                    continue

            significant = None
            if cmap.get("sig"):
                significant = as_bool(r.get(cmap["sig"]))
            elif cmap.get("p"):
                pv = safe_float(r.get(cmap["p"]))
                if pv is not None:
                    significant = pv < p_threshold
            elif cmap.get("mlogp"):
                mp = safe_float(r.get(cmap["mlogp"]))
                if mp is not None:
                    significant = mp > -math.log10(p_threshold)

            if significant is False:
                counters["drop_not_significant"] += 1
                continue

            in_mhc = None
            if cmap.get("mhc"):
                in_mhc = as_bool(r.get(cmap["mhc"]))

            if in_mhc is None and chr_v == "6" and pos_i is not None:
                in_mhc = 25_000_000 <= pos_i <= 34_000_000

            if in_mhc is True:
                counters["drop_mhc"] += 1
                continue

            ea = str(r.get(cmap["ea"], "")).strip().upper()
            oa = str(r.get(cmap["oa"], "")).strip().upper()
            if not ea or not oa or ea == oa:
                counters["drop_bad_alleles"] += 1
                continue

            protein = str(r.get(cmap["protein"], "")).strip() if cmap.get("protein") else ""
            gene = str(r.get(cmap["gene"], "")).strip() if cmap.get("gene") else ""
            uniprot = str(r.get(cmap["uniprot"], "")).strip() if cmap.get("uniprot") else ""

            if not protein:
                protein = gene if gene else uniprot
            if not protein:
                counters["drop_missing_protein"] += 1
                continue

            eaf = safe_float(r.get(cmap["eaf"])) if cmap.get("eaf") else None

            rows.append({
                "protein_id": protein,
                "gene_symbol": gene,
                "uniprot": uniprot,
                "rsid": rsid,
                "source_variant": raw_variant,
                "snp_source": snp_source,
                "variant_id_hg19_source": variant_id_source,
                "query_variant": query_variant,
                "query_key_type": query_key_type,
                "query_key_recovery_source": query_key_recovery_source,
                "chr_hg19": chr_v,
                "pos_hg19": pos_i if pos_i is not None else "",
                "effect_allele_exposure": ea,
                "other_allele_exposure": oa,
                "eaf_exposure": eaf,
                "beta_exposure": bx,
                "se_exposure": sex,
                "f_stat": fstat,
            })
            counters["eligible_rows_pre_dedup"] += 1

    # Preserve the validated canonical instrument rows. Do not silently drop
    # instruments just because they lack an rsID; chr:pos is a valid OpenGWAS key.
    pair_counts = defaultdict(int)
    for r in rows:
        pair_counts[(r["protein_id"], r["query_variant"])] += 1
    counters["duplicate_protein_query_pairs"] = sum(v - 1 for v in pair_counts.values() if v > 1)
    counters["eligible_rows"] = len(rows)
    counters["eligible_proteins"] = len({r["protein_id"] for r in rows})
    counters["unique_query_variants"] = len({r["query_variant"] for r in rows})
    counters["rsid_query_rows"] = sum(r["query_key_type"] == "rsid" for r in rows)
    counters["chrpos_query_rows"] = sum(r["query_key_type"] == "chrpos" for r in rows)
    counters["unique_rsid_queries"] = len({r["query_variant"] for r in rows if r["query_key_type"] == "rsid"})
    counters["unique_chrpos_queries"] = len({r["query_variant"] for r in rows if r["query_key_type"] == "chrpos"})

    qc = {
        "exposure_path": str(path),
        "resolved_path": str(path.resolve()),
        "columns": fieldnames,
        "canonical_map": cmap,
        "p_threshold": p_threshold,
        **dict(counters),
        "no_query_key_rows": no_query_key_rows,
    }
    return rows, qc


def request_json(path: str, params: list[tuple[str, Any]], jwt: str,
                 retry_max: int = 5, timeout: int = 180) -> tuple[Any, dict]:
    url = API_BASE + path + "?" + urllib.parse.urlencode(params, doseq=True)
    last_error = None

    for attempt in range(retry_max):
        req = urllib.request.Request(
            url,
            data=b"",
            method="POST",
            headers={
                "Authorization": f"Bearer {jwt}",
                "Accept": "application/json",
                "User-Agent": "skin-beauty-mr/1.0",
            },
        )
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                raw = resp.read().decode("utf-8", errors="replace")
                headers = {k.lower(): v for k, v in resp.headers.items()}
                return json.loads(raw), headers
        except urllib.error.HTTPError as e:
            body = e.read().decode("utf-8", errors="replace")
            last_error = RuntimeError(f"HTTP {e.code} {path}: {body[:500]}")
            if e.code == 429:
                wait = e.headers.get("Retry-After")
                try:
                    wait_s = min(max(float(wait), 1.0), 600.0)
                except Exception:
                    wait_s = 60.0
                time.sleep(wait_s)
                continue
            if e.code in {500, 502, 503, 504}:
                time.sleep(min(5 * (2 ** attempt), 90))
                continue
            raise last_error
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            last_error = e
            time.sleep(min(5 * (2 ** attempt), 90))

    raise RuntimeError(f"OpenGWAS request failed after {retry_max} attempts: {last_error}")


def flatten_rows(obj: Any) -> list[dict]:
    if isinstance(obj, list):
        return [x for x in obj if isinstance(x, dict)]
    if isinstance(obj, dict):
        for k in ("data", "results", "result"):
            v = obj.get(k)
            if isinstance(v, list):
                return [x for x in v if isinstance(x, dict)]
        # Some APIs may return a dict keyed by integer-like record ids.
        vals = [v for v in obj.values() if isinstance(v, dict)]
        if vals:
            return vals
    return []


def cache_key(variants: list[str], outcomes: list[str]) -> str:
    s = "|".join(outcomes) + "::" + "|".join(variants)
    return hashlib.sha256(s.encode()).hexdigest()[:20]


def query_batches(unique_variants: list[str], out_dir: Path, jwt: str,
                  batch_size: int, retry_max: int, pause: float) -> tuple[list[dict], dict]:
    cache_dir = out_dir / "query_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    errors = []
    allowance = {}
    outcomes = list(OUTCOMES.keys())

    # OpenGWAS /associations enforces:
    #   N(id) * N(variant) <= 64
    # Therefore 2 outcomes imply <=32 variants/request.
    api_assoc_product_limit = 64
    requested_batch_size = max(1, int(batch_size))
    max_variants_per_request = max(
        1, api_assoc_product_limit // max(1, len(outcomes))
    )
    effective_batch_size = min(
        requested_batch_size, max_variants_per_request
    )

    if effective_batch_size != requested_batch_size:
        print(
            "OpenGWAS batch-size cap: "
            f"requested={requested_batch_size}, "
            f"n_outcomes={len(outcomes)}, "
            f"effective={effective_batch_size} "
            f"(N_id*N_variant <= {api_assoc_product_limit})"
        )

    n_batches = math.ceil(
        len(unique_variants) / effective_batch_size
    )

    for bi in range(n_batches):
        batch = unique_variants[
            bi * effective_batch_size:
            (bi + 1) * effective_batch_size
        ]
        key = cache_key(batch, outcomes)
        cache_file = cache_dir / f"assoc_{bi:04d}_{key}.json"

        if cache_file.exists():
            try:
                obj = json.loads(cache_file.read_text(encoding="utf-8"))
                rows = flatten_rows(obj)
                all_rows.extend(rows)
                print(f"[{bi+1}/{n_batches}] cache variants={len(batch)} rows={len(rows)}")
                continue
            except Exception:
                pass

        params = [("variant", x) for x in batch]
        params += [("id", x) for x in outcomes]
        params += [("proxies", 0)]

        try:
            obj, headers = request_json("/associations", params, jwt, retry_max=retry_max)
            cache_file.write_text(json.dumps(obj, ensure_ascii=False), encoding="utf-8")
            rows = flatten_rows(obj)
            all_rows.extend(rows)
            allowance = {
                k: v for k, v in headers.items()
                if k.startswith("x-allowance") or k == "retry-after"
            }
            print(
                f"[{bi+1}/{n_batches}] api variants={len(batch)} rows={len(rows)} "
                f"allowance_remaining={allowance.get('x-allowance-remaining','NA')}"
            )
        except Exception as e:
            errors.append({
                "batch": bi,
                "n_variants": len(batch),
                "first_variant": batch[0] if batch else "",
                "last_variant": batch[-1] if batch else "",
                "error": str(e),
            })
            print(f"[{bi+1}/{n_batches}] ERROR {e}", file=sys.stderr)

        if pause > 0:
            time.sleep(pause)

    qc = {
        "requested_unique_variants": len(unique_variants),
        "n_outcomes": len(outcomes),
        "api_assoc_product_limit": api_assoc_product_limit,
        "requested_batch_size": requested_batch_size,
        "effective_batch_size": effective_batch_size,
        "batch_size": effective_batch_size,
        "n_batches": n_batches,
        "raw_association_rows": len(all_rows),
        "n_failed_batches": len(errors),
        "errors": errors,
        "last_allowance_headers": allowance,
    }
    return all_rows, qc


def getv(d: dict, *names: str):
    nmap = {norm(k): k for k in d.keys()}
    for n in names:
        if norm(n) in nmap:
            return d.get(nmap[norm(n)])
    return None


def canonical_outcome_row(r: dict) -> Optional[dict]:
    oid = str(getv(r, "id", "outcome", "gwas_id") or "").strip()
    rsid_raw = str(getv(r, "rsid", "snp", "variant") or "").strip().lower()
    rsid = rsid_raw if re.fullmatch(r"rs\d+", rsid_raw, flags=re.I) else ""
    if oid not in OUTCOMES:
        return None
    beta = safe_float(getv(r, "beta", "b"))
    se = safe_float(getv(r, "se", "standard_error", "stderr"))
    if beta is None or se is None or se <= 0:
        return None
    chr_raw = str(getv(r, "chr", "chromosome") or "").strip()
    chr_raw = re.sub(r"^chr", "", chr_raw, flags=re.I)
    pos_raw = safe_float(getv(r, "position", "pos"))
    pos_int = int(pos_raw) if pos_raw is not None and pos_raw > 0 else None
    chrpos = f"{chr_raw}:{pos_int}" if chr_raw and pos_int is not None else ""
    if not rsid and not chrpos:
        return None
    return {
        "outcome_id": oid,
        "analysis_label": OUTCOMES[oid],
        "trait": str(getv(r, "trait") or ""),
        "rsid": rsid,
        "chr_outcome": chr_raw,
        "pos_outcome": pos_int if pos_int is not None else "",
        "chrpos_outcome": chrpos,
        "effect_allele_outcome": str(getv(r, "ea", "effect_allele") or "").upper(),
        "other_allele_outcome": str(getv(r, "nea", "other_allele", "oa") or "").upper(),
        "eaf_outcome": safe_float(getv(r, "eaf", "effect_allele_frequency")),
        "beta_outcome_raw": beta,
        "se_outcome": se,
        "p_outcome": safe_float(getv(r, "p", "pval", "p_value")),
        "n_outcome": safe_float(getv(r, "n", "sample_size")),
    }


def is_palindromic(a: str, b: str) -> bool:
    return len(a) == 1 and len(b) == 1 and {a, b} in ({"A", "T"}, {"C", "G"})


def complement(a: str) -> Optional[str]:
    if len(a) != 1 or a not in COMPLEMENT:
        return None
    return COMPLEMENT[a]


def harmonize_one(ex: dict, oy: dict) -> tuple[bool, str, Optional[float], Optional[float]]:
    ea = ex["effect_allele_exposure"].upper()
    oa = ex["other_allele_exposure"].upper()
    ye = oy["effect_allele_outcome"].upper()
    yo = oy["other_allele_outcome"].upper()

    if not ye or not yo:
        return False, "outcome_alleles_missing", None, None

    # Cross-ancestry primary analysis: drop all palindromic SNPs rather than
    # inferring orientation from ancestry-dependent allele frequencies.
    if is_palindromic(ea, oa):
        return False, "palindromic_dropped_cross_ancestry", None, None

    by = oy["beta_outcome_raw"]
    eaf_y = oy["eaf_outcome"]

    if ea == ye and oa == yo:
        return True, "direct", by, eaf_y
    if ea == yo and oa == ye:
        return True, "swap", -by, (1 - eaf_y if eaf_y is not None else None)

    cea = complement(ea)
    coa = complement(oa)
    if cea is not None and coa is not None:
        if cea == ye and coa == yo:
            return True, "strand_direct", by, eaf_y
        if cea == yo and coa == ye:
            return True, "strand_swap", -by, (1 - eaf_y if eaf_y is not None else None)

    return False, "allele_mismatch", None, None


def p_two_sided_z(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def gammaincc(a: float, x: float) -> float:
    """Regularized upper incomplete gamma Q(a,x), standard-library implementation."""
    if a <= 0 or x < 0:
        return float("nan")
    if x == 0:
        return 1.0
    eps = 3e-14
    itmax = 1000
    fpmin = 1e-300

    if x < a + 1.0:
        # Series for P(a,x), return 1-P.
        ap = a
        summ = 1.0 / a
        delta = summ
        for _ in range(itmax):
            ap += 1.0
            delta *= x / ap
            summ += delta
            if abs(delta) < abs(summ) * eps:
                break
        log_pref = -x + a * math.log(x) - math.lgamma(a)
        p = summ * math.exp(log_pref)
        return min(max(1.0 - p, 0.0), 1.0)

    # Continued fraction for Q(a,x).
    b = x + 1.0 - a
    c = 1.0 / fpmin
    d = 1.0 / max(b, fpmin)
    h = d
    for i in range(1, itmax + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < fpmin:
            d = fpmin
        c = b + an / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    log_pref = -x + a * math.log(x) - math.lgamma(a)
    q = math.exp(log_pref) * h
    return min(max(q, 0.0), 1.0)


def chi2_sf(x: float, df: int) -> Optional[float]:
    if df <= 0 or x < 0:
        return None
    return gammaincc(df / 2.0, x / 2.0)


def weighted_median(values: list[float], weights: list[float]) -> Optional[float]:
    pairs = sorted((v, w) for v, w in zip(values, weights) if w > 0 and math.isfinite(v) and math.isfinite(w))
    if not pairs:
        return None
    total = sum(w for _, w in pairs)
    acc = 0.0
    for v, w in pairs:
        acc += w
        if acc >= total / 2.0:
            return v
    return pairs[-1][0]


def wls_egger(rows: list[dict]) -> dict:
    k = len(rows)
    if k < 3:
        return {
            "egger_beta": None, "egger_se": None, "egger_p": None,
            "egger_intercept": None, "egger_intercept_se": None, "egger_intercept_p": None,
        }

    sw = sx = sy = sxx = sxy = 0.0
    for r in rows:
        w = 1.0 / (r["se_outcome"] ** 2)
        x = r["beta_exposure"]
        y = r["beta_outcome_aligned"]
        sw += w
        sx += w * x
        sy += w * y
        sxx += w * x * x
        sxy += w * x * y

    det = sw * sxx - sx * sx
    if det <= 0:
        return {
            "egger_beta": None, "egger_se": None, "egger_p": None,
            "egger_intercept": None, "egger_intercept_se": None, "egger_intercept_p": None,
        }

    intercept = (sxx * sy - sx * sxy) / det
    slope = (sw * sxy - sx * sy) / det
    q = 0.0
    for r in rows:
        w = 1.0 / (r["se_outcome"] ** 2)
        resid = r["beta_outcome_aligned"] - intercept - slope * r["beta_exposure"]
        q += w * resid * resid

    phi = max(1.0, q / max(k - 2, 1))
    var_intercept = phi * sxx / det
    var_slope = phi * sw / det
    se_i = math.sqrt(max(var_intercept, 0))
    se_s = math.sqrt(max(var_slope, 0))

    return {
        "egger_beta": slope,
        "egger_se": se_s,
        "egger_p": p_two_sided_z(slope / se_s) if se_s > 0 else None,
        "egger_intercept": intercept,
        "egger_intercept_se": se_i,
        "egger_intercept_p": p_two_sided_z(intercept / se_i) if se_i > 0 else None,
    }


def mr_group(rows: list[dict]) -> dict:
    rows = [r for r in rows if r["beta_exposure"] != 0 and r["se_outcome"] > 0]
    k = len(rows)
    if k == 0:
        return {}

    ratios = [r["beta_outcome_aligned"] / r["beta_exposure"] for r in rows]
    ratio_ses = [r["se_outcome"] / abs(r["beta_exposure"]) for r in rows]
    ratio_w = [1.0 / (s * s) for s in ratio_ses]

    base = {
        "n_instruments": k,
        "min_f": min(r["f_stat"] for r in rows),
        "mean_f": sum(r["f_stat"] for r in rows) / k,
        "weighted_median_beta": weighted_median(ratios, ratio_w),
    }

    if k == 1:
        b = ratios[0]
        se = ratio_ses[0]
        base.update({
            "primary_method": "Wald",
            "primary_beta": b,
            "primary_se": se,
            "primary_p": p_two_sided_z(b / se) if se > 0 else None,
            "ivw_beta": None,
            "ivw_se_fixed": None,
            "ivw_se_mre": None,
            "ivw_p_mre": None,
            "cochran_q": None,
            "cochran_q_df": None,
            "cochran_q_p": None,
        })
        base.update(wls_egger(rows))
        return base

    denom = 0.0
    numer = 0.0
    for r in rows:
        w = 1.0 / (r["se_outcome"] ** 2)
        bx = r["beta_exposure"]
        by = r["beta_outcome_aligned"]
        denom += w * bx * bx
        numer += w * bx * by

    if denom <= 0:
        return {}

    bivw = numer / denom
    se_fixed = math.sqrt(1.0 / denom)
    q = sum(
        (1.0 / (r["se_outcome"] ** 2))
        * (r["beta_outcome_aligned"] - bivw * r["beta_exposure"]) ** 2
        for r in rows
    )
    qdf = k - 1
    phi = max(1.0, q / qdf) if qdf > 0 else 1.0
    se_mre = se_fixed * math.sqrt(phi)
    p_mre = p_two_sided_z(bivw / se_mre) if se_mre > 0 else None

    base.update({
        "primary_method": "IVW_MRE",
        "primary_beta": bivw,
        "primary_se": se_mre,
        "primary_p": p_mre,
        "ivw_beta": bivw,
        "ivw_se_fixed": se_fixed,
        "ivw_se_mre": se_mre,
        "ivw_p_mre": p_mre,
        "cochran_q": q,
        "cochran_q_df": qdf,
        "cochran_q_p": chi2_sf(q, qdf),
    })
    base.update(wls_egger(rows))
    return base


def bh_adjust(pvals: list[Optional[float]]) -> list[Optional[float]]:
    idx = [(i, p) for i, p in enumerate(pvals) if p is not None and math.isfinite(p)]
    m = len(idx)
    out = [None] * len(pvals)
    if m == 0:
        return out
    idx.sort(key=lambda x: x[1])
    prev = 1.0
    for rank_rev, (i, p) in enumerate(reversed(idx), start=1):
        rank = m - rank_rev + 1
        q = min(prev, p * m / rank)
        prev = q
        out[i] = min(max(q, 0.0), 1.0)
    return out


def write_tsv(path: Path, rows: list[dict], fieldnames: Optional[list[str]] = None):
    if fieldnames is None:
        keys = []
        seen = set()
        for r in rows:
            for k in r:
                if k not in seen:
                    keys.append(k)
                    seen.add(k)
        fieldnames = keys
    with open_text(path, "wt") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/srv/is-analysis/IS_Analysis_V3")
    ap.add_argument(
        "--exposure",
        default="/srv/is-analysis/results/metabolic_resilience/stage1_pqtl/instruments/UKBPPP_CIS_INSTRUMENTS_PRIMARY_NON_MHC.tsv.gz",
    )
    ap.add_argument(
        "--out-dir",
        default="/srv/is-analysis/IS_Analysis_V3/results/skin_beauty/stage1_mr",
    )
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--retry-max", type=int, default=5)
    ap.add_argument("--pause", type=float, default=0.15)
    ap.add_argument("--p-threshold", type=float, default=5e-8)
    ap.add_argument(
        "--pilot-variants",
        type=int,
        default=0,
        help="If >0, deterministically sample this many unique query variants (rsID or chr:pos) instead of full query.",
    )
    ap.add_argument("--seed", type=int, default=20260927)
    ap.add_argument(
        "--exposure-only", action="store_true",
        help="Validate exposure schema/counts and stop before OpenGWAS calls.",
    )
    args = ap.parse_args()

    root = Path(args.root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    exposure = choose_exposure(args.exposure, root)
    print("===== STAGE 1A: EXPOSURE QC =====")
    print("exposure =", exposure)
    full_ex_rows, ex_qc = read_exposure(exposure, args.p_threshold)

    # Reproducibility gate is always evaluated on the FULL canonical exposure
    # before any pilot subsampling.
    full_unique_query_variants = sorted({r["query_variant"] for r in full_ex_rows})
    full_proteins = {r["protein_id"] for r in full_ex_rows}
    ex_qc["full_eligible_rows"] = len(full_ex_rows)
    ex_qc["full_eligible_proteins"] = len(full_proteins)
    ex_qc["full_unique_query_variants"] = len(full_unique_query_variants)

    # Always write query-key diagnostics before evaluating the gate.  This is
    # intentionally audit-first: a failed gate must still leave enough evidence
    # to explain exactly which canonical instruments were not queryable.
    no_query_key_rows = ex_qc.pop("no_query_key_rows", [])
    fallback_rows = [r for r in full_ex_rows if r["query_key_type"] == "chrpos"]
    write_tsv(out_dir / "STAGE1_CHRPOS_FALLBACK_INSTRUMENTS.tsv", fallback_rows)
    write_tsv(out_dir / "STAGE1_NO_QUERY_KEY_INSTRUMENTS.tsv", no_query_key_rows)
    ex_qc["no_query_key_audit_rows"] = len(no_query_key_rows)
    ex_qc["chrpos_fallback_audit_rows"] = len(fallback_rows)

    expected_gate = {}
    if exposure.name == "UKBPPP_CIS_INSTRUMENTS_PRIMARY_NON_MHC.tsv.gz":
        expected_gate = {
            "expected_rows": 10010,
            "observed_rows": len(full_ex_rows),
            "expected_proteins": 1936,
            "observed_proteins": len(full_proteins),
            "rsid_query_rows": ex_qc.get("rsid_query_rows", 0),
            "chrpos_query_rows": ex_qc.get("chrpos_query_rows", 0),
            "drop_no_query_key": ex_qc.get("drop_no_query_key", 0),
            "duplicate_protein_query_pairs": ex_qc.get("duplicate_protein_query_pairs", 0),
            "pass": (
                len(full_ex_rows) == 10010
                and len(full_proteins) == 1936
                and ex_qc.get("drop_no_query_key", 0) == 0
            ),
        }
        ex_qc["canonical_reproducibility_gate"] = expected_gate
        print("canonical_reproducibility_gate =", json.dumps(expected_gate))
        print(
            "query_key_audit =",
            json.dumps({
                "chrpos_fallback_rows": len(fallback_rows),
                "no_query_key_rows": len(no_query_key_rows),
                "no_query_key_file": str(out_dir / "STAGE1_NO_QUERY_KEY_INSTRUMENTS.tsv"),
            })
        )
        if not expected_gate["pass"]:
            (out_dir / "STAGE1_EXPOSURE_QC.json").write_text(
                json.dumps(ex_qc, indent=2, ensure_ascii=False), encoding="utf-8"
            )
            print("EXPOSURE_GATE=FAIL_CANONICAL_COUNT_MISMATCH")
            return 5

    # Preflight ends here: no JWT and no OpenGWAS allowance are required.
    if args.exposure_only:
        ex_qc["mode"] = "preflight"
        ex_qc["queried_unique_variants"] = 0
        (out_dir / "STAGE1_EXPOSURE_QC.json").write_text(
            json.dumps(ex_qc, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        write_tsv(out_dir / "STAGE1_ELIGIBLE_EXPOSURE.tsv.gz", full_ex_rows)
        print(
            f"mode=preflight eligible_rows={len(full_ex_rows)} "
            f"proteins={len(full_proteins)} unique_query_variants={len(full_unique_query_variants)} "
            f"rsid_rows={ex_qc.get('rsid_query_rows',0)} chrpos_rows={ex_qc.get('chrpos_query_rows',0)}"
        )
        print("EXPOSURE_GATE=PASS")
        return 0

    jwt = os.getenv("OPENGWAS_JWT", "").strip()
    if not jwt:
        print("ERROR: OPENGWAS_JWT is not configured.", file=sys.stderr)
        return 2

    # Pilot sampling happens only AFTER the full canonical gate has passed.
    unique_query_variants = list(full_unique_query_variants)
    ex_rows = list(full_ex_rows)
    if args.pilot_variants > 0 and args.pilot_variants < len(unique_query_variants):
        rng = random.Random(args.seed)
        unique_query_variants = sorted(rng.sample(unique_query_variants, args.pilot_variants))
        selected = set(unique_query_variants)
        ex_rows = [r for r in full_ex_rows if r["query_variant"] in selected]
        mode = "pilot"
    else:
        mode = "full"

    ex_qc["mode"] = mode
    ex_qc["queried_unique_variants"] = len(unique_query_variants)
    ex_qc["queried_exposure_rows"] = len(ex_rows)
    ex_qc["queried_proteins"] = len({r["protein_id"] for r in ex_rows})
    (out_dir / "STAGE1_EXPOSURE_QC.json").write_text(
        json.dumps(ex_qc, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    write_tsv(out_dir / "STAGE1_ELIGIBLE_EXPOSURE.tsv.gz", ex_rows)

    print(
        f"mode={mode} queried_rows={len(ex_rows)} "
        f"queried_proteins={len({r['protein_id'] for r in ex_rows})} "
        f"queried_unique_variants={len(unique_query_variants)}"
    )

    print("\n===== STAGE 1B: OPENGWAS ASSOCIATIONS =====")
    raw_rows, query_qc = query_batches(
        unique_query_variants, out_dir, jwt, args.batch_size, args.retry_max, args.pause
    )

    successful_batches = query_qc["n_batches"] - query_qc["n_failed_batches"]
    query_qc["successful_batches"] = successful_batches

    if query_qc["n_batches"] > 0 and successful_batches == 0:
        query_qc["query_gate"] = "FAIL_ALL_BATCHES"
        (out_dir / "STAGE1_QUERY_QC.json").write_text(
            json.dumps(query_qc, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        failure_summary = {
            "stage": "skin_beauty_stage1_mr",
            "mode": mode,
            "status": "QUERY_GATE_FAIL_ALL_BATCHES",
            "eligible_exposure_rows": len(ex_rows),
            "eligible_proteins": len({r["protein_id"] for r in ex_rows}),
            "unique_variants_requested": len(unique_query_variants),
            "failed_api_batches": query_qc["n_failed_batches"],
            "successful_api_batches": 0,
            "note": (
                "No outcome-overlap or MR statistics are interpretable "
                "because every OpenGWAS association batch failed."
            ),
        }
        (out_dir / "STAGE1_SUMMARY.json").write_text(
            json.dumps(failure_summary, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        print("\nSTAGE1_GATE=FAIL_ALL_API_BATCHES")
        return 3

    oy_rows = []
    for r in raw_rows:
        x = canonical_outcome_row(r)
        if x is not None:
            oy_rows.append(x)

    # Deduplicate exact returned records, retaining chr:pos identity for
    # instruments that were queried without rsIDs.
    oy_best = {}
    for r in oy_rows:
        key = (r["outcome_id"], r["rsid"], r.get("chrpos_outcome", ""), r["effect_allele_outcome"], r["other_allele_outcome"])
        oy_best.setdefault(key, r)
    oy_rows = list(oy_best.values())

    write_tsv(out_dir / "STAGE1_OUTCOME_ASSOCIATIONS.tsv.gz", oy_rows)

    overlap = {}
    requested_rsid = {v for v in unique_query_variants if re.fullmatch(r"rs\d+", v, flags=re.I)}
    requested_chrpos = set(unique_query_variants) - requested_rsid
    for oid, label in OUTCOMES.items():
        subset = [r for r in oy_rows if r["outcome_id"] == oid]
        returned_rsids = {r["rsid"] for r in subset if r["rsid"]}
        returned_chrpos = {r.get("chrpos_outcome", "") for r in subset if r.get("chrpos_outcome", "")}
        matched_rsid = requested_rsid & returned_rsids
        matched_chrpos = requested_chrpos & returned_chrpos
        overlap[oid] = {
            "analysis_label": label,
            "requested_variants": len(unique_query_variants),
            "requested_rsid_queries": len(requested_rsid),
            "requested_chrpos_queries": len(requested_chrpos),
            "matched_rsid_queries": len(matched_rsid),
            "matched_chrpos_queries": len(matched_chrpos),
            "matched_query_variants": len(matched_rsid) + len(matched_chrpos),
            "overlap_fraction": (len(matched_rsid) + len(matched_chrpos)) / len(unique_query_variants) if unique_query_variants else 0,
        }

    query_qc["canonical_association_rows"] = len(oy_rows)
    query_qc["overlap"] = overlap
    query_qc["query_gate"] = (
        "PASS"
        if query_qc["n_failed_batches"] == 0
        else "PARTIAL_RETRY_REQUIRED"
    )
    (out_dir / "STAGE1_QUERY_QC.json").write_text(
        json.dumps(query_qc, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    if query_qc["n_failed_batches"] > 0:
        print(
            f"WARNING: {query_qc['n_failed_batches']} API batches failed. "
            "Rerun the same command; successful batches are cached.",
            file=sys.stderr,
        )

    print("\n===== STAGE 1C: HARMONISATION =====")
    oy_by_rsid = defaultdict(list)
    oy_by_chrpos = defaultdict(list)
    for oy in oy_rows:
        if oy["rsid"]:
            oy_by_rsid[(oy["outcome_id"], oy["rsid"])].append(oy)
        if oy.get("chrpos_outcome"):
            oy_by_chrpos[(oy["outcome_id"], oy["chrpos_outcome"])].append(oy)

    harmonized = []
    harm_counts = defaultdict(int)
    for ex in ex_rows:
        for oid, label in OUTCOMES.items():
            if ex["query_key_type"] == "rsid":
                candidates = oy_by_rsid.get((oid, ex["query_variant"]), [])
            else:
                candidates = oy_by_chrpos.get((oid, ex["query_variant"]), [])

            if not candidates:
                harm_counts[f"{oid}:missing_outcome"] += 1
                continue

            chosen = None
            rejected = []
            for oy in candidates:
                ok, status, by, eaf_y_aligned = harmonize_one(ex, oy)
                if ok:
                    chosen = (oy, status, by, eaf_y_aligned)
                    break
                rejected.append(status)

            if chosen is None:
                status = rejected[0] if rejected else "no_allele_compatible_candidate"
                harm_counts[f"{oid}:{status}"] += 1
                continue

            oy, status, by, eaf_y_aligned = chosen
            harm_counts[f"{oid}:{status}"] += 1
            row = dict(ex)
            row.update(oy)
            row["harmonization_status"] = status
            row["beta_outcome_aligned"] = by
            row["eaf_outcome_aligned"] = eaf_y_aligned
            if ex["eaf_exposure"] is not None and eaf_y_aligned is not None:
                row["eaf_abs_diff"] = abs(ex["eaf_exposure"] - eaf_y_aligned)
            else:
                row["eaf_abs_diff"] = None
            harmonized.append(row)

    write_tsv(out_dir / "STAGE1_HARMONIZED.tsv.gz", harmonized)

    print("harmonized rows =", len(harmonized))

    print("\n===== STAGE 1D: PROTEIN-LEVEL MR =====")
    grouped = defaultdict(list)
    eligible_counts = defaultdict(int)
    for ex in ex_rows:
        for oid in OUTCOMES:
            eligible_counts[(oid, ex["protein_id"])] += 1
    for r in harmonized:
        grouped[(r["outcome_id"], r["protein_id"])].append(r)

    results = []
    for (oid, protein), rows in grouped.items():
        mr = mr_group(rows)
        if not mr:
            continue
        first = rows[0]
        rec = {
            "outcome_id": oid,
            "analysis_label": OUTCOMES[oid],
            "protein_id": protein,
            "gene_symbol": first.get("gene_symbol", ""),
            "uniprot": first.get("uniprot", ""),
            "n_eligible_exposure_instruments": eligible_counts[(oid, protein)],
            **mr,
        }
        rec["instrument_availability_fraction"] = (
            rec["n_instruments"] / rec["n_eligible_exposure_instruments"]
            if rec["n_eligible_exposure_instruments"] else None
        )
        rec["direction"] = (
            "higher_protein_higher_wrinkle_score"
            if rec["primary_beta"] > 0
            else "higher_protein_lower_wrinkle_score"
        )
        results.append(rec)

    # FDR within each outcome.
    for oid in OUTCOMES:
        idx = [i for i, r in enumerate(results) if r["outcome_id"] == oid]
        qs = bh_adjust([results[i]["primary_p"] for i in idx])
        m = len(idx)
        for j, i in enumerate(idx):
            results[i]["fdr_outcome"] = qs[j]
            results[i]["fdr_outcome_lt_0_05"] = (qs[j] is not None and qs[j] < 0.05)
            results[i]["bonferroni_outcome"] = (
                results[i]["primary_p"] < 0.05 / m
                if m and results[i]["primary_p"] is not None
                else False
            )

    global_q = bh_adjust([r["primary_p"] for r in results])
    for r, q in zip(results, global_q):
        r["fdr_global"] = q
        r["fdr_global_lt_0_05"] = (q is not None and q < 0.05)

    results.sort(
        key=lambda r: (
            r["outcome_id"],
            r["primary_p"] if r["primary_p"] is not None else 2.0,
        )
    )
    write_tsv(out_dir / "STAGE1_PROTEIN_MR.tsv.gz", results)

    candidates = [
        r for r in results
        if (r.get("fdr_outcome_lt_0_05") or r.get("bonferroni_outcome"))
    ]
    write_tsv(out_dir / "STAGE1_MR_CANDIDATES.tsv", candidates)

    summary = {
        "stage": "skin_beauty_stage1_mr",
        "mode": mode,
        "exposure": str(exposure),
        "eligible_exposure_rows": len(ex_rows),
        "eligible_proteins": len({r["protein_id"] for r in ex_rows}),
        "unique_variants_queried": len(unique_query_variants),
        "outcome_overlap": overlap,
        "harmonized_rows": len(harmonized),
        "harmonization_counts": dict(harm_counts),
        "proteins_tested": {
            oid: sum(r["outcome_id"] == oid for r in results) for oid in OUTCOMES
        },
        "fdr05_candidates": {
            oid: sum(
                r["outcome_id"] == oid and bool(r.get("fdr_outcome_lt_0_05"))
                for r in results
            )
            for oid in OUTCOMES
        },
        "bonferroni_candidates": {
            oid: sum(
                r["outcome_id"] == oid and bool(r.get("bonferroni_outcome"))
                for r in results
            )
            for oid in OUTCOMES
        },
        "failed_api_batches": query_qc["n_failed_batches"],
        "primary_harmonization_policy": (
            "exact/swap/non-palindromic strand harmonisation; "
            "all palindromic SNPs dropped for cross-ancestry primary MR; proxies=0"
        ),
    }
    (out_dir / "STAGE1_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\n===== STAGE 1 SUMMARY =====")
    print(json.dumps(summary, indent=2, ensure_ascii=False))

    if query_qc["n_failed_batches"] > 0:
        print("\nSTAGE1_GATE=PARTIAL_RETRY_REQUIRED")
        return 3
    if not harmonized:
        print("\nSTAGE1_GATE=FAIL_NO_HARMONIZED_VARIANTS")
        return 4

    print("\nSTAGE1_GATE=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
