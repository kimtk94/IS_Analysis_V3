#!/usr/bin/env python3
"""Validate OASIS source SNP identity, chromosome-window coverage, and AIS allele schema.

Strictly DOES NOT harmonize signed OASIS beta to AIS until source beta allele
semantics are proven. Existing rsID/liftover map is a prior candidate reference,
NOT independently newly verified remapping.
"""
import argparse
import csv
import hashlib
import json
from math import isfinite
from pathlib import Path

from audit_g0022_oasis_browser_qtl import parse

CS = (
    ("rs11066015", "12:112168009:G:A", "chr12:111730205:G:A", 0.2724730523),
    ("rs671", "12:112241766:G:A", "chr12:111803962:G:A", 0.3748385486),
    ("rs11066132", "12:112468206:C:T", "chr12:112030402:C:T", 0.2977596915),
    ("rs77768175", "12:112736118:A:G", "chr12:112298314:A:G", 0.03778),
)
GENES = ("ALDH2", "BRAP", "RPH3A")
CELLS = ("Mono-L1", "B_Activated-L2")
AIS_STUDY = "GCST90104545"
AIS_TRAIT = "EAS_AIS"


def split_variant(key, style):
    parts = key.split(":")
    if len(parts) != 4:
        raise ValueError(f"Variant ID must have four fields: {key}")
    chrom, pos, ref, alt = parts
    chrom = chrom.removeprefix("chr")
    if not chrom.isdigit() or not pos.isdigit() or int(pos) <= 0:
        raise ValueError(f"Invalid coordinates: {key}")
    if not ref or not alt or ref == alt or not all(c in "ACGT" for c in (ref + alt)):
        raise ValueError(f"Invalid allele representation: {key}")
    return chrom, int(pos), ref, alt


def allele_set_result(gwas_id, oasis_id, expected_hg38):
    gc, gp, gref, galt = split_variant(gwas_id, "GRCh37")
    oc, op, oref, oalt = split_variant(oasis_id, "GRCh38")
    ec, ep, eref, ealt = split_variant(expected_hg38, "GRCh38")
    if (oc, op) != (ec, ep) or oc != gc:
        return "MISMATCH_BUILD_OR_POSITION"
    if (oref, oalt) == (gref, galt) == (eref, ealt):
        return "EXACT_REF_ALT_STRING_MATCH_AFTER_PRIOR_LIFTOVER"
    if (oref, oalt) == (galt, gref) and (eref, ealt) == (gref, galt):
        return "SWAPPED_REF_ALT_REQUIRES_SIGN_FLIP_IF_ALT_EFFECT_VERIFIED"
    if {oref, oalt} == {gref, galt} == {eref, ealt}:
        return "ALLELE_SET_MATCH_BUT_UNRESOLVED_DIRECTION"
    return "ALLELE_MISMATCH_DO_NOT_HARMONIZE"


def load_ais_four(path):
    records = {}
    with Path(path).open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["study"] != AIS_STUDY or row["trait"] != AIS_TRAIT:
                continue
            key = row["variant"]
            if key in records:
                raise ValueError(f"Duplicate AIS variant {key}")
            if row["ref"] != key.split(":")[2]:
                raise ValueError("GWAS ref inconsistently named")
            original_allele = row["original_effect_allele"]
            ref, alt = key.split(":")[2:]
            if row["effect_allele"] != "ALT" or original_allele not in (ref, alt):
                raise ValueError(f"GWAS ALT alignment not proven for {key}")
            # ALT_beta was already source-audited against the original source
            # effect_allele and flipped as necessary. Never require that the
            # ORIGINAL source effect allele itself is ALT; for HECTD4 rs77768175,
            # source used REF A, while audited ALT is G.
            beta, se, eaf = (float(row["ALT_beta"]), float(row["se"]),
                             float(row["ALT_EAF"]))
            if not all(isfinite(x) for x in (beta, se, eaf)) or se <= 0 or not 0 <= eaf <= 1:
                raise ValueError("Invalid AIS effect, SE, EAF")
            records[key] = {"study": AIS_STUDY, "trait": AIS_TRAIT,
                            "gwas_beta_alt": beta, "gwas_se_alt": se,
                            "gwas_alt_eaf": eaf, "gwas_p": float(row["p"])}
    expected = {x[1] for x in CS}
    if set(records) != expected:
        raise ValueError(f"Wrong AIS four-variant table. Missing {expected-set(records)}")
    return records


def audit_page(gene, cell, content, ais):
    text = content.decode("utf-8", "replace")
    try:
        rows = parse(text, gene)
        state = "BOKEH_PLOTTED_VARIANTS_AVAILABLE"
    except ValueError as e:
        if "Missing embedded Bokeh data" not in str(e):
            raise
        rows = []
        state = "NO_BOKEH_DATA_NOT_A_NEGATIVE_ASSOCIATION"
    pos = [int(r["variant_id_hg38"].split(":")[1]) for r in rows
           if r["variant_id_hg38"].startswith("chr12:")]
    if rows and len(pos) != len(rows):
        raise ValueError(f"Unexpected chromosome or variant schema in {gene}/{cell}")
    chr_range = (min(pos), max(pos)) if pos else (None, None)
    audit = {"gene": gene, "celltype": cell, "source_status": state,
             "bokeh_plotted_snp_count": len(rows),
             "non_significant_p_greater_than_0_05": sum(float(x["pval_nominal"]) > .05 for x in rows),
             "observed_chr12_hg38_min": chr_range[0], "observed_chr12_hg38_max": chr_range[1],
             "source_complete_all_tested_cis_verified": False}
    indexed = {}
    for row in rows:
        indexed.setdefault(row["rsid"], []).append(row)
    annotations = []
    for rsid, gwas_id, expected_hg38, pip in CS:
        c, p38, ref38, alt38 = split_variant(expected_hg38, "GRCh38")
        nearby = min((abs(p - p38) for p in pos), default=None)
        hits = indexed.get(rsid, [])
        if len(hits) > 1:
            raise ValueError(f"Duplicate QTL rsID in gene-cell page {gene}/{cell}/{rsid}")
        status = "NOT_ASSESSED_NO_SOURCE"
        allele = None
        row = hits[0] if hits else None
        if rows and not hits:
            status = "NOT_PRESENT_IN_BOKEH_RENDERED_PAGE_NOT_BIOLOGICAL_NEGATIVE"
        if row is not None:
            allele = allele_set_result(gwas_id, row["variant_id_hg38"], expected_hg38)
            if allele != "EXACT_REF_ALT_STRING_MATCH_AFTER_PRIOR_LIFTOVER":
                status = "BLOCKED_ALLELE_IDENTITY_MISMATCH"
            else:
                status = "EXACT_RSID_POSITION_AND_REF_ALT_MATCH"
        rec = {"gene": gene, "celltype": cell, "rsid": rsid,
               "gwas_variant_grch37": gwas_id, "prior_rsid_mapped_grch38": expected_hg38,
               "ais_exploratory_cs_pip": pip,
               "within_observed_bokeh_min_max_range": chr_range[0] <= p38 <= chr_range[1] if pos else None,
               "nearest_bokeh_snp_distance_bp": nearby,
               "source_coverage_status": status,
               "allele_identity_check": allele,
               "actual_qtl_variant_hg38": row["variant_id_hg38"] if row else None,
               "oasis_qtl_beta_not_harmonized": float(row["effect_size"]) if row else None,
               "oasis_qtl_se": float(row["effect_size_SE"]) if row else None,
               "oasis_qtl_p": float(row["pval_nominal"]) if row else None,
               "oasis_qtl_maf": float(row["maf"]) if row else None,
               "gwas_beta_alt": ais[gwas_id]["gwas_beta_alt"],
               "gwas_se_alt": ais[gwas_id]["gwas_se_alt"],
               "gwas_alt_eaf": ais[gwas_id]["gwas_alt_eaf"],
               "qtl_effect_allele_author_verified": False,
               "effect_sign_harmonization_status": "BLOCKED_QTL_EFFECT_ALLELE_SEMANTICS_UNVERIFIED",
               "qtl_gwas_colocalization_valid": False,
               "causal_gene_validated": False}
        if row is not None:
            pval, std, maf = rec["oasis_qtl_p"], rec["oasis_qtl_se"], rec["oasis_qtl_maf"]
            if not (0 <= pval <= 1 and std > 0 and 0 <= maf <= .5):
                raise ValueError("OASIS malformed statistical fields")
        annotations.append(rec)
    audit["source_identified_cs_markers"] = sum(
        r["source_coverage_status"] == "EXACT_RSID_POSITION_AND_REF_ALT_MATCH" for r in annotations)
    audit["unobserved_bokeh_markers"] = sum(
        r["source_coverage_status"] == "NOT_PRESENT_IN_BOKEH_RENDERED_PAGE_NOT_BIOLOGICAL_NEGATIVE"
        for r in annotations)
    return audit, annotations


def write_tsv(path, rows):
    if not rows:
        raise ValueError("No output")
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run(args):
    ais = load_ais_four(args.gwas)
    manifest, records = [], []
    for gene in GENES:
        for cell in CELLS:
            path = args.browser_dir / f"OASIS_{gene}_{cell}_original.html"
            content = path.read_bytes()
            audit, rows = audit_page(gene, cell, content, ais)
            audit["source_file"] = str(path)
            audit["sha256"] = hashlib.sha256(content).hexdigest()
            manifest.append(audit)
            records.extend(rows)
    summary = {
        "audit": "IS_G0022_OASIS_AIS_SOURCE_COVERAGE_AND_ALLELE_IDENTITY_V1",
        "gwas_study": AIS_STUDY, "gwas_trait": AIS_TRAIT,
        "gwas_provenance": str(args.gwas),
        "gwas_sha256": hashlib.sha256(args.gwas.read_bytes()).hexdigest(),
        "browser_sources": manifest,
        "sources_reviewed": len(manifest),
        "sc_gene_cell_pages_with_bokeh": sum(x["bokeh_plotted_snp_count"] > 0 for x in manifest),
        "plotted_rs11066015_gene_cell_matches": sum(x["source_coverage_status"] == "EXACT_RSID_POSITION_AND_REF_ALT_MATCH" for x in records),
        "other_3_cs_snp_test_status": "NOT_PRESENT_IN_SOURCE_BROWSER_SUBSET_NOT_RAW_GENOTYPE_NEGATIVE",
        "gwas_rs11066015_alt_beta": ais["12:112168009:G:A"]["gwas_beta_alt"],
        "reference_build_mapping_coverage": "PREVIOUS_SOURCE_MAPPING_NOT_NEW_INDEPENDENT_LIFTOVER",
        "effect_allele_harmonization_performed": False,
        "valid_colocalizations": 0,
        "causal_gene_confirmed": False,
        "scientific_claim": "PARTIAL_OASIS_MARGINAL_QTL_WITH_REFERENCE_IDENTITY_NOT_CAUSAL",
        "caveat": ("Only displayed Bokeh SNP rows, not independently verified complete tested cis universe. "
                   "Three AIS-CS variants absent from plot; no evidence of genotype absence. "
                   "MAF is not equivalent to oriented GWAS ALT EAF; QTL beta effect-allele convention unknown. "
                   "No valid GWAS-QTL coloc, MR or causal transcript mediation.")
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir / "G0022_OASIS_FOUR_CS_ALLELE_COVERAGE.tsv", records)
    write_tsv(args.outdir / "G0022_OASIS_GENE_CELL_SOURCE_COVERAGE.tsv", manifest)
    (args.outdir / "G0022_OASIS_ALLELE_COVERAGE_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False)+"\n")
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--browser-dir", required=True, type=Path)
    p.add_argument("--gwas", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args=p.parse_args()
    print(json.dumps(run(args),indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
