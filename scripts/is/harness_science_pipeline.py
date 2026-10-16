#!/usr/bin/env python3
"""IS science-stage gate: audited P1 provenance, independent P2 ID/AF checks,
P3 fail-closed coloc eligibility and P4 gene evidence joins.

No source TSV mutation, arbitrary scripts, inferred cohort ancestry, or
automatic colocalization on unmatched LD. Python standard library only.
"""
from __future__ import annotations
import collections
import csv
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile

TRUTH = {"true", "1", "yes"}
NUCLEOTIDES = {"A", "C", "G", "T"}
SNP_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT):([0-9]+):([ACGT]+):([ACGT]+)$", re.I)
ROOT = Path("/srv/is-analysis/results/is")
OUT = ROOT / "github_actions_runs" / "evidence"
INPUTS = OUT / "science_inputs"
P1_REGISTRY = INPUTS / "P1_APPROVED_QTL_LD.tsv"
P2_REF = INPUTS / "P2_APPROVED_GRCH37_REFERENCE.json"
P3_REGISTRY = INPUTS / "P3_APPROVED_COLOC_INPUTS.tsv"
P1_COLUMNS = (
    "region_id", "QTL_dataset", "gene_id", "qtl_source_cohort",
    "source_ancestry", "ld_source_cohort", "ld_ancestry", "ld_matrix_path",
    "ld_sha256", "verified_by", "verified_at", "evidence_url"
)
P3_COLUMNS = (
    "region_id", "QTL_dataset", "gene_id", "input_path", "input_sha256",
    "qtl_sample_size", "gwas_sample_size", "build", "variant_alleles_verified",
    "matched_ld_matrix_verified", "reviewed_by", "reviewed_at"
)


def stamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def safe_path(path_string: str, *, root: Path) -> Path | None:
    """Enforce absolute, existing nonsymlink regular files under configured root."""
    try:
        candidate = Path(path_string)
        if not candidate.is_absolute() or not candidate.is_file():
            return None
        resolved = candidate.resolve(strict=True)
        if not resolved.is_relative_to(root.resolve(strict=True)):
            return None
        if candidate.is_symlink():
            return None
        return resolved
    except (OSError, ValueError):
        return None


def file_sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def registry_rows(path: Path, columns: tuple[str, ...]):
    if not path.is_file():
        return [], "MISSING_REGISTRY"
    if path.stat().st_size > 10_000_000:
        raise ValueError("Input registry exceeds safe limit")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Registry header is missing or duplicated")
        missing = set(columns).difference(reader.fieldnames)
        if missing:
            raise ValueError("Missing registry columns: " + ",".join(sorted(missing)))
        rows = list(reader)
        if any(None in r or None in r.values() for r in rows):
            raise ValueError("Malformed registry row")
    return rows, "PRESENT"


def verified_p1_registry(records, *, root=ROOT):
    """A verified LD cohort demands explicit provenance AND actual file hash."""
    approved = {}
    decisions = collections.Counter()
    for row in records:
        key = (row["region_id"], row["QTL_dataset"], row["gene_id"])
        if key in approved:
            raise ValueError("Duplicate approved P1 key")
        cohort = row["qtl_source_cohort"].strip()
        ancestry = row["source_ancestry"].strip()
        if not (cohort and ancestry and row["verified_by"].strip() and row["verified_at"].strip()
                and row["evidence_url"].strip()):
            decisions["MISSING_COHORT_PROVENANCE"] += 1
            continue
        if cohort != row["ld_source_cohort"].strip() or ancestry != row["ld_ancestry"].strip():
            decisions["LD_COHORT_OR_ANCESTRY_MISMATCH"] += 1
            continue
        p = safe_path(row["ld_matrix_path"], root=root)
        if p is None or not re.fullmatch(r"[a-f0-9]{64}", row["ld_sha256"]):
            decisions["LD_FILE_OR_HASH_MISSING"] += 1
            continue
        if file_sha256(p) != row["ld_sha256"]:
            decisions["LD_SHA256_MISMATCH"] += 1
            continue
        # Passing this is only *provenance verified*. Row/column LD integrity and
        # allele order remain independent P3 gates.
        approved[key] = {"ld_sha256":row["ld_sha256"],"cohort":cohort,
                         "ancestry":ancestry,"path":str(p)}
        decisions["PROVENANCE_REVIEWED_LD_HASH_MATCH"] += 1
    return approved, dict(decisions)


def parse_variant(value):
    value = str(value or "").strip()
    match = SNP_RE.fullmatch(value)
    if not match:
        return None
    chrom, pos, ref, alt = match.groups()
    chrom = chrom.upper().removeprefix("CHR")
    if chrom.isdigit():
        chrom = str(int(chrom))
    ref = ref.upper()
    alt = alt.upper()
    if int(pos) < 1 or ref == alt or not set(ref+alt).issubset(NUCLEOTIDES):
        return None
    return chrom, int(pos), ref, alt


def allele_pair_equal(left, right):
    return left is not None and right is not None and left == right


def palindromic(parsed):
    return bool(parsed and len(parsed[2]) == 1 and len(parsed[3]) == 1
        and frozenset((parsed[2], parsed[3])) in
        (frozenset(("A", "T")), frozenset(("C", "G"))))


def numeric_af(value):
    if value is None or str(value).strip() in {"", "NA", "nan", "None", "."}:
        return None
    try:
        v = float(value)
    except (ValueError, TypeError):
        return None
    return v if math.isfinite(v) and 0 <= v <= 1 else None


class ApprovedFasta:
    """Indexed GRCh37 FASTA, read-only; requires manually approved source manifest."""

    def __init__(self, manifest_path=P2_REF, allowed_root=Path("/srv/is-analysis")):
        self.status = "BLOCKED_NO_APPROVED_GRCH37_FASTA"
        self.ref = None
        self.fai = {}
        if not manifest_path.is_file():
            return
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
        if data.get("assembly") != "GRCh37" or not (data.get("approved_by") and data.get("reference_source")):
            self.status = "BLOCKED_REFERENCE_APPROVAL_INCOMPLETE"
            return
        path = safe_path(str(data.get("fasta_path", "")), root=allowed_root)
        index = safe_path(str(data.get("fai_path", "")), root=allowed_root)
        if not path or not index or index != Path(str(path)+".fai"):
            self.status = "BLOCKED_REFERENCE_OR_FAI_MISSING"
            return
        with index.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise ValueError("Malformed FASTA .fai index")
                chrom = fields[0].upper().removeprefix("CHR")
                if chrom.isdigit():
                    chrom = str(int(chrom))
                self.fai[chrom] = tuple(int(v) for v in fields[1:5])
        self.ref = path
        self.status = "APPROVED_REFERENCE_AVAILABLE"

    def reference(self, chrom, pos, size):
        if self.ref is None or chrom not in self.fai or not 1 <= size <= 1000:
            return None
        length, offset, bases, width = self.fai[chrom]
        if not 1 <= pos <= length or pos + size - 1 > length or bases <= 0 or width < bases:
            return None
        start = pos - 1
        try:
            with self.ref.open("rb") as handle:
                result = bytearray()
                for i in range(size):
                    absolute = start + i
                    handle.seek(offset + (absolute//bases)*width + absolute%bases)
                    result.extend(handle.read(1))
                bases_read = result.decode("ascii").upper()
                return bases_read if set(bases_read).issubset(NUCLEOTIDES) else None
        except (OSError, UnicodeError):
            return None


def independent_p2_qc(variants_path: Path, fasta=None, max_rows=2_000_000):
    """Independent string-position/REF/ALT and EAF consistency checks.

    SNP matches, AF deltas and approved-FASTA REF comparisons are descriptive;
    they cannot independently establish QTL ancestry, genome liftOver correctness,
    genomic build, signed-LD validity or colocalization.
    """
    fa = fasta if fasta is not None else ApprovedFasta()
    counters = collections.Counter()
    grouped = collections.defaultdict(collections.Counter)
    seen_fasta = {}
    by_gene = collections.defaultdict(collections.Counter)
    with variants_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"region_id","QTL_dataset","gene_id",
            "mapped_variant_GRCh37","BBJ_exact_SNP","EAS_AIS_exact_SNP",
            "QTL_alt_eaf_ac_div_an","BBJ_alt_eaf","EAS_AIS_alt_eaf"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("P2 input missing required columns: "+str(sorted(required-set(reader.fieldnames or []))))
        for row in reader:
            if None in row or None in row.values():
                raise ValueError("P2 malformed row")
            counters["total_rows"] += 1
            if counters["total_rows"] > max_rows:
                raise RuntimeError("P2 row cap exceeded")
            key = (row["region_id"], row["QTL_dataset"], row["gene_id"])
            group = grouped[key]
            group["rows"] += 1
            mapped = parse_variant(row["mapped_variant_GRCh37"])
            bbj = parse_variant(row["BBJ_exact_SNP"])
            eas = parse_variant(row["EAS_AIS_exact_SNP"])
            if mapped:
                counters["mapped_variant_parseable"] += 1
            else:
                counters["mapped_variant_unparseable"] += 1
            if allele_pair_equal(mapped, bbj):
                group["bbj_allele_key_equal"] += 1
                counters["bbj_allele_key_equal"] += 1
            if allele_pair_equal(mapped, eas):
                group["eas_allele_key_equal"] += 1
                counters["eas_allele_key_equal"] += 1
            dual_match = allele_pair_equal(mapped, bbj) and allele_pair_equal(mapped, eas)
            if dual_match:
                group["dual_allele_key_equal"] += 1
                counters["dual_allele_key_equal"] += 1
                by_gene[row["gene_id"]]["dual_allele_key_equal"] += 1
            if palindromic(mapped):
                counters["palindromic_rows"] += 1
            else:
                if dual_match:
                    group["nonpal_dual_allele_key_equal"] += 1
                    by_gene[row["gene_id"]]["nonpal_dual_allele_key_equal"] += 1
            q_af = numeric_af(row["QTL_alt_eaf_ac_div_an"])
            e_af = numeric_af(row["EAS_AIS_alt_eaf"])
            b_af = numeric_af(row["BBJ_alt_eaf"])
            if q_af is not None and e_af is not None:
                counters["qtl_eas_af_comparable"] += 1
                if abs(q_af-e_af)>0.2:
                    counters["qtl_eas_af_delta_gt_0_2_review"] += 1
                    group["qtl_eas_af_delta_gt_0_2_review"] += 1
                elif dual_match:
                    group["allele_key_equal_and_qtl_eas_af_delta_le_0_2"] += 1
            if b_af is not None and e_af is not None:
                counters["bbj_eas_af_comparable"] += 1
                if abs(b_af-e_af)>0.2:
                    counters["bbj_eas_af_delta_gt_0_2_review"] += 1
            if fa.ref is not None and mapped:
                loc = mapped[:2]
                if loc not in seen_fasta:
                    seen_fasta[loc] = fa.reference(*loc, len(mapped[2]))
                if seen_fasta[loc] is None:
                    counters["fasta_ref_unavailable"] += 1
                elif mapped[2] == seen_fasta[loc]:
                    counters["fasta_ref_equal_rows"] += 1
                    if dual_match and not palindromic(mapped):
                        group["fasta_and_nonpal_allele_key_equal"] += 1
                else:
                    counters["fasta_ref_mismatch_rows"] += 1
            elif mapped:
                counters["fasta_not_checked_rows"] += 1
    if not counters["total_rows"]:
        raise ValueError("P2 no SNP observations")
    return {
        "schema":"is-independent-p2-v1","generated_at":stamp(),"status":"AUDITED_WITH_GATES",
        "fasta_status":fa.status,"counters":dict(counters),
        "groups":[{"region_id":k[0],"QTL_dataset":k[1],"gene_id":k[2],**dict(v)}
                  for k,v in sorted(grouped.items())],
        "by_gene":{k:dict(v) for k,v in sorted(by_gene.items())},
        "caveats":[
            "String allele equality does not independently validate liftOver/reference assembly.",
            "EAF delta >0.2 is a review flag, not statistical evidence of strand inversion.",
            "Unapproved or missing FASTA means genome reference validation is blocked.",
        ],
    }


def evaluate_p1_gate(readiness_path: Path, approved):
    rows, summaries = [], collections.Counter()
    with readiness_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required={"region_id","QTL_dataset","gene_id"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("P1 readiness schema mismatch")
        for r in reader:
            if None in r or None in r.values():
                raise ValueError("P1 malformed row")
            key=(r["region_id"],r["QTL_dataset"],r["gene_id"])
            approved_ld=key in approved
            condition="PROVENANCE_REVIEWED_LD_HASH_MATCH" if approved_ld else "BLOCKED_NO_VERIFIED_COHORT_MATCHED_LD"
            summaries[condition]+=1
            rows.append({"region_id":key[0],"QTL_dataset":key[1],"gene_id":key[2],
                         "p1_status":condition})
    if not rows:
        raise ValueError("P1 empty readiness")
    return {"schema":"is-p1-source-gates-v1","generated_at":stamp(),"counts":dict(summaries),"rows":rows}


def evaluate_p3_gate(p1, p2, approved_ld, p3_inputs, p2_fasta_approved=False):
    """Never invoke coloc without individually reviewed, versioned full inputs.

    LD file hash, P2 FASTA REF match and reviewed analysis input metadata are
    prerequisites. R/coloc execution deliberately remains disabled until the
    sample-size/LD matrix/variant-order validation runner is independently tested.
    """
    records=[]
    counts=collections.Counter()
    p2group={(g["region_id"],g["QTL_dataset"],g["gene_id"]):g for g in p2["groups"]}
    registry={ (r["region_id"],r["QTL_dataset"],r["gene_id"]):r for r in p3_inputs }
    for row in p1["rows"]:
        key=(row["region_id"],row["QTL_dataset"],row["gene_id"])
        group=p2group.get(key,{})
        reasons=[]
        if key not in approved_ld:
            reasons.append("BLOCKED_COHORT_MATCHED_LD")
        if not p2_fasta_approved or group.get("fasta_and_nonpal_allele_key_equal",0)<1:
            reasons.append("BLOCKED_INDEPENDENT_GRCH37_FASTA_ALLELE")
        candidate=registry.get(key)
        if not candidate:
            reasons.append("BLOCKED_MISSING_REVIEWED_COLOC_INPUT")
        else:
            if (candidate["build"]!="GRCh37" or
                candidate["variant_alleles_verified"].strip().lower() not in TRUTH or
                candidate["matched_ld_matrix_verified"].strip().lower() not in TRUTH or
                not candidate["reviewed_by"].strip() or not candidate["reviewed_at"].strip()):
                reasons.append("BLOCKED_COLOC_INPUT_REVIEW")
            path=safe_path(candidate["input_path"],root=ROOT)
            if (not path or not re.fullmatch(r"[a-f0-9]{64}",candidate["input_sha256"])
                    or (path and file_sha256(path)!=candidate["input_sha256"])):
                reasons.append("BLOCKED_COLOC_INPUT_HASH")
            for name in ("qtl_sample_size","gwas_sample_size"):
                try:
                    if int(candidate[name])<=0:reasons.append("BLOCKED_SAMPLE_SIZE")
                except (ValueError,TypeError):
                    reasons.append("BLOCKED_SAMPLE_SIZE")
        # Even when prereqs are met the new coloc executor still needs verified
        # same-allele signed LD matrix and susie/coloc R implementation tests.
        if not reasons:
            reasons.append("AWAITING_VALIDATED_COLOC_EXECUTOR")
        status=";".join(sorted(set(reasons)))
        for reason in set(reasons):
            counts[reason]+=1
        records.append({"region_id":key[0],"QTL_dataset":key[1],"gene_id":key[2],
                        "status":status,"coloc_executed":False})
    return {"schema":"is-p3-coloc-execution-gate-v1","generated_at":stamp(),
            "counts":dict(counts),"rows":records,"coloc_executed":0,
            "disclaimer":"Not a completed coloc analysis. P3 R executor is gated, not simulated."}


def join_p4_candidates(candidate_detail, p2, p1, p3):
    p2_by_gene=collections.defaultdict(collections.Counter)
    for r in p2["groups"]:
        z=p2_by_gene[r["gene_id"]]
        z["independently_parsed_allele_rows"]+=r.get("dual_allele_key_equal",0)
        z["nonpal_dual_key_rows"]+=r.get("nonpal_dual_allele_key_equal",0)
        z["fasta_validated_nonpal_rows"]+=r.get("fasta_and_nonpal_allele_key_equal",0)
    p1_by_gene=collections.defaultdict(collections.Counter)
    for r in p1["rows"]:
        p1_by_gene[r["gene_id"]][r["p1_status"]]+=1
    p3_by_gene=collections.defaultdict(collections.Counter)
    for r in p3["rows"]:
        p3_by_gene[r["gene_id"]]["assayed_groups"]+=1
    results=[]
    status_counts=collections.Counter()
    for candidate in candidate_detail["candidates"]:
        gene=candidate["gene_id"]
        e=p2_by_gene.get(gene,{})
        label = ("HAS_INDEPENDENT_SNP_ID_OVERLAP" if e.get("nonpal_dual_key_rows",0)>0
                 else "NO_CONFIRMED_NONPAL_SNP_ID_OVERLAP")
        status_counts[label]+=1
        results.append({
            "gene_id":gene,"gene_symbol":candidate["gene_symbol"],
            "source_nominal_p":candidate["nominal_min_p"],
            "nominal_variant_rows":candidate["nominal_variant_rows"],
            "dual_allele_key_equal":e.get("independently_parsed_allele_rows",0),
            "nonpal_dual_allele_key_equal":e.get("nonpal_dual_key_rows",0),
            "fasta_validated_nonpal_rows":e.get("fasta_validated_nonpal_rows",0),
            "p1_registry_matched_ld":p1_by_gene.get(gene,{}).get("PROVENANCE_REVIEWED_LD_HASH_MATCH",0),
            "p3_coloc_executed":False,"evidence_tier":label,
        })
    return {"schema":"is-p4-cross-evidence-v1","generated_at":stamp(),
            "candidate_count":len(results),"tier_counts":dict(status_counts),
            "candidates":results,
            "interpretation":"Independent SNP identifier QC crossed with nominal-QTL candidates; no causal inference."}


def atomic_json(path, obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile("w",encoding="utf-8",dir=path.parent,prefix=".science_",delete=False) as out:
        json.dump(obj,out,ensure_ascii=False,indent=2,sort_keys=True)
        out.write("\n")
        temp=Path(out.name)
    try:
        os.replace(temp,path)
    finally:
        temp.unlink(missing_ok=True)


def advance(variants_path, readiness_path, candidate_detail, output=OUT, *, registry=P1_REGISTRY, reference=P2_REF,
            p3_registry=P3_REGISTRY, root=ROOT):
    p1_rows,p1_state=registry_rows(registry,P1_COLUMNS)
    approved,provenance_decisions=verified_p1_registry(p1_rows,root=root)
    p1=evaluate_p1_gate(readiness_path,approved)
    fasta=ApprovedFasta(reference)
    p2=independent_p2_qc(variants_path,fasta=fasta)
    p3_rows,p3_state=registry_rows(p3_registry,P3_COLUMNS)
    p3=evaluate_p3_gate(p1,p2,approved,p3_rows,p2_fasta_approved=fasta.ref is not None)
    p4=join_p4_candidates(candidate_detail,p2,p1,p3)
    report={
        "schema":"is-science-orchestrator-v1","generated_at":stamp(),"integrity":"PASS",
        "source_provenance_state":p1_state,"p1_registry_decisions":provenance_decisions,
        "p2_reference_state":fasta.status,"p3_input_registry_state":p3_state,
        "p1": {"assayed_groups":len(p1["rows"]),"approved_matched_ld":len(approved),
               "blocked":p1["counts"].get("BLOCKED_NO_VERIFIED_COHORT_MATCHED_LD",0)},
        "p2": {"scanned_snp_records":p2["counters"]["total_rows"],
               "dual_id_equal":p2["counters"].get("dual_allele_key_equal",0),
               "nonpalindromic_records":p2["counters"].get("dual_allele_key_equal",0)-p2["counters"].get("palindromic_rows",0),
               "fasta_ref_equal_rows":p2["counters"].get("fasta_ref_equal_rows",0),
               "fasta_ref_mismatch_rows":p2["counters"].get("fasta_ref_mismatch_rows",0),
               "af_review_flags":p2["counters"].get("qtl_eas_af_delta_gt_0_2_review",0)},
        "p3": {"coloc_executed":0,"blocked_gates":p3["counts"]},
        "p4": {"candidate_count":p4["candidate_count"],"tier_counts":p4["tier_counts"]},
        "scientific_result":"NO_NEW_COLOC_OR_CAUSAL_GENE_EVIDENCE",
    }
    # Write detailed individual-level files only to server-local evidence folder.
    # Scientific inputs themselves are never mutated.
    science=output/"science"
    for name,value in [
        ("p1_cohort_ld_review.json",p1),("p2_independent_variant_qc.json",p2),
        ("p3_coloc_gates.json",p3),("p4_candidate_cross_evidence.json",p4),
        ("science_pipeline_summary.json",report)
    ]:
        atomic_json(science/name,value)
    return report
