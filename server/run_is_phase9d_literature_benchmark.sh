#!/usr/bin/env bash

ROOT="/srv/is-analysis"
P9C="$ROOT/results/is/stage5_functional/phase9c_convergence"

OUT="$ROOT/results/is/stage5_functional/phase9d_literature_benchmark"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_phase9d_literature/$RUN_ID"
STATUS="$OUT/STEP_STATUS.tsv"

mkdir -p "$OUT" "$LOGDIR"

printf "step\tstatus\trc\n" > "$STATUS"

run_step() {
  local name="$1"
  shift

  local LOG="$LOGDIR/${name}.log"

  echo
  echo "================================================================================"
  echo "$name"
  echo "================================================================================"

  (
    "$@"
  ) > >(tee "$LOG") 2>&1

  local RC=$?

  if [[ $RC -eq 0 ]]; then
    printf "%s\tPASS\t%s\n" "$name" "$RC" >> "$STATUS"
  else
    printf "%s\tFAIL\t%s\n" "$name" "$RC" >> "$STATUS"
  fi

  return 0
}


# =============================================================================
# 01. PREFLIGHT
# =============================================================================

preflight() {

  echo "===== PHASE9C ====="

  cat "$P9C/STEP_STATUS.tsv"

  FAIL_N="$(
    awk -F'\t' '
      NR>1 && $2!="PASS" {n++}
      END {print n+0}
    ' "$P9C/STEP_STATUS.tsv"
  )"

  echo "PHASE9C_FAIL_N=$FAIL_N"

  if [[ "$FAIL_N" -ne 0 ]]; then
    return 1
  fi

  for f in \
    "$P9C/COLOC_ABF_MASTER_ANNOTATED_V2.tsv" \
    "$P9C/COLOC_SUSIE_SIGNAL_PAIRS.tsv" \
    "$P9C/RS671_FUNCTIONAL_ANNOTATION.tsv" \
    "$P9C/L003_STATISTICAL_FREEZE.tsv"
  do

    if [[ -s "$f" ]]; then
      echo "FOUND $f"
    else
      echo "MISSING $f"
      return 1
    fi

  done

  echo "PREFLIGHT=PASS"
}


# =============================================================================
# 02. REPAIR / EXPAND FUNCTIONAL CONVERGENCE MASTER
# =============================================================================

repair_convergence_master() {

python3 - \
  "$P9C" \
  "$OUT" <<'PY'

from pathlib import Path
import csv
import math
import sys

P9C = Path(sys.argv[1])
OUT = Path(sys.argv[2])

ABF = P9C / "COLOC_ABF_MASTER_ANNOTATED_V2.tsv"
SUSIE = P9C / "COLOC_SUSIE_SIGNAL_PAIRS.tsv"

targets = [
    ("FGF5",     "BBJ_IS_L001", "CORE", "REGULATORY_EXPRESSION_PROTEIN"),
    ("ALDH2",    "BBJ_IS_L003", "CORE", "CODING_PROTEIN_METABOLIC"),
    ("SH3PXD2A", "BBJ_IS_L002", "CORE", "CELL_SPECIFIC_REGULATORY"),
    ("COL4A2",   "BBJ_IS_L004", "CORE", "VASCULAR_STRUCTURAL_REGULATORY"),
    ("COL4A1",   "BBJ_IS_L004", "CORE", "VASCULAR_STRUCTURAL_REGULATORY"),
    ("C4orf22",  "BBJ_IS_L001", "SECONDARY_SCREEN", "MOLECULAR_COLOC_SCREEN"),
    ("CALHM2",   "BBJ_IS_L002", "SECONDARY_SCREEN", "MOLECULAR_COLOC_SCREEN"),
    ("NEURL1",   "BBJ_IS_L002", "SECONDARY_SCREEN", "MOLECULAR_COLOC_SCREEN"),
    ("INA",      "BBJ_IS_L002", "SECONDARY_SCREEN", "MOLECULAR_COLOC_SCREEN"),
]

with ABF.open() as f:
    abf = list(csv.DictReader(f, delimiter="\t"))

with SUSIE.open() as f:
    susie = list(csv.DictReader(f, delimiter="\t"))


# dataset_key -> exact GTEx tissue label
dataset_to_tissue = {}

for r in abf:
    key = r.get("dataset_key", "")
    lab = r.get("tissue_label", "")

    if key and lab:
        dataset_to_tissue[key] = lab


def number(x):

    try:
        v = float(x)
    except Exception:
        return None

    return v if math.isfinite(v) else None


def best_abf(gene, metric):

    rows = [
        r for r in abf
        if r.get("gene_symbol") == gene
        and r.get("status") == "PASS"
        and number(r.get(metric)) is not None
    ]

    if not rows:
        return {}

    return max(
        rows,
        key=lambda r: float(r[metric])
    )


def best_susie(gene):

    rows = [
        r for r in susie
        if r.get("gene_symbol") == gene
        and number(r.get("PP.H4.abf")) is not None
    ]

    if not rows:
        return {}

    return max(
        rows,
        key=lambda r: float(r["PP.H4.abf"])
    )


out = []

for gene, locus, role, branch in targets:

    h4 = best_abf(
        gene,
        "PP.H4"
    )

    h3 = best_abf(
        gene,
        "PP.H3"
    )

    ss = best_susie(gene)

    susie_dataset = ss.get(
        "dataset_key",
        ""
    )

    row = {
        "gene": gene,
        "locus": locus,
        "role": role,
        "mechanism_branch": branch,

        "best_abf_h4":
            h4.get("PP.H4", ""),

        "h3_at_best_h4":
            h4.get("PP.H3", ""),

        "best_abf_dataset":
            h4.get("dataset_key", ""),

        "best_abf_tissue":
            h4.get("tissue_label", ""),

        "best_abf_tissue_class":
            h4.get("tissue_class", ""),

        "best_abf_qtl_n":
            h4.get("qtl_n", ""),

        "max_abf_h3":
            h3.get("PP.H3", ""),

        "h4_at_max_h3":
            h3.get("PP.H4", ""),

        "max_h3_dataset":
            h3.get("dataset_key", ""),

        "max_h3_tissue":
            h3.get("tissue_label", ""),

        "best_susie_h4":
            ss.get("PP.H4.abf", ""),

        "best_susie_dataset":
            susie_dataset,

        "best_susie_tissue":
            dataset_to_tissue.get(
                susie_dataset,
                susie_dataset.replace(
                    "GTEx_V8__",
                    ""
                ).replace(
                    "_",
                    " "
                )
                if susie_dataset
                else ""
            ),

        "best_susie_n_model":
            ss.get("n_model", ""),

        "best_susie_shared_variants":
            ss.get("n_shared", ""),

        "interpretation": "",
    }

    if gene == "FGF5":

        row["interpretation"] = (
            "ABF_AND_MULTISIGNAL_COLOC_SUPPORT"
        )

    elif gene == "ALDH2":

        row["interpretation"] = (
            "RS671_DOMINANT_EAS_SIGNAL;"
            "STRONG_H3_IN_VASCULAR_TISSUE;"
            "EXPRESSION_SHARED_CAUSAL_SIGNAL_NOT_ESTABLISHED"
        )

    elif gene == "SH3PXD2A":

        row["interpretation"] = (
            "MODERATE_BULK_EQTL_SUPPORT;"
            "CELL_SPECIFIC_REGULATION_PRIORITY"
        )

    elif gene in {"COL4A2", "COL4A1"}:

        row["interpretation"] = (
            "LIMITED_BULK_EQTL_COLOC;"
            "VASCULAR_CELL_REGULATORY_PRIORITY"
        )

    else:

        row["interpretation"] = (
            "SECONDARY_SCREEN_REQUIRES_ORTHOGONAL_VALIDATION"
        )

    out.append(row)


p = OUT / "IS_FUNCTIONAL_CONVERGENCE_MASTER_R1.tsv"

with p.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(out[0].keys())
    )

    w.writeheader()
    w.writerows(out)


print("===== REPAIRED MASTER =====")

for x in out:

    print(
        x["gene"],
        "ABF_H4=",
        x["best_abf_h4"],
        "tissue=",
        x["best_abf_tissue"],
        "MAX_H3=",
        x["max_abf_h3"],
        "H3_tissue=",
        x["max_h3_tissue"],
        "SuSiE_H4=",
        x["best_susie_h4"],
        "SuSiE_tissue=",
        x["best_susie_tissue"],
    )

PY

}


# =============================================================================
# 03. LITERATURE SOURCE MANIFEST
# =============================================================================

build_literature_manifest() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys

OUT = Path(sys.argv[1])

rows = [

    {
        "source_id": "GIGASTROKE_2022",
        "first_author": "Mishra",
        "year": "2022",
        "journal": "Nature",
        "title":
            "Stroke genetics informs drug discovery and risk prediction across ancestries",
        "doi":
            "10.1038/s41586-022-05165-3",
        "method":
            "Cross-ancestry GWAS; ancestry-specific SuSiE fine-mapping; MENTR; TWAS/PWAS; pQTL-MR",
        "ld_reference":
            "UK Biobank ~420000 EUR; BioBank Japan ~170000 EAS for fine-mapping",
        "direct_relevance":
            "FGF5; ALDH2/rs671; SH3PXD2A",
    },

    {
        "source_id": "YAO_2025",
        "first_author": "Yao",
        "year": "2025",
        "journal": "Stroke",
        "title":
            "Proteome-Wide Genetic Study in East Asians and Europeans Identified Multiple Therapeutic Targets for Ischemic Stroke",
        "doi":
            "10.1161/STROKEAHA.125.050982",
        "method":
            "Olink cis-pQTL; ancestry-specific two-sample MR; Bayesian colocalization; downstream PheWAS/druggability",
        "ld_reference":
            "Ancestry-specific molecular genetics",
        "direct_relevance":
            "FGF5; ALDH2",
    },

    {
        "source_id": "SEO_2024",
        "first_author": "Seo",
        "year": "2024",
        "journal": "Computational Biology and Chemistry",
        "title":
            "Bayesian colocalization of GWAS and eQTL signals reveals cell type-specific genes and regulatory variants for susceptibility to subtypes of ischemic stroke",
        "doi":
            "10.1016/j.compbiolchem.2024.108086",
        "method":
            "MEGASTROKE/GIGASTROKE GWAS integrated with single-cell brain/blood eQTL using Bayesian colocalization",
        "ld_reference":
            "",
        "direct_relevance":
            "Method benchmark for cell-type-specific regulatory follow-up",
    },

    {
        "source_id": "RANNIKMAE_2017",
        "first_author": "Rannikmae",
        "year": "2017",
        "journal": "Neurology",
        "title":
            "COL4A2 is associated with lacunar ischemic stroke and deep ICH",
        "doi":
            "10.1212/WNL.0000000000004560",
        "method":
            "European ancestry meta-analysis of common variants in cerebral small-vessel-disease genes",
        "ld_reference":
            "",
        "direct_relevance":
            "COL4A2",
    },

    {
        "source_id": "REID_2025",
        "first_author": "Reid",
        "year": "2025",
        "journal": "Neuron",
        "title":
            "Human brain vascular multi-omics elucidates disease-risk associations",
        "doi":
            "10.1016/j.neuron.2025.07.001",
        "method":
            "MultiVINE-seq: simultaneous RNA and chromatin accessibility profiling in vascular/perivascular/immune cells from human brain",
        "ld_reference":
            "",
        "direct_relevance":
            "COL4A2 vascular regulatory mechanism",
    },
]


p = OUT / "LITERATURE_SOURCE_MANIFEST.tsv"

with p.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(rows[0].keys())
    )

    w.writeheader()
    w.writerows(rows)


for x in rows:
    print(
        x["source_id"],
        x["doi"]
    )

PY

}


# =============================================================================
# 04. LITERATURE BENCHMARK RESULTS
# =============================================================================

build_benchmark() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys

OUT = Path(sys.argv[1])

rows = [

    {
        "gene": "FGF5",
        "source_id": "GIGASTROKE_2022",
        "published_method":
            "EAS/cross-ancestry GWAS with ancestry-matched fine-mapping; secondary MTAG",
        "published_result":
            "FGF5 locus associated with any stroke in East Asian secondary analysis",
        "our_result":
            "BBJ L001; GTEx cerebellar-hemisphere ABF H4 ~0.779; multi-signal SuSiE H4 ~0.752",
        "relationship":
            "CONSISTENT_LOCUS_PLUS_NEW_EXPRESSION_LAYER",
        "independence_note":
            "Our BBJ discovery and published EAS analyses may share underlying Japanese resources; do not label as fully independent replication",
    },

    {
        "gene": "FGF5",
        "source_id": "YAO_2025",
        "published_method":
            "EAS cis-pQTL MR plus colocalization",
        "published_result":
            "EAS ischemic stroke OR 1.14 per SD genetically predicted FGF5; P=2.1e-7; PP4=0.992",
        "our_result":
            "FGF5 brain eQTL ABF and SuSiE molecular-colocalization support",
        "relationship":
            "ORTHOGONAL_MOLECULAR_CONVERGENCE",
        "independence_note":
            "EAS ischemic-stroke outcome uses BBJ; molecular layer is orthogonal but outcome is not an independent replication of our BBJ signal",
    },

    {
        "gene": "ALDH2",
        "source_id": "GIGASTROKE_2022",
        "published_method":
            "Ancestry-specific SuSiE fine-mapping",
        "published_result":
            "rs671 (ALDH2) included among nonsynonymous variants in fine-mapped credible sets",
        "our_result":
            "rs671 is BBJ L003 lead; GRCh37 validated; JPT sensitivity shows no robust independent secondary GWS signal",
        "relationship":
            "DIRECT_VARIANT_LEVEL_CONVERGENCE",
        "independence_note":
            "Our fine-mapping remains limited by external EAS/JPT LD versus GIGASTROKE BBJ-scale ancestry-matched LD",
    },

    {
        "gene": "ALDH2",
        "source_id": "YAO_2025",
        "published_method":
            "EAS cis-pQTL MR plus colocalization",
        "published_result":
            "EAS ischemic stroke OR 1.60 per SD genetically predicted ALDH2; P=8.2e-11; PP4=0.991",
        "our_result":
            "Bulk eQTL shared-causal support is weak; max H3 is near 1 in vascular tissue, while rs671 is missense p.Glu504Lys",
        "relationship":
            "SUPPORTS_CODING_PROTEIN_BRANCH_NOT_SIMPLE_BULK_EQTL_MODEL",
        "independence_note":
            "Do not interpret pQTL result as proof that rs671 effect is mediated solely by plasma ALDH2 concentration",
    },

    {
        "gene": "SH3PXD2A",
        "source_id": "GIGASTROKE_2022",
        "published_method":
            "EUR/EAS SuSiE plus MENTR regulatory prediction",
        "published_result":
            "19 variants overlap between EUR/EAS credible sets; promoter variant predicted to modulate SH3PXD2A expression in macrophages",
        "our_result":
            "Bulk GTEx best H4 ~0.351 in tibial artery; no successful multi-signal molecular-coloc confirmation",
        "relationship":
            "CELL_SPECIFIC_FOLLOWUP_STRONGLY_JUSTIFIED",
        "independence_note":
            "Bulk tissue non-confirmation does not refute macrophage-specific regulation",
    },

    {
        "gene": "SH3PXD2A",
        "source_id": "SEO_2024",
        "published_method":
            "Single-cell brain/blood eQTL Bayesian colocalization",
        "published_result":
            "Stroke eQTL-colocalization was strongly cell-type dependent; different subtype signals localized to specific neural cell types",
        "our_result":
            "Bulk-tissue H4 is only moderate",
        "relationship":
            "METHOD_SUPPORT_FOR_CELL_TYPE_RESOLVED_QTL",
        "independence_note":
            "Seo et al. did not directly replicate SH3PXD2A as one of their five reported eGenes",
    },

    {
        "gene": "COL4A2",
        "source_id": "RANNIKMAE_2017",
        "published_method":
            "European case-control meta-analysis",
        "published_result":
            "rs9515201 associated with lacunar ischemic stroke; OR 1.17, P=6.62e-8",
        "our_result":
            "L004 positional/fine-mapping candidate; bulk GTEx best H4 ~0.331 in putamen",
        "relationship":
            "DIRECT_LOCUS_BIOLOGICAL_SUPPORT",
        "independence_note":
            "Published association is European and lacunar-stroke focused; population/phenotype differ from BBJ overall ischemic stroke",
    },

    {
        "gene": "COL4A2",
        "source_id": "REID_2025",
        "published_method":
            "Human brain vascular RNA plus chromatin multi-omics",
        "published_result":
            "Stroke SNP rs9515201 lies in cRE active in smooth muscle cells, pericytes and perivascular fibroblasts; regulatory links support COL4A2",
        "our_result":
            "Bulk-tissue molecular colocalization is limited",
        "relationship":
            "STRONG_SUPPORT_FOR_VASCULAR_CELL_REGULATORY_BRANCH",
        "independence_note":
            "Functional mapping is not equivalent to causal mediation; requires locus-specific validation",
    },

    {
        "gene": "COL4A1",
        "source_id": "REID_2025",
        "published_method":
            "Human brain vascular multi-omics",
        "published_result":
            "COL4A1/COL4A2 basement-membrane biology provides locus context, with strongest mapped common-risk evidence centered on COL4A2",
        "our_result":
            "COL4A1 bulk GTEx H4 weaker than COL4A2",
        "relationship":
            "SECONDARY_GENE_WITHIN_L004",
        "independence_note":
            "Do not merge COL4A1 and COL4A2 into a single gene-level causal claim",
    },
]


p = OUT / "IS_LITERATURE_BENCHMARK.tsv"

with p.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(rows[0].keys())
    )

    w.writeheader()
    w.writerows(rows)


for x in rows:

    print(
        x["gene"],
        x["source_id"],
        x["relationship"]
    )

PY

}


# =============================================================================
# 05. BUILD GENE COMPARISON MASTER
# =============================================================================

build_gene_master() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys

OUT = Path(sys.argv[1])

with (
    OUT /
    "IS_FUNCTIONAL_CONVERGENCE_MASTER_R1.tsv"
).open() as f:

    functional = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )

with (
    OUT /
    "IS_LITERATURE_BENCHMARK.tsv"
).open() as f:

    literature = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


lit_by_gene = {}

for x in literature:

    lit_by_gene.setdefault(
        x["gene"],
        []
    ).append(x)


rows = []

for x in functional:

    gene = x["gene"]

    lit = lit_by_gene.get(
        gene,
        []
    )

    rows.append({

        **x,

        "literature_source_count":
            len(lit),

        "literature_sources":
            ";".join(
                y["source_id"]
                for y in lit
            ),

        "literature_relationships":
            ";".join(
                y["relationship"]
                for y in lit
            ),

        "next_priority":
            (
                "PQTL_MR_BP_MEDIATION"
                if gene == "FGF5"

                else
                "CODING_PROTEIN_METABOLIC_FUNCTION"
                if gene == "ALDH2"

                else
                "MACROPHAGE_MONOCYTE_VASCULAR_CELL_QTL_ATAC"
                if gene == "SH3PXD2A"

                else
                "VASCULAR_SMC_PERICYTE_FIBROBLAST_ATAC"
                if gene in {
                    "COL4A2",
                    "COL4A1",
                }

                else
                "ORTHOGONAL_VALIDATION"
            ),
    })


p = OUT / "GENE_EVIDENCE_LITERATURE_MASTER.tsv"

with p.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(rows[0].keys())
    )

    w.writeheader()
    w.writerows(rows)


print(
    "GENES=",
    len(rows)
)

for x in rows:

    print(
        x["gene"],
        "sources=",
        x["literature_source_count"],
        "next=",
        x["next_priority"]
    )

PY

}


# =============================================================================
# 06. REPORT
# =============================================================================

build_report() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys

OUT = Path(sys.argv[1])


def read(name):

    with (
        OUT / name
    ).open() as f:

        return list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )


genes = read(
    "GENE_EVIDENCE_LITERATURE_MASTER.tsv"
)

bench = read(
    "IS_LITERATURE_BENCHMARK.tsv"
)


L = []

L.append(
    "# IS Phase9D Literature Benchmark\n\n"
)

L.append(
    "## Functional convergence master\n\n"
)

L.append(
    "|Gene|ABF H4|Best tissue|Max H3|Max-H3 tissue|SuSiE H4|Next|\n"
)

L.append(
    "|---|---:|---|---:|---|---:|---|\n"
)


for x in genes:

    L.append(
        f"|{x['gene']}|"
        f"{x['best_abf_h4']}|"
        f"{x['best_abf_tissue']}|"
        f"{x['max_abf_h3']}|"
        f"{x['max_h3_tissue']}|"
        f"{x['best_susie_h4']}|"
        f"{x['next_priority']}|\n"
    )


L.append(
    "\n## Literature benchmark\n\n"
)

for x in bench:

    L.append(
        f"### {x['gene']} — {x['source_id']}\n\n"
    )

    L.append(
        f"- Published method: {x['published_method']}\n"
    )

    L.append(
        f"- Published result: {x['published_result']}\n"
    )

    L.append(
        f"- Current study: {x['our_result']}\n"
    )

    L.append(
        f"- Relationship: **{x['relationship']}**\n"
    )

    L.append(
        f"- Caveat: {x['independence_note']}\n\n"
    )


L.append(
    "## Working thesis structure\n\n"
)

L.append(
    "1. FGF5: convergent GWAS-expression-protein pathway.\n"
)

L.append(
    "2. ALDH2/rs671: ancestry-specific coding/protein/metabolic mechanism.\n"
)

L.append(
    "3. SH3PXD2A: cell-specific immune/vascular regulatory mechanism.\n"
)

L.append(
    "4. COL4A2: cerebrovascular mural/stromal regulatory mechanism.\n"
)

L.append(
    "\nThe next analytical layer is cell-type-resolved regulatory validation rather than additional bulk-GTEx screening.\n"
)


(
    OUT /
    "PHASE9D_LITERATURE_BENCHMARK.md"
).write_text(
    "".join(L)
)


print(
    "".join(L)
)

PY

}


# =============================================================================
# 07. READINESS
# =============================================================================

build_readiness() {

cat > "$OUT/PHASE9D_READINESS.tsv" <<'EOF'
component	status
PHASE9D_LITERATURE_BENCHMARK	COMPLETE
PHASE9C_TISSUE_LABEL_REPAIR	COMPLETE
FGF5_INTERPRETATION	CONVERGENT_REGULATORY_EXPRESSION_PROTEIN
ALDH2_INTERPRETATION	RS671_DOMINANT_CODING_PROTEIN_METABOLIC
SH3PXD2A_INTERPRETATION	CELL_SPECIFIC_REGULATORY_PRIORITY
COL4A2_INTERPRETATION	VASCULAR_CELL_REGULATORY_PRIORITY
COL4A1_INTERPRETATION	SECONDARY_L004_GENE
L003_STATISTICAL_BRANCH	FROZEN
BULK_EQTL_STAGE	SUFFICIENT_FOR_PRIORITIZATION
NEXT_STAGE	PHASE10_FUNCTIONAL_CELLTYPE_ATAC
EOF

cat "$OUT/PHASE9D_READINESS.tsv"

}


# =============================================================================
# RUN
# =============================================================================

run_step \
  "01_preflight" \
  preflight

run_step \
  "02_repair_convergence_master" \
  repair_convergence_master

run_step \
  "03_literature_manifest" \
  build_literature_manifest

run_step \
  "04_literature_benchmark" \
  build_benchmark

run_step \
  "05_gene_master" \
  build_gene_master

run_step \
  "06_report" \
  build_report

run_step \
  "07_readiness" \
  build_readiness


echo
echo "================================================================================"
echo "PHASE9D COMPLETE"
echo "================================================================================"

echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"

echo
echo "===== REPAIRED MASTER ====="

column -t -s $'\t' \
  "$OUT/IS_FUNCTIONAL_CONVERGENCE_MASTER_R1.tsv" \
  2>/dev/null \
  || cat "$OUT/IS_FUNCTIONAL_CONVERGENCE_MASTER_R1.tsv"

echo
echo "===== READINESS ====="

column -t -s $'\t' \
  "$OUT/PHASE9D_READINESS.tsv" \
  2>/dev/null \
  || cat "$OUT/PHASE9D_READINESS.tsv"

echo
echo "===== OUTPUT FILES ====="

find "$OUT" \
  -maxdepth 1 \
  -type f \
  -printf '%f\t%s bytes\n' \
  | sort

echo
echo "OUTPUT_ROOT=$OUT"
echo "LOG_ROOT=$LOGDIR"
