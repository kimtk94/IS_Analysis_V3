#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

P9D="$ROOT/results/is/stage5_functional/phase9d_literature_benchmark"

OUT="$ROOT/results/is/stage5_functional/phase9e_celltype_resource_audit"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGDIR="$ROOT/logs/is_phase9e_celltype/$RUN_ID"
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

  echo "===== PHASE9D STATUS ====="

  if [[ ! -s "$P9D/STEP_STATUS.tsv" ]]; then
    echo "PHASE9D_STATUS_MISSING"
    return 1
  fi

  cat "$P9D/STEP_STATUS.tsv"

  FAIL_N="$(
    awk -F'\t' '
      NR>1 && $2!="PASS" {n++}
      END {print n+0}
    ' "$P9D/STEP_STATUS.tsv"
  )"

  echo
  echo "PHASE9D_FAIL_N=$FAIL_N"

  if [[ "$FAIL_N" -ne 0 ]]; then
    echo "PHASE9D_NOT_CLEAN"
    return 1
  fi


  echo
  echo "===== GENE MASTER ====="

  MASTER="$P9D/GENE_EVIDENCE_LITERATURE_MASTER.tsv"

  if [[ -s "$MASTER" ]]; then
    echo "FOUND $MASTER"
  else
    echo "GENE_MASTER_MISSING"
    return 1
  fi


  echo
  echo "===== PYTHON ====="

  python3 --version

  python3 - <<'PY'
mods = [
    "numpy",
    "pandas",
    "scipy",
    "anndata",
    "scanpy",
]

for x in mods:

    try:
        m = __import__(x)
        print(
            x,
            "FOUND",
            getattr(
                m,
                "__version__",
                ""
            )
        )
    except Exception as e:
        print(
            x,
            "MISSING",
            type(e).__name__
        )
PY

  echo
  echo "PREFLIGHT=PASS"
}


# =============================================================================
# 02. LOCAL SINGLE-CELL / ATAC / SPATIAL INVENTORY
# =============================================================================

inventory_local_resources() {

python3 - "$ROOT" "$OUT" <<'PY'

from pathlib import Path
import csv
import os
import sys


ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])

search_roots = [

    ROOT / "data/is",

    ROOT / "data",

    ROOT / "results/is",

    ROOT / "resources",
]


extensions = {

    ".h5ad",
    ".h5",
    ".loom",
    ".rds",
    ".rda",
    ".rdata",
    ".mtx",
    ".h5seurat",
    ".bed",
    ".bed.gz",
    ".fragments.tsv.gz",
}


accessions = [

    "GSE174574",
    "GSE189432",
    "GSE225948",
    "GSE234052",
    "GSE247474",
    "GSE250245",
    "GSE266033",
]


keywords = [

    "stroke",
    "ischemic",
    "ischaemic",
    "mca",
    "mcao",
    "brain",
    "vascular",
    "endothelial",
    "pericyte",
    "smooth_muscle",
    "smoothmuscle",
    "fibroblast",
    "microglia",
    "macrophage",
    "singlecell",
    "single_cell",
    "scrna",
    "snrna",
    "snatac",
    "atac",
    "spatial",
]


rows = []
seen = set()


def classify(path):

    low = str(path).lower()

    if "atac" in low:
        modality = "ATAC"

    elif (
        "spatial" in low
        or "visium" in low
    ):
        modality = "SPATIAL"

    elif (
        "scrna" in low
        or "snrna" in low
        or path.suffix.lower()
        in {
            ".h5ad",
            ".loom",
        }
    ):
        modality = "SINGLE_CELL_RNA"

    else:
        modality = "UNKNOWN"


    if any(
        x in low
        for x in [
            "stroke",
            "ischemic",
            "ischaemic",
            "mcao",
        ]
    ):

        disease_context = "STROKE_RELATED"

    else:
        disease_context = "REFERENCE_OR_UNKNOWN"


    if any(
        x in low
        for x in [
            "human",
            "homo",
            "hca",
        ]
    ):

        species_hint = "HUMAN"

    elif any(
        x in low
        for x in [
            "mouse",
            "murine",
            "mcao",
        ]
    ):

        species_hint = "MOUSE"

    elif "rat" in low:

        species_hint = "RAT"

    else:
        species_hint = "UNKNOWN"


    return (
        modality,
        disease_context,
        species_hint,
    )


for root in search_roots:

    if not root.exists():
        continue

    for current, dirs, files in os.walk(root):

        cur = Path(current)

        # Skip environments/caches.
        low_cur = str(cur).lower()

        if any(
            x in low_cur
            for x in [
                "/.git",
                "/venv",
                "/env/",
                "/site-packages",
                "/node_modules",
            ]
        ):
            dirs[:] = []
            continue

        for fn in files:

            p = cur / fn
            low = str(p).lower()

            suffix_match = any(
                low.endswith(x)
                for x in extensions
            )

            accession_match = any(
                x.lower() in low
                for x in accessions
            )

            keyword_match = any(
                x in low
                for x in keywords
            )

            if not (
                suffix_match
                and (
                    accession_match
                    or keyword_match
                )
            ):
                continue


            key = str(p)

            if key in seen:
                continue

            seen.add(key)

            modality, disease, species = classify(p)

            acc = ";".join(
                x
                for x in accessions
                if x.lower() in low
            )

            try:
                size = p.stat().st_size
            except Exception:
                size = ""

            rows.append({

                "path":
                    str(p),

                "bytes":
                    size,

                "accession_hint":
                    acc,

                "modality_hint":
                    modality,

                "disease_context_hint":
                    disease,

                "species_hint":
                    species,
            })


rows.sort(
    key=lambda x:
        (
            x["modality_hint"],
            x["path"],
        )
)


outfile = (
    OUT /
    "LOCAL_CELLTYPE_RESOURCE_INVENTORY.tsv"
)


with outfile.open(
    "w",
    newline=""
) as f:

    fields = [
        "path",
        "bytes",
        "accession_hint",
        "modality_hint",
        "disease_context_hint",
        "species_hint",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(rows)


print(
    "LOCAL_RESOURCE_FILES=",
    len(rows)
)


for x in rows[:200]:

    print(
        x["modality_hint"],
        x["disease_context_hint"],
        x["species_hint"],
        x["accession_hint"],
        x["bytes"],
        x["path"],
    )

PY

}


# =============================================================================
# 03. ACCESSION INVENTORY
# =============================================================================

audit_accessions() {

python3 - "$OUT" "$ROOT" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])
ROOT = Path(sys.argv[2])


datasets = [

    {
        "accession": "GSE174574",
        "species_class": "ANIMAL",
        "modality_expected": "sc/snRNA",
        "analysis_role":
            "STROKE_DISEASE_CELLSTATE_VALIDATION",
        "primary_use":
            "ischemic injury cell-state evidence",
    },

    {
        "accession": "GSE189432",
        "species_class": "ANIMAL_OR_STROKE_REFERENCE",
        "modality_expected": "sc/snRNA",
        "analysis_role":
            "STROKE_DISEASE_VALIDATION",
        "primary_use":
            "candidate-gene disease-state validation",
    },

    {
        "accession": "GSE225948",
        "species_class": "ANIMAL",
        "modality_expected": "sc/snRNA",
        "analysis_role":
            "STROKE_DISEASE_VALIDATION",
        "primary_use":
            "stroke cell-state validation",
    },

    {
        "accession": "GSE234052",
        "species_class": "ANIMAL",
        "modality_expected": "scRNA",
        "analysis_role":
            "PERICYTE_STROKE_VALIDATION",
        "primary_use":
            "COL4A2/COL4A1 mural-cell disease response",
    },

    {
        "accession": "GSE247474",
        "species_class": "ANIMAL",
        "modality_expected": "scRNA",
        "analysis_role":
            "ASTROCYTE_STROKE_VALIDATION",
        "primary_use":
            "brain injury response context",
    },

    {
        "accession": "GSE250245",
        "species_class": "RAT",
        "modality_expected": "snRNA",
        "analysis_role":
            "STROKE_DISEASE_VALIDATION",
        "primary_use":
            "cross-species disease-state validation",
    },

    {
        "accession": "GSE266033",
        "species_class": "ANIMAL",
        "modality_expected": "sc/snRNA",
        "analysis_role":
            "RECOVERY_VASCULAR_VALIDATION",
        "primary_use":
            "stroke recovery/cerebrovascular context",
    },
]


inventory = (
    OUT /
    "LOCAL_CELLTYPE_RESOURCE_INVENTORY.tsv"
)

local_rows = []

if inventory.exists():

    with inventory.open() as f:

        local_rows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )


for d in datasets:

    accession = d["accession"]

    hits = [
        x
        for x in local_rows
        if accession
        in (
            x.get(
                "accession_hint",
                ""
            )
        )
        or accession.lower()
        in x.get(
            "path",
            ""
        ).lower()
    ]

    d["local_file_n"] = len(
        hits
    )

    d["local_status"] = (
        "LOCAL_CANDIDATE_FOUND"
        if hits
        else "NOT_FOUND_LOCALLY"
    )

    d["local_paths"] = ";".join(
        x["path"]
        for x in hits[:20]
    )


outfile = (
    OUT /
    "STROKE_DATASET_ACCESSION_AUDIT.tsv"
)


with outfile.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            datasets[0].keys()
        )
    )

    w.writeheader()
    w.writerows(
        datasets
    )


for x in datasets:

    print(
        x["accession"],
        x["analysis_role"],
        x["local_status"],
        "n=",
        x["local_file_n"],
    )

PY

}


# =============================================================================
# 04. TARGET → CELL-TYPE PLAN
# =============================================================================

build_target_plan() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])


rows = [

    {
        "gene": "FGF5",
        "locus": "BBJ_IS_L001",
        "genetic_model":
            "REGULATORY_EXPRESSION_PROTEIN",
        "priority_cell_types":
            "vascular smooth muscle;pericyte;endothelial;brain glia",
        "human_reference_goal":
            "localize expression and accessible regulatory context",
        "stroke_disease_goal":
            "test injury-state regulation without claiming causal mediation",
        "priority_modality":
            "human sc/snRNA;human vascular ATAC;pQTL",
    },

    {
        "gene": "ALDH2",
        "locus": "BBJ_IS_L003",
        "genetic_model":
            "RS671_CODING_PROTEIN_METABOLIC",
        "priority_cell_types":
            "vascular cells;neurons;astrocytes;microglia",
        "human_reference_goal":
            "localize ALDH2 expression only",
        "stroke_disease_goal":
            "contextualize disease-state expression;not primary causal test",
        "priority_modality":
            "coding annotation;protein/pQTL;scRNA supportive",
    },

    {
        "gene": "SH3PXD2A",
        "locus": "BBJ_IS_L002",
        "genetic_model":
            "CELL_SPECIFIC_REGULATORY",
        "priority_cell_types":
            "macrophage;monocyte;microglia;endothelial;vascular immune",
        "human_reference_goal":
            "identify relevant expressing/regulatory cell type",
        "stroke_disease_goal":
            "test disease-associated immune-cell regulation",
        "priority_modality":
            "human cell-type eQTL;sc/snRNA;ATAC",
    },

    {
        "gene": "COL4A2",
        "locus": "BBJ_IS_L004",
        "genetic_model":
            "VASCULAR_STRUCTURAL_REGULATORY",
        "priority_cell_types":
            "vascular smooth muscle;pericyte;perivascular fibroblast;endothelial",
        "human_reference_goal":
            "confirm mural/stromal localization and regulatory accessibility",
        "stroke_disease_goal":
            "test vascular injury/remodelling response",
        "priority_modality":
            "human vascular scRNA;snATAC;stroke pericyte scRNA",
    },

    {
        "gene": "COL4A1",
        "locus": "BBJ_IS_L004",
        "genetic_model":
            "SECONDARY_VASCULAR_STRUCTURAL",
        "priority_cell_types":
            "vascular smooth muscle;pericyte;perivascular fibroblast",
        "human_reference_goal":
            "compare localization with COL4A2",
        "stroke_disease_goal":
            "distinguish COL4A1 vs COL4A2 disease-response patterns",
        "priority_modality":
            "human vascular scRNA;snATAC",
    },

    {
        "gene": "CALHM2",
        "locus": "BBJ_IS_L002",
        "genetic_model":
            "SECONDARY_MOLECULAR_COLOC",
        "priority_cell_types":
            "brain cell types;vascular cells",
        "human_reference_goal":
            "resolve strong cerebellar H4 versus arterial H3 discordance",
        "stroke_disease_goal":
            "orthogonal validation only",
        "priority_modality":
            "human scRNA;cell-type eQTL",
    },

    {
        "gene": "NEURL1",
        "locus": "BBJ_IS_L002",
        "genetic_model":
            "SECONDARY_MOLECULAR_COLOC",
        "priority_cell_types":
            "neuronal populations;brain cell types",
        "human_reference_goal":
            "resolve cerebellar signal",
        "stroke_disease_goal":
            "orthogonal validation only",
        "priority_modality":
            "human scRNA;cell-type eQTL",
    },

    {
        "gene": "C4orf22",
        "locus": "BBJ_IS_L001",
        "genetic_model":
            "SECONDARY_MOLECULAR_COLOC",
        "priority_cell_types":
            "brain cell types",
        "human_reference_goal":
            "resolve cerebellar signal",
        "stroke_disease_goal":
            "orthogonal validation only",
        "priority_modality":
            "human scRNA;cell-type eQTL",
    },

    {
        "gene": "INA",
        "locus": "BBJ_IS_L002",
        "genetic_model":
            "SECONDARY_MOLECULAR_COLOC",
        "priority_cell_types":
            "neurons",
        "human_reference_goal":
            "test neuronal localization",
        "stroke_disease_goal":
            "orthogonal validation only",
        "priority_modality":
            "human neuronal scRNA/snRNA",
    },
]


with (
    OUT /
    "CELLTYPE_TARGET_PLAN.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=list(
            rows[0].keys()
        )
    )

    w.writeheader()
    w.writerows(rows)


for x in rows:

    print(
        x["gene"],
        "->",
        x["priority_cell_types"],
        "|",
        x["priority_modality"],
    )

PY

}


# =============================================================================
# 05. H5AD HEADER AUDIT — NO FULL MATRIX LOAD
# =============================================================================

audit_h5ad() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])

inventory = (
    OUT /
    "LOCAL_CELLTYPE_RESOURCE_INVENTORY.tsv"
)


try:
    import anndata as ad

    ANNDATA_OK = True

except Exception:

    ANNDATA_OK = False


rows = []


if inventory.exists():

    with inventory.open() as f:

        resources = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

else:

    resources = []


h5ads = [

    Path(x["path"])

    for x in resources

    if x["path"].lower().endswith(
        ".h5ad"
    )
]


print(
    "H5AD_FILES=",
    len(h5ads)
)

print(
    "ANNDATA_AVAILABLE=",
    ANNDATA_OK
)


for p in h5ads:

    row = {

        "path":
            str(p),

        "status":
            "",

        "n_obs":
            "",

        "n_vars":
            "",

        "obs_columns":
            "",

        "candidate_celltype_columns":
            "",

        "candidate_condition_columns":
            "",

        "target_genes_present":
            "",
    }


    if not ANNDATA_OK:

        row["status"] = (
            "ANNDATA_NOT_INSTALLED"
        )

        rows.append(row)
        continue


    try:

        a = ad.read_h5ad(
            p,
            backed="r"
        )

        row["n_obs"] = (
            a.n_obs
        )

        row["n_vars"] = (
            a.n_vars
        )

        obs = list(
            a.obs.columns
        )

        row["obs_columns"] = ";".join(
            obs[:100]
        )


        cell_candidates = [

            x for x in obs

            if any(
                k in x.lower()

                for k in [
                    "cell_type",
                    "celltype",
                    "cell type",
                    "subclass",
                    "supercluster",
                    "cluster",
                    "annotation",
                ]
            )
        ]


        condition_candidates = [

            x for x in obs

            if any(
                k in x.lower()

                for k in [
                    "condition",
                    "disease",
                    "stroke",
                    "case",
                    "status",
                    "time",
                    "group",
                    "treatment",
                ]
            )
        ]


        row[
            "candidate_celltype_columns"
        ] = ";".join(
            cell_candidates
        )


        row[
            "candidate_condition_columns"
        ] = ";".join(
            condition_candidates
        )


        targets = [
            "FGF5",
            "ALDH2",
            "SH3PXD2A",
            "COL4A2",
            "COL4A1",
            "CALHM2",
            "NEURL1",
            "C4orf22",
            "INA",
        ]


        var_set = set(
            map(
                str,
                a.var_names
            )
        )


        row[
            "target_genes_present"
        ] = ";".join(
            g
            for g in targets
            if g in var_set
        )


        row["status"] = (
            "HEADER_AUDIT_PASS"
        )


        if getattr(
            a,
            "file",
            None
        ) is not None:

            try:
                a.file.close()
            except Exception:
                pass


    except Exception as e:

        row["status"] = (
            "AUDIT_ERROR:"
            + type(e).__name__
        )


    rows.append(row)


outfile = (
    OUT /
    "H5AD_HEADER_AUDIT.tsv"
)


with outfile.open(
    "w",
    newline=""
) as f:

    fields = [
        "path",
        "status",
        "n_obs",
        "n_vars",
        "obs_columns",
        "candidate_celltype_columns",
        "candidate_condition_columns",
        "target_genes_present",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(rows)


for x in rows:

    print(
        x["status"],
        x["n_obs"],
        x["n_vars"],
        x["candidate_celltype_columns"],
        x["candidate_condition_columns"],
        x["target_genes_present"],
        x["path"],
    )

PY

}


# =============================================================================
# 06. RESOURCE → BRANCH MATCHING
# =============================================================================

match_resources() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])


with (
    OUT /
    "LOCAL_CELLTYPE_RESOURCE_INVENTORY.tsv"
).open() as f:

    inv = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )


rules = {

    "FGF5":
        [
            "brain",
            "vascular",
            "atac",
        ],

    "ALDH2":
        [
            "brain",
            "vascular",
        ],

    "SH3PXD2A":
        [
            "vascular",
            "microglia",
            "macrophage",
            "stroke",
            "atac",
        ],

    "COL4A2":
        [
            "vascular",
            "pericyte",
            "smooth",
            "fibroblast",
            "atac",
            "stroke",
        ],

    "COL4A1":
        [
            "vascular",
            "pericyte",
            "smooth",
            "fibroblast",
            "atac",
            "stroke",
        ],

}


rows = []


for gene, keys in rules.items():

    for x in inv:

        low = x["path"].lower()

        score = sum(
            k in low
            for k in keys
        )

        if score == 0:
            continue

        rows.append({

            "gene":
                gene,

            "resource_path":
                x["path"],

            "resource_modality":
                x["modality_hint"],

            "disease_context":
                x[
                    "disease_context_hint"
                ],

            "accession_hint":
                x[
                    "accession_hint"
                ],

            "match_score":
                score,

            "match_terms":
                ";".join(
                    k
                    for k in keys
                    if k in low
                ),
        })


rows.sort(
    key=lambda x:
        (
            x["gene"],
            -int(
                x["match_score"]
            ),
            x["resource_path"],
        )
)


with (
    OUT /
    "TARGET_RESOURCE_MATCH.tsv"
).open(
    "w",
    newline=""
) as f:

    fields = [
        "gene",
        "resource_path",
        "resource_modality",
        "disease_context",
        "accession_hint",
        "match_score",
        "match_terms",
    ]

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields,
    )

    w.writeheader()
    w.writerows(rows)


for x in rows[:100]:

    print(
        x["gene"],
        "score=",
        x["match_score"],
        x["resource_modality"],
        x["accession_hint"],
        x["resource_path"],
    )

PY

}


# =============================================================================
# 07. READINESS
# =============================================================================

build_readiness() {

python3 - "$OUT" <<'PY'

from pathlib import Path
import csv
import sys


OUT = Path(sys.argv[1])


def read(name):

    p = OUT / name

    if not p.exists():
        return []

    with p.open() as f:

        return list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )


inv = read(
    "LOCAL_CELLTYPE_RESOURCE_INVENTORY.tsv"
)

h5 = read(
    "H5AD_HEADER_AUDIT.tsv"
)

geo = read(
    "STROKE_DATASET_ACCESSION_AUDIT.tsv"
)


h5_ready = [

    x for x in h5

    if x.get(
        "status"
    ) == "HEADER_AUDIT_PASS"
]


stroke_local = [

    x for x in geo

    if x.get(
        "local_status"
    ) == "LOCAL_CANDIDATE_FOUND"
]


atac_local = [

    x for x in inv

    if x.get(
        "modality_hint"
    ) == "ATAC"
]


spatial_local = [

    x for x in inv

    if x.get(
        "modality_hint"
    ) == "SPATIAL"
]


rows = [

    (
        "PHASE9E_RESOURCE_AUDIT",
        "COMPLETE"
    ),

    (
        "LOCAL_SINGLECELL_FILES",
        str(
            len(inv)
        )
    ),

    (
        "H5AD_HEADER_AUDIT_PASS",
        str(
            len(h5_ready)
        )
    ),

    (
        "STROKE_ACCESSIONS_LOCAL",
        str(
            len(stroke_local)
        )
    ),

    (
        "LOCAL_ATAC_CANDIDATES",
        str(
            len(atac_local)
        )
    ),

    (
        "LOCAL_SPATIAL_CANDIDATES",
        str(
            len(spatial_local)
        )
    ),

    (
        "SH3PXD2A_NEXT",
        "HUMAN_IMMUNE_VASCULAR_CELL_LOCALIZATION"
    ),

    (
        "COL4A2_COL4A1_NEXT",
        "HUMAN_VASCULAR_MURAL_CELL_ATAC_LOCALIZATION"
    ),

    (
        "FGF5_NEXT",
        "HUMAN_BRAIN_VASCULAR_EXPRESSION_LOCALIZATION"
    ),

    (
        "ALDH2_NEXT",
        "CODING_PROTEIN_PRIMARY_SC_SUPPORTIVE_ONLY"
    ),

    (
        "DISEASE_SC_ANALYSIS",
        (
            "READY_FOR_TARGETED_ANALYSIS"
            if stroke_local
            else "DATA_ACQUISITION_REQUIRED"
        )
    ),

    (
        "HUMAN_REFERENCE_SC_ANALYSIS",
        (
            "READY_FOR_TARGETED_ANALYSIS"
            if h5_ready
            else "DATA_ACQUISITION_REQUIRED"
        )
    ),

    (
        "ATAC_ANALYSIS",
        (
            "READY_FOR_TARGETED_ANALYSIS"
            if atac_local
            else "DATA_ACQUISITION_REQUIRED"
        )
    ),

    (
        "NEXT_STAGE",
        "PHASE9F_TARGETED_CELLTYPE_ANALYSIS"
    ),
]


with (
    OUT /
    "PHASE9E_READINESS.tsv"
).open(
    "w",
    newline=""
) as f:

    w = csv.writer(
        f,
        delimiter="\t"
    )

    w.writerow([
        "component",
        "status"
    ])

    w.writerows(rows)


for a, b in rows:
    print(
        f"{a} = {b}"
    )

PY

}


# =============================================================================
# RUN
# =============================================================================

run_step \
  "01_preflight" \
  preflight

run_step \
  "02_local_resource_inventory" \
  inventory_local_resources

run_step \
  "03_accession_audit" \
  audit_accessions

run_step \
  "04_target_plan" \
  build_target_plan

run_step \
  "05_h5ad_header_audit" \
  audit_h5ad

run_step \
  "06_target_resource_match" \
  match_resources

run_step \
  "07_readiness" \
  build_readiness


# =============================================================================
# FINAL
# =============================================================================

echo
echo "================================================================================"
echo "PHASE9E CELL-TYPE RESOURCE AUDIT COMPLETE"
echo "================================================================================"

echo
echo "===== STEP STATUS ====="

column -t -s $'\t' \
  "$STATUS" \
  2>/dev/null \
  || cat "$STATUS"


echo
echo "===== ACCESSION AUDIT ====="

column -t -s $'\t' \
  "$OUT/STROKE_DATASET_ACCESSION_AUDIT.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/STROKE_DATASET_ACCESSION_AUDIT.tsv"


echo
echo "===== H5AD AUDIT ====="

column -t -s $'\t' \
  "$OUT/H5AD_HEADER_AUDIT.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/H5AD_HEADER_AUDIT.tsv"


echo
echo "===== TARGET PLAN ====="

column -t -s $'\t' \
  "$OUT/CELLTYPE_TARGET_PLAN.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/CELLTYPE_TARGET_PLAN.tsv"


echo
echo "===== READINESS ====="

column -t -s $'\t' \
  "$OUT/PHASE9E_READINESS.tsv" \
  2>/dev/null \
  || cat \
       "$OUT/PHASE9E_READINESS.tsv"


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

