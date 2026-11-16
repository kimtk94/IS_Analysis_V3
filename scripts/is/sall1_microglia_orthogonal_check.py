#!/usr/bin/env python3
"""Orthogonal SALL1/microglia marker crosscheck inside provisional CellTypist cells.

Based on Zhang et al., Mol Neurodegener 2026 DOI 10.1186/s13024-026-01003-6.
Detectable transcripts only: does NOT implement the paper's SALL1-high
classifier, establish microglia specificity in humans, or make cell-type calls.
All outputs are public GEO sample-level aggregates without barcodes.
"""
from collections import Counter

MICROGLIA_GENES = ("SALL1", "TMEM119", "P2RY12")
PERIPHERAL_MY_GENES = ("FCN1", "CD14", "S100A8")
MARKERS = MICROGLIA_GENES + PERIPHERAL_MY_GENES
MODEL_LABELS = ("Monocytes", "Macrophages", "Mono-mac", "T cells", "DC")
MYELOID_LABELS = ("Monocytes", "Macrophages", "Mono-mac")


def check_binary_markers(markers, labels):
    """marker arrays -> coexpression metrics, with no inferred cell identities."""
    labels = list(labels)
    n = len(labels)
    if n == 0:
        raise ValueError("No QC cells")
    if set(markers) != set(MARKERS):
        raise ValueError("Expected all six SALL1/peripheral marker genes")
    for gene in MARKERS:
        if len(markers[gene]) != n:
            raise ValueError("Marker array length differs from model labels: "+gene)
        if any(v not in (True,False,0,1) for v in markers[gene]):
            raise ValueError("Marker evidence is not binary: "+gene)
    result = []
    for label in (*MODEL_LABELS, "ALL_PREDICTED_MY_MARKERS"):
        idx = [
            i for i, pred in enumerate(labels)
            if pred == label or (label=="ALL_PREDICTED_MY_MARKERS" and pred in MYELOID_LABELS)
        ]
        count=len(idx)
        if count==0:
            result.append({"model_provisional_label":label,"n_barcodes":0,
                           "status":"NO_PREDICTED_CELLS","rates":None,
                           "counts":None})
            continue
        present=lambda name,i:bool(markers[name][i])
        n_sall1=sum(present("SALL1",i) for i in idx)
        n_sall1_and_canon=sum(present("SALL1",i) and
                              (present("TMEM119",i) or present("P2RY12",i))
                              for i in idx)
        n_strict_triple=sum(all(present(g,i) for g in MICROGLIA_GENES) for i in idx)
        n_canonical_pair=sum(present("TMEM119",i) and present("P2RY12",i) for i in idx)
        n_periph_two=sum(sum(present(g,i) for g in PERIPHERAL_MY_GENES)>=2 for i in idx)
        n_ambiguous=sum(present("SALL1",i) and
                        (present("TMEM119",i) or present("P2RY12",i)) and
                        sum(present(g,i) for g in PERIPHERAL_MY_GENES)>=2 for i in idx)
        counts={
            "SALL1_any":n_sall1,
            "SALL1_plus_one_canonical":n_sall1_and_canon,
            "SALL1_TMEM119_P2RY12_triple":n_strict_triple,
            "TMEM119_P2RY12_pair_regardless_SALL1":n_canonical_pair,
            "at_least_two_FCN1_CD14_S100A8":n_periph_two,
            "both_SALL1_canonical_and_peripheral":n_ambiguous,
        }
        for val in counts.values():
            if not 0<=val<=count:
                raise ValueError("Invalid co-detection count")
        result.append({
            "model_provisional_label":label,"n_barcodes":count,
            "status":"DESCRIPTIVE_MARKER_EVIDENCE_ONLY",
            "counts":counts,
            "rates":{key:round(v/count,6) for key,v in counts.items()},
        })
    # Mark clearly that SALL1-positive cells cannot establish resident lineage.
    return {
        "scope":"PUBLIC_GEO_PER_GSM_ORTHOGONAL_BINARY_RNA_DETECTION",
        "source_doi":"10.1186/s13024-026-01003-6",
        "paper_SALL1_high_classifier_reproduced":False,
        "paper_empirical_sensitivity_specificity_applies_to_this_dataset":False,
        "RNA_dropout_possible":True,
        "celltypes_validated":False,
        "myeloid_resident_versus_infiltrating_resolved":False,
        "n_total_barcodes":n,
        "model_categories":result,
    }


def ann_data_binary_markers(adata, labels):
    """Return sample-level grouped marker evidence; source row ordering must match predictions."""
    import numpy as np
    features=list(adata.var_names)
    if len(set(features))!=len(features):
        raise ValueError("Duplicate 10x feature names")
    names={}
    for gene in MARKERS:
        ix=[i for i,name in enumerate(features) if name==gene]
        if len(ix)!=1:
            raise ValueError("Required marker feature missing or ambiguous: "+gene)
        names[gene]=ix[0]
    mark={}
    for gene,ix in names.items():
        m=adata.X[:,ix]
        dense=m.toarray().ravel() if hasattr(m,"toarray") else np.asarray(m).ravel()
        if not np.all(np.isfinite(dense)) or (dense<0).any():
            raise ValueError("Invalid normalized marker expression: "+gene)
        mark[gene]=(dense>0).tolist()
    return check_binary_markers(mark,labels)


if __name__=="__main__":
    raise SystemExit("Import in CellTypist pilot; this module does not download protected data.")
