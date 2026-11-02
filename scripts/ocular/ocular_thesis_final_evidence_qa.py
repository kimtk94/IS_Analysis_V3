#!/usr/bin/env python3
"""Fail-closed, source-linked QA of actual Ocular TM thesis figures and claims.

Use only already-produced public aggregate GitHub Action artifacts. This is an
integrity/provenance check, not a new inferential test, peer review, or validation
of biological causality. The original author DEG80 identifiers are not required.
"""
import argparse
import csv
import hashlib
import json
import math
import struct
import tempfile
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

RUNS = {
    "figures": {"run": 41208972011, "artifact": 13045363678},
    "strata": {"run": 41389639386, "artifact": 13116138830},
    "strain": {"run": 41423487666, "artifact": 13130831706},
    "design": {"run": 41430965326, "artifact": 13133928461},
}
PINNED_FIGURE_INPUTS = {
    "three_core": "926b7899639d9b759d03688b911b89becf6b776103927a226e4f301415a86180",
    "three_emt": "edac2a23c47e40fde419214c8c45fed6f1bd33861a72f8c37c789db2f7d88e1d",
    "cutoff_panel": "ffaf1f6376f8c1fb7e041b9890a3723e3396da3a1f365438556f07eea6794e1e",
    "dual_panel": "6b42e89ef03d8284e6875310a4a33fa193abe37fbfd50055eefb77f9d332777e",
}
CORE = {"MYOC": 1, "TGFBI": -1, "MMP2": -1}
FIGURES = (
    "Fig1_three_study_core_transcript_signs",
    "Fig2_BeamA_cutoff_and_matched55_annotation",
    "Fig3_EMT_194_signed_cross_experiment_concordance",
)


def find_one(directory, name):
    matches = list(directory.rglob(name))
    if len(matches) != 1:
        raise ValueError("EXPECTED_EXACTLY_ONE_FILE_" + name + "_" + str(len(matches)))
    return matches[0]


def load_json(directory, name):
    return json.loads(find_one(directory, name).read_text(encoding="utf8"))


def read_tsv(directory, name):
    with find_one(directory, name).open(newline="", encoding="utf8") as f:
        data = list(csv.DictReader(f, delimiter="\t"))
    if not data:
        raise ValueError("EMPTY_TSV_" + name)
    return data


def eq(actual, expected, key):
    if actual != expected:
        raise ValueError("CLAIM_DRIFT_" + key + "_" + str(actual) + "_EXPECTED_" + str(expected))


def close(actual, expected, key, tol=1e-5):
    if not (math.isfinite(float(actual)) and abs(float(actual)-expected)<tol):
        raise ValueError("NUMBER_DRIFT_" + key + "_" + str(actual))


def check_png(path):
    data = path.read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("INVALID_PNG_MAGIC_" + path.name)
    off, dimensions, iend = 8, None, False
    while off<len(data):
        if off+12>len(data):
            raise ValueError("TRUNCATED_PNG_" + path.name)
        n=struct.unpack(">I",data[off:off+4])[0]
        chunk=data[off+4:off+8]
        body=data[off+8:off+8+n]
        crc=struct.unpack(">I",data[off+8+n:off+12+n])[0] if off+12+n<=len(data) else -1
        if zlib.crc32(chunk+body)&0xffffffff != crc:
            raise ValueError("PNG_CHUNK_CRC_" + path.name)
        if chunk==b"IHDR":
            dimensions=struct.unpack(">II",body[:8])
        off+=n+12
        if chunk==b"IEND":
            iend=True
            break
    if not iend or off!=len(data) or dimensions is None or min(dimensions)<600:
        raise ValueError("TRUNCATED_OR_TOO_SMALL_PNG_" + path.name)
    return dimensions


def check_svg(path):
    raw=path.read_bytes()
    root=ET.fromstring(raw)
    if not root.tag.endswith("svg") or len(raw)<5000:
        raise ValueError("NOT_A_COMPLETE_FIGURE_SVG_" + path.name)
    return len(raw)


def figure_integrity(d):
    meta=load_json(d, "FIGURE_SOURCE_MANIFEST.json")
    eq(meta.get("status"), "PASS_FROZEN_ARTIFACT_FIGURES", "figure_source_status")
    eq(meta.get("all_three_emt_signed"),124,"figure_3_emt_old_124")
    eq(meta.get("core_three_3studies_signed"),["MYOC","TGFBI","MMP2"],"figure_1_core")
    if not meta.get("ACTA2_7day_direction_discordant") or not meta.get("same_55_donor_sensitivity"):
        raise ValueError("FIGURE_STUDY_CLAIM_DRIFT")
    hashes=meta.get("source_file_hashes",{})
    for name, sha in PINNED_FIGURE_INPUTS.items():
        eq(hashes.get(name,{}).get("sha256"),sha,"pinned_figure_source_"+name)
    vals={}
    for x in FIGURES:
        png=find_one(d,x+".png")
        svg=find_one(d,x+".svg")
        vals[x]={"png_width_height":check_png(png),"svg_bytes":check_svg(svg)}
    legend=find_one(d,"FIGURE_LEGENDS_AND_INTERPRETATION.md").read_text(encoding="utf8")
    for token in ("Figure 1", "Figure 2", "Figure 3", "Supplementary Table S1"):
        if token not in legend:
            raise ValueError("MISSING_LEGEND_"+token)
    s1=read_tsv(d,"Supplementary_Table_1_evidence_and_filters.tsv")
    eq(len(s1),16,"s1_candidate_count")
    by={x["gene_symbol"]:x for x in s1}
    eq(len(by),16,"s1_unique_symbols")
    if by["ERFE"]["scPred55_filter_status"].lower() not in ("filtered", "not_tested", "filtered_out", "not_tested_in_both"):
        # The source marks absent testing with explicit filter text; avoid guessing FDR.
        if not any(t in by["ERFE"]["scPred55_filter_status"].lower() for t in ("filter","not_test")):
            raise ValueError("ERFE_MISSING_EXPLICIT_FILTER_CAVEAT")
    return vals,s1


def check_strata(d):
    q=load_json(d,"ocular_stratified_sensitivity_manifest.json")
    eq(q["status"],"PASS_REAL_AGGREGATE_STRATIFIED_SENSITIVITY","strata_source")
    models={x["comparison"]:x for x in q["models"]}
    eq(models["score"]["shared_tested"],16050,"score_shared")
    eq(models["label"]["shared_tested"],13975,"label_shared")
    eq(models["label"]["same_direction_n"],12634,"label_concordance")
    eq(models["score"]["same_direction_n"],15851,"score_concordance")
    close(models["score"]["pearson_log2FC"],.99911597,"score_cor")
    close(models["label"]["pearson_log2FC"],.94990493,"label_cor")
    rows=read_tsv(d,"ocular_effect_stratum_sign_and_fdr.tsv")
    eq(len(rows),16,"stratum_16_rows")
    by={(x["comparison"],x["stratum_type"],x["stratum"]):x for x in rows}
    z=by[("label","baseline_abs_log2FC","|LFC| < 0.25")]
    eq(int(z["n_shared_tested"]),5349,"small_effect_gene_n")
    close(z["sign_concordance_percent"],78.258,"small_effect_sign_pct",.001)
    assert q["confidence_or_cell_annotation_sensitivity_is_not_independent_replication"]
    fig=check_png(find_one(d,"ocular_effect_stratified_sign_concordance.png"))
    return {"label_small_effect_sign_percent":float(z["sign_concordance_percent"]),
            "label_small_effect_n":5349,"supplemental_strata_figure_dimensions":fig}


def check_strain(d,s1):
    q=load_json(d,"GSE301517_strain_LOO_audit_manifest.json")
    eq(q["status"],"PASS_DESCRIPTIVE_7DAY_THREE_STRAIN_LEAVE_ONE_OUT_AUDIT","strain_source")
    p=q["panel_16"]
    eq((p["genes"],p["three_strains_unanimous_n"],p["7day_mean_matches_acute_n"],
        p["all_leave_one_strain_out_mean_match_acute_n"]), (16,7,10,8),"panel_three_criteria")
    emt=q["emt_194"]
    eq(emt["genes"],194,"emt_194")
    eq(emt["all_three_study_signed"],124,"emt_124")
    eq(emt["three_study_and_7day_3of3"],90,"emt_strict_90")
    eq(emt["7day_3of3_strains_same_as_acute"],103,"emt_strains_103")
    eq(emt["distribution_n_7day_strains_same_as_acute"],{"0":12,"1":25,"2":54,"3":103},"emt_distribution")
    rows=read_tsv(d,"GSE301517_16gene_three_strain_LOO_source_audit.tsv")
    eq(len(rows),16,"strain_panel_rows")
    by={x["gene_symbol"]:x for x in rows}
    eq(len(by),16,"strain_gene_symbol_duplicates")
    archived_s1={x["gene_symbol"]:x for x in s1}
    for gene, expected in CORE.items():
        x=by[gene]
        eq(x["three_of_three_strains_same_direction_as_acute"],"True","strain_unanimity_"+gene)
        eq(x["seven_day_all_three_leave_one_strain_out_same_direction_as_acute"],"True","strain_loo_"+gene)
        for col in ("TM62","TM70","TM72"):
            e=float(x[col+"_7day_log2ratio"])
            if e*expected<=0:
                raise ValueError("CORE_STRAIN_SIGN_MISMATCH_"+gene+"_"+col)
        if float(archived_s1[gene]["other_lab_500nM_7d_CPM_plus1_log2ratio"])*expected<=0:
            raise ValueError("ORIGINAL_FIGURE_SIGN_VS_NEW_7DAY_STRATUM_"+gene)
        close(x["seven_day_3strain_CPM_pseudocount_log2ratio"],
              archived_s1[gene]["other_lab_500nM_7d_CPM_plus1_log2ratio"],
              "mean_strain_vs_fig1_"+gene,.00004)
        if gene=="MMP2":
            close(x["TM72_7day_log2ratio"],-.16197,"mmp2_TM72",1e-5)
    ac=by["ACTA2"]
    eq(ac["seven_day_mean_matches_acute"],"False","acta2_context_discordance")
    if float(ac["TM62_7day_log2ratio"])<0 or float(ac["TM70_7day_log2ratio"])>=0:
        raise ValueError("ACTA2_SOURCE_STRAIN_CAVEAT_CHANGED")
    return {"three_study_signs":124,"strict_3_study_3strain_signs":90,
            "strict_percent":round(100*90/194,3),"core_strain_Loo":3}


def check_design(d):
    q=load_json(d,"GSE313736_technical_confounding_design_and_pool_decomposition.json")
    eq(q["status"],"PASS_AGGREGATE_IDENTIFIABILITY_AND_SELECTION_DECOMPOSITION","design_source")
    a=q["model_design_exact_rank"]
    for k,n in {
        "n_original_10x_libraries":18,"model_pool_plus_treatment_rank":10,
        "model_library_fixed_effect_plus_treatment_design_columns":19,
        "model_library_fixed_effect_plus_treatment_rank":18,
        "new_rank_from_treatment_after_unrestricted_library_effects":0,
    }.items():
        eq(a[k],n,"rank_"+k)
    if a["unrestricted_library_fixed_effect_model_drug_effect_identifiable"]:
        raise ValueError("FALSE_CAUSAL_IDENTIFIABILITY_CLAIM")
    p=q["pool_level_fraction_decomposition"]
    close(p["symmetric_within_pool_fraction_contribution_pp"],.8686384100016398,"within_pp")
    close(p["symmetric_pool_arm_weight_mix_contribution_pp"],1.5377008448624505,"mix_pp")
    close(p["pooled_BeamA_Dex_minus_Vehicle_percentage_points"],2.406339254864087,"total_pp")
    close(p["symmetric_within_pool_fraction_contribution_pp"]+
          p["symmetric_pool_arm_weight_mix_contribution_pp"],
          p["pooled_BeamA_Dex_minus_Vehicle_percentage_points"],"decomposition_sum",1e-8)
    close(p["paired_donor_BeamA_fraction_median_delta_pp"],-.421,"donor_median")
    if q["causal_dex_effect_vs_original_library_nuisance"] != "NOT_IDENTIFIABLE_FROM_SOURCE":
        raise ValueError("UNJUSTIFIED_CAUSAL_RESULT")
    return {"pool_within_pp":p["symmetric_within_pool_fraction_contribution_pp"],
            "pool_weight_pp":p["symmetric_pool_arm_weight_mix_contribution_pp"],
            "global_delta_pp":p["pooled_BeamA_Dex_minus_Vehicle_percentage_points"]}


def check_manuscript(path):
    text=path.read_text(encoding="utf8")
    required_sections={
      "Abstract": "## Abstract", "Source Methods": "### 2.1",
      "Identifiability Methods": "### 2.3a", "Stratified Results":"### 3.4a",
      "Technical Confounding Results":"### 3.4b", "Strain Results":"### 3.5a",
      "Rank interpretation":"### 4.3a", "Clinical boundary":"## 5. Conclusion",
      "Evidence GO/HOLD/STOP":"## 7. Explicit manuscript validation gates"}
    for label, phrase in required_sections.items():
        if phrase not in text:
            raise ValueError("MANUSCRIPT_MISSING_SECTION_"+label)
    ab=text.split("## Abstract",1)[1].split("## 1. Introduction",1)[0]
    for token in ("57 donor", "55", "13,975", "MYOC", "TGFBI", "MMP2",
                  "90 of 194", "18", "library"):
        if token not in ab:
            raise ValueError("ABSTRACT_STALE_OR_MISSING_"+token)
    for token in ("90/194", "1.537701", "19-column", "rank 18",
                  "NON-BLOCKING", "NOT", "CPM", "no clinical"):
        if token.lower() not in text.lower():
            raise ValueError("MANUSCRIPT_CLAIM_GATE_MISSING_"+token)
    if "80 gene" in ab:
        raise ValueError("AUTHOR80_MUST_NOT_BE_CLAIMED_IN_ABSTRACT")
    return {"required_section_count":len(required_sections),
            "abstract_updated_with_strict_90_and_alias":True}


def selftest():
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/"x.png"
        # Test CRC detection on a minimal valid PNG; size is deliberately too small.
        sig=b"\x89PNG\r\n\x1a\n"
        body=struct.pack(">IIBBBBB",2,2,8,2,0,0,0)
        def chunk(name, value):
            return struct.pack(">I",len(value))+name+value+struct.pack(">I",zlib.crc32(name+value)&0xffffffff)
        p.write_bytes(sig+chunk(b"IHDR",body)+chunk(b"IEND",b""))
        try:
            check_png(p)
            raise AssertionError("Bad small png accepted")
        except ValueError as e:
            assert "TOO_SMALL" in str(e)
        eq(18,18,"identity")
        try:
            eq(19,18,"drift")
            raise AssertionError("count drift not caught")
        except ValueError as e:
            assert "CLAIM_DRIFT" in str(e)
    print("OCULAR_THESIS_EVIDENCE_QA_SELFTEST_PASS")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--figures",type=Path)
    parser.add_argument("--strata",type=Path)
    parser.add_argument("--strain",type=Path)
    parser.add_argument("--design",type=Path)
    parser.add_argument("--manuscript",type=Path)
    parser.add_argument("--out",type=Path)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    if args.self_test:
        selftest();return
    if not all((args.figures,args.strata,args.strain,args.design,args.manuscript,args.out)):
        parser.error("all four actual archived source directories and manuscript/out paths are required")
    fig,s1=figure_integrity(args.figures)
    strata=check_strata(args.strata)
    strain=check_strain(args.strain,s1)
    design=check_design(args.design)
    manuscript=check_manuscript(args.manuscript)
    artifact_hashes={}
    for key, dirpath in (("figures",args.figures),("strata",args.strata),
                          ("strain",args.strain),("design",args.design)):
        names={}
        for f in sorted(dirpath.rglob("*")):
            if f.is_file():
                names[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
        artifact_hashes[key]=names
    report={
        "status":"PASS_ORIGINAL_FIGURES_LATEST_REAL_CI_AND_MANUSCRIPT_INTEGRITY",
        "source_actions":RUNS,
        "source_file_sha256":artifact_hashes,
        "figure_png_and_svg_integrity":fig,
        "stratified_effects":strata,
        "independent_seven_day_strain_directions":strain,
        "design_and_capture_weight":design,
        "manuscript_structure_claims":manuscript,
        "author_gene80_verified":False,"author_contact_requested":False,
        "original_library_drug_effect_separably_identified":False,
        "peer_review_or_original_experiment_replication":False,
        "article_ready":False,
        "limitations":"CI confirms the alignment of frozen public artifacts, NOT novelty, blind image layout quality, source-independent replication, protein function, clinical DEX-IOP validation or peer review."
    }
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"OCULAR_FINAL_FIGURE_SOURCE_CLAIM_QA.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    lines=[
        "# GSE313736 ocular manuscript: actual four-archive source/figure QA",
        "",
        "**Computational source and claim integrity: PASS. Journal submission readiness: HOLD.**",
        "",
        "| Audit | Real archived check | Verdict |",
        "|---|---|---|",
        "| Main Fig 1–3 | 3 PNG CRC and dimensions, 3 SVG XML, legends and 16-gene S1 | PASS asset integrity only |",
        "| Fig 3 vs strict follow-up | source 124/194, strict 3-strain criterion 90/194 | PASS with supplement disclosure |",
        "| Stratum supplement | 13,975 common genes; 78.258% sign concordance in 5,349 weak-effect genes | PASS |",
        "| Seven-day independent study | 3 strain-level signs for MYOC/TGFBI/MMP2; no gene-wise 7d FDR | PASS descriptive |",
        "| Source treatment/lib alias | 19 columns/rank 18, DEX incremental rank 0 | HOLD causality |",
        "| Cell-weight decomposition | +0.868638 within +1.537701 mix = +2.406339 pp | PASS arithmetic, not causal |",
        "| Submitted/peer-reviewed manuscript | Not established by the artifact audit | HOLD |",
        "",
        "**DO NOT** call the external 3-strain technical-replicate count seven biological donors, say 124/194 passed per-strain unanimity, or infer that 63.9% of DEX signal is batch artifact.",
        "The figure PNG CRC/SVG parse does not establish visual absence of overlaps; editorial visual QA remains a separate manual gate.",
        "",
        "Exact source artifact IDs, SHA-256 filenames and counts are in the JSON."
    ]
    (args.out/"OCULAR_FINAL_FIGURE_SOURCE_CLAIM_QA.md").write_text(
        "\n".join(lines)+"\n",encoding="utf8")
    print("OCULAR_THESIS_FOUR_ARCHIVE_QA_PASS",json.dumps({
        "main_figure_count":len(FIGURES),"supplementary_figures":1,
        "emt_mean":124,"emt_strict":90,"rank_library_treatment":18,
        "pool_weights_pp":design["pool_weight_pp"],
        "manuscript_checked":True,
        "submission_ready":False},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
