#!/usr/bin/env python3
"""Ocular P9: auditable gene panel from immutable public P5b/P7a/P8a aggregates.

No raw donor matrices; published cohort comparisons are descriptive. Never call
cross-study concordance proof of dose effect, protein change, or IOP prediction.
"""
import argparse,csv,gzip,json,math
from pathlib import Path
from collections import Counter

PRIMARY_REFERENCE=["MYOC","TGFBI","MMP2","ACTA2"]
CONTEXTUAL=["SDC1","EGR3","LPAR1","TFPI"]
ECM_SECONDARY=["LUM","IGFBP3","VCAM1","NTM"]
TECHNICAL_ONLY=["FKBP5"]
EXPLORATORY_LOW=["CHAC1","ERFE","PARM1"]
ALIASES={"CTGF":"CCN2"}
EMT_SIGNIFICANT_BOTH=49
EMT_SIGNIFICANT_SAME=46
EXPECTED_OPPOSITE={"ERFE","CHAC1","PARM1","ADAMTS15","EGR3","LPAR1","ABHD3","DDIT4L","TFPI","STK17B","RAB32","SDC1"}
def readtsv(p):
    with p.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def readgz(p):
    with gzip.open(p,"rt",encoding="utf8",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def number(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None
def bools(v):
    if v in (True,"True","TRUE","true"):return True
    if v in (False,"False","FALSE","false"):return False
    raise ValueError("INVALID_BOOLEAN_VALUE "+str(v))
def write_json(p,d):
    p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")
def write_tsv(p,records):
    if not records:raise ValueError("EMPTY_OUTPUT")
    with p.open("w",encoding="utf8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t");wr.writeheader();wr.writerows(records)
def summarize(o):
    a=readtsv(o.p7_panel)
    b=readtsv(o.p7_emt)
    c=readtsv(o.p8_both)
    d=readgz(o.primary)
    if (len(a),len(b),len(c),len(d))!=(24,194,49,16050):
        raise ValueError("VERIFIED_P5_P7_P8_ARTIFACT_DIMENSION_DRIFT")
    bypanel={x["symbol"]:x for x in a}
    byemt={x["gene_symbol"]:x for x in b}
    byboth={x["gene_symbol"]:x for x in c}
    if len(bypanel)!=24 or len(byemt)!=194 or len(byboth)!=49:
        raise ValueError("DUPLICATED_SYMBOLS_IN_SOURCE")
    if {x["gene_symbol"] for x in c}!={x["gene_symbol"] for x in b if bools(x["both_study_FDR_lt_0_05"])}:
        raise ValueError("P8_BOTH_FDR_49_VS_P7_SOURCE_MISMATCH")
    if sum(bools(x["sign_agrees"]) for x in c)!=EMT_SIGNIFICANT_SAME:
        raise ValueError("49_SIGNIFICANT_EMT_WITH_46_SAME_EXPECTED")
    allprimary={x["gene_symbol"]:x for x in d}
    if set(PRIMARY_REFERENCE+CONTEXTUAL+ECM_SECONDARY+TECHNICAL_ONLY+EXPLORATORY_LOW)-set(allprimary):
        raise ValueError("P9_REQUIRED_GENE_MISSING_IN_PRIMARY")
    for sym,x in bypanel.items():
        if sym=="CTGF":
            if bools(x["present_both"]):
                raise ValueError("CTGF_MUST_BE_ALIAS_NOT_SECOND_LOCUS")
            continue
        if not bools(x["present_both"]):raise ValueError("P7_PANEL_MISSING "+sym)
        p=number(x["primary_lfc"])
        q=number(allprimary[sym]["logFC"])
        if p is None or q is None or abs(p-q)>0.00002:
            raise ValueError("PRIMARY_PANEL_LOGFC_VS_P5B_DRIFT "+sym)
    # 12 originally selected opposite-sign genes should not become validation positives.
    if not EXPECTED_OPPOSITE.issubset(bypanel):raise ValueError("P7_ORIGINAL_12_DRIFT")
    if any(bools(bypanel[s]["same_signed_effect"]) for s in EXPECTED_OPPOSITE):
        raise ValueError("FALSE_CROSS_STUDY_12_REPLICATION_LABEL")
    for sym in PRIMARY_REFERENCE:
        x=bypanel[sym]
        if not bools(x["same_signed_effect"]) or not number(x["primary_FDR"])<.05 or not number(x["external_FDR"])<.05:
            raise ValueError("CORE_4_CROSS_STUDY_FDR_SIGN_INVALID "+sym)
        if min(abs(number(x["primary_lfc"])),abs(number(x["external_lfc"])))<0.3:
            raise ValueError("CORE_4_CROSS_STUDY_MIN_EFFECT_CHANGED "+sym)
    for sym in ECM_SECONDARY:
        if sym not in byboth or not bools(byboth[sym]["sign_agrees"]):
            raise ValueError("ECM_SECONDARY_49_LIST_MISMATCH "+sym)
    return bypanel,byemt,byboth,allprimary

def collect_panel(a):
    pan,emt,both,prim=summarize(a)
    selected=(PRIMARY_REFERENCE+CONTEXTUAL+ECM_SECONDARY+TECHNICAL_ONLY+EXPLORATORY_LOW)
    out=[]
    for symbol in selected:
        row=pan.get(symbol)
        x=prim[symbol]
        isprimary=number(x["logFC"])
        express=number(x["logCPM"])
        fdr=number(x["FDR"])
        if any(v is None for v in (isprimary,express,fdr)):
            raise ValueError("SOURCE_EXPRESSION_NUMERIC_MISSING "+symbol)
        if row:
            ext=number(row["external_lfc"])
            external_fdr=number(row["external_FDR"])
            primary=number(row["primary_lfc"])
            sign_agree=bools(row["same_signed_effect"])
            if primary is None or abs(primary-isprimary)>0.00002:
                raise ValueError("PRIMARY_SOURCE_EFFECT_MISMATCH "+symbol)
        elif symbol in both:
            ext=number(both[symbol]["external_log2fc"])
            external_fdr=number(both[symbol]["external_FDR"])
            sign_agree=bools(both[symbol]["sign_agrees"])
        else:
            # FKBP5 positive reference in primary source, but not externally
            # quantified in our matched P7/P8 artifact; do not invent value.
            ext=None;external_fdr=None;sign_agree=None
        bothsig=(fdr<.05 and external_fdr is not None and external_fdr<.05)
        if symbol in PRIMARY_REFERENCE:
            group="CORE_CROSS_STUDY_CONCORDANT"
            biological_role="MYOC_STEROID_RESPONSE_REFERENCE" if symbol=="MYOC" else "ECM_CELL_STATE_TRANSCRIPT"
            if not(bothsig and sign_agree):raise ValueError("CORE_STATUS_DRIFT")
        elif symbol in CONTEXTUAL:
            group="SECONDARY_DISCORDANT_CONTEXT_TEST"
            biological_role="EXPLORATORY_CONDITION_SENSITIVE"
            if not(bothsig and sign_agree is False):raise ValueError("CONTEXTUAL_SIGN_DRIFT "+symbol)
        elif symbol in ECM_SECONDARY:
            group="SECONDARY_CONCORDANT_ECM_PRIORITY"
            biological_role="EXPLORATORY_ECM_OR_CELL_ADHESION"
            if not(bothsig and sign_agree):raise ValueError("EXTERNAL_ECM_SIGN_FDR_DRIFT "+symbol)
        elif symbol in TECHNICAL_ONLY:
            group="PRIMARY_ONLY_GLUCORTICOID_REFERENCE"
            biological_role="GR_SIGNALING_ASSAY_CHECK"
            if ext is not None:raise ValueError("PRIMARY_ONLY_EXTERNAL_VALUE_INVENTED")
        else:
            group="LOW_EXPRESSION_EXPLORATORY_NOT_CORE"
            biological_role="LOW_ABUNDANCE_SENSITIVE"
            if not express<2:raise ValueError("LOW_EXPRESSION_TIER_DRIFT "+symbol)
        out.append({"gene_symbol":symbol,"gene_id":x["gene_id"],
            "panel_group":group,"biological_role":biological_role,
            "primary_500nM_log2fc":round(isprimary,6),
            "primary_mean_logCPM":round(express,6),
            "primary_BH_FDR":fdr,
            "external_100nM_log2fc":round(ext,6) if ext is not None else "",
            "external_BH_FDR":external_fdr if external_fdr is not None else "",
            "same_direction_in_two_nonmatched_experiments":sign_agree if sign_agree is not None else "",
            "external_numeric_crosswalk_verified":ext is not None,
            "both_studies_BH_FDR_lt_0_05":bothsig,
            "EMT_member":symbol in emt,
            "core_qpcr_recommended":symbol in PRIMARY_REFERENCE,
            "orthogonal_protein_assay_requires_feasibility":symbol in ("MYOC","TGFBI","MMP2","ACTA2"),
            "validated_steroid_induced_delta_IOP_predictor":False})
    if len(out)!=16 or len({x["gene_id"] for x in out})!=16:
        raise ValueError("SELECTED_PANEL_EXPECTED_16_DISTINCT_GENES")
    write_tsv(a.out/"P9a_16_gene_evidence_tiered_validation_panel.tsv",out)
    statuses=Counter(x["panel_group"] for x in out)
    data={"status":"PASS_P9_CORE_CANDIDATE_EVIDENCE_RECONCILED",
      "n_distinct_panel_genes":len(out),"group_counts":dict(statuses),
      "core_four":PRIMARY_REFERENCE,
      "contextual_discordant_four":CONTEXTUAL,
      "secondary_ECM_four":ECM_SECONDARY,
      "primary_only_control":TECHNICAL_ONLY,
      "low_abundance_exploratory_three":EXPLORATORY_LOW,
      "original_opposite_12_candidate_set_is_posthoc":True,
      "CTGF_alias_of_CCN2":True,
      "CTGF_CCn2_is_single_locus_not_missing_independent_gene":True,
      "external_FKBP5_gene_level_value_not_present_in_supplied_P7_P8_crosswalk":True,
      "clinical_predictive_or_protein_functional_validation":False,
      "primary_treatment_and_library_batch_confounded":True,
      "interpretation":"Only four CORE genes have same-sign and gene-wise BH-FDR<.05 in both different 48h TM studies; neither dose-specific randomized contrast nor ECM protein changes were validated."}
    write_json(a.out/"P9a_core_panel_reconciliation_gate.json",data)
    print("P9_PANEL",json.dumps(data,ensure_ascii=False),flush=True)

def audit(a):
    pan,emt,both,prim=summarize(a)
    selected=Counter()
    # Cross-study cell state marker response / interpretability gate
    markers=["MYOC","TGFBI","MMP2","ACTA2","FN1","COL1A1","CCN2","THBS1","TIMP1","SDC1"]
    rows=[]
    for sym in markers:
        x=pan[sym]
        p=number(x["primary_FDR"]);q=number(x["external_FDR"])
        ext=number(x["external_lfc"]);pfc=number(x["primary_lfc"])
        p_active=p is not None and p<.05
        q_active=q is not None and q<.05
        same=(pfc*ext>0) if pfc is not None and ext is not None else False
        status=("TWO_STUDIES_SIGNED_AND_FDR" if p_active and q_active and same
            else "TWO_STUDIES_SIGNIFICANT_OPPOSITE" if p_active and q_active
            else "PRIMARY_ONLY_SIGNIFICANT" if p_active
            else "NOT_SIGNIFICANT_IN_EITHER_SOURCE")
        selected[status]+=1
        rows.append({"gene":sym,"primary_fc":pfc if pfc is not None else "",
          "external_fc":ext if ext is not None else "",
          "source_primary_BH_FDR":p if p is not None else "",
          "external_BH_FDR":q if q is not None else "",
          "status":status,"gene_is_not_biological_replicate":True,
          "unknown_protein_ECM_state":True})
    write_tsv(a.out/"P9b_preselected_TM_marker_assay_risk_register.tsv",rows)
    if not all(x["status"]=="TWO_STUDIES_SIGNED_AND_FDR" for x in rows if x["gene"] in PRIMARY_REFERENCE):
        raise ValueError("CORE_PHENOTYPE_REFERENCE_CHANGED")
    if not any(x["gene"]=="SDC1" and x["status"]=="TWO_STUDIES_SIGNIFICANT_OPPOSITE" for x in rows):
        raise ValueError("ORIGINAL_DISCORDANT_NEGATIVE_CONTROL_LOST")
    if not any(x["gene"]=="CCN2" and x["status"]=="PRIMARY_ONLY_SIGNIFICANT" for x in rows):
        raise ValueError("CCN2_EXPRESSION_VS_EXTERNAL_NON_SIG_STATUS_CHANGED")
    # CTGF and CCN2 are aliases, not an additional separate RNA expression target.
    aliases=[]
    for label,canonical in ALIASES.items():
        arow=pan.get(label);brow=pan.get(canonical)
        if arow is None or brow is None:raise ValueError("ALIAS_EVIDENCE_MISSING")
        aliases.append({"input_alias":label,"canonical_gene_symbol":canonical,
          "alias_would_double_count_distinct_genes":True,
          "source_CTGF_entry_has_actual_quantified_gene":bools(arow["present_both"]),
          "source_CCN2_entry_has_actual_quantified_gene":bools(brow["present_both"]),
          "recommended_qPCR_label":canonical})
    write_tsv(a.out/"P9b_alias_normalization_CTGF_CCN2.tsv",aliases)
    summary={"status":"PASS_P9_MARKER_POSITIVE_NEGATIVE_AND_ALIAS_AUDIT",
      "marker_panel_tested":len(markers),
      "marker_status_counts":dict(selected),
      "MYOC_is_dose_phenotype_pos_control_not_proof_of_batch_balance":True,
      "TGFBI_MMP2_ACTA2_are_pathway_readouts_not_proof_of_functional_ecm":True,
      "SDC1_is_direction_discordant_warning":True,
      "CCN2_and_CTGF_are_alias_same_gene":True,
      "require_housekeeping_stability_on_proposed_sample_conditions":True,
      "feature_selection_was_posthoc":True}
    write_json(a.out/"P9b_assay_risk_gate.json",summary)
    print("P9_AUDIT",json.dumps(summary,ensure_ascii=False),flush=True)

def protocol(a):
    pan,emt,both,prim=summarize(a)
    # This is prospective protocol DESIGN only, not a completed wet-lab experiment.
    blocks=[]
    for objective,doses,times,readouts in (
      ("PRIMARY_RANDOMIZED_48H",["0 nM vehicle","100 nM DEX","500 nM DEX"],["48 h"],
       ["RT-qPCR: MYOC,TGFBI,MMP2,ACTA2","ELISA/Western: MYOC,TGFBI,MMP2 when feasible","MMP2 zymography, extracellular matrix imaging and contractility/outflow surrogate if feasible"]),
      ("SECONDARY_TIME_COURSE_OPTIONAL",["0 nM vehicle","100 nM DEX","500 nM DEX"],["1 h","6 h","48 h"],
       ["RT-qPCR: four core + EGR3,LPAR1,SDC1,TFPI","MMP2 activity and selected ECM protein if feasible"])):
        blocks.append({"objective":objective,"doses":doses,"timepoints":times,
          "treatment_time_conditions":len(doses)*len(times),
          "assays":readouts,"donor_is_biological_unit":True,
          "plate_library_batch_must_contain_randomized_dose_conditions":True,
          "number_of_donors_to_fix_after_pilot_variance_not_invented":True})
    plan={
      "status":"PASS_EXPLICIT_FUTURE_PROTOCOL_DESIGN_NOT_EXPERIMENT_PERFORMED",
      "core_genes":PRIMARY_REFERENCE,
      "context_study_genes":CONTEXTUAL,
      "technical_pilot_reference":"FKBP5",
      "negative_or_contextual_controls":["SDC1","FN1","COL1A1","CCN2"],
      "primary_hypothesis":"Dexamethasone changes TM MYOC expression and ECM-related TGFBI/MMP2/ACTA2 consistently in independently randomized donor cultures",
      "biological_unit":"one independent human donor TM culture",
      "dose_factor_levels":3,
      "primary_timepoint_hours":48,
      "randomization":"each donor split into dose arms; balance treatment order and processing plate within donor",
      "technical_replicates_are_not_biological_replicates":True,
      "housekeeping_reference_genes_need_assay_stability_validation":True,
      "negative_control":"vehicle matched for solvent and incubation time; mycoplasma and endotoxin status documented",
      "positive_reference":"MYOC and FKBP5 response, but positive marker does not eliminate all batch/culture confounding",
      "planned_statistical_contrast":"donor-blocked dose effect with 0/100/500nM as randomized fixed factor; testing within 48h primary design; correct multiplicity",
      "dose_reversal_test":"pre-specify 100 vs 500 nM contrast and direction interaction using common culture/time; do not infer from cross-study discordance",
      "primary_endpoint":"relative RNA expression by validated qPCR; no patient steroid delta-IOP endpoint",
      "functional_secondary_endpoints":["secreted MYOC","MMP2 activity","ECM deposition","collagen/fibronectin imaging","contractility/outflow surrogate, if technically available"],
      "pilot_required_for_variance_and_power":True,
      "sample_size_claim":"NOT_DETERMINED_NO_VARIANCE_DATA",
      "analytical_blinding":"blind image/activity quantification to assigned condition",
      "condition_blocks":blocks,
      "hard_stop_without_controls":["no donor-matched vehicle","no documented microbial culture QC","separate non-randomized libraries fully aliased with condition","no pre-specified multiple-testing adjustment"],
      "required_provenance":["sample donor codes pseudonymized","dose/time/solvent","passage and viability","mycoplasma test","library/batch/order","raw qPCR Ct with reference-gene QC","protein assay controls"],
      "not_clinical_IOP_prediction":True,
      "published_lab_precedent":[
       {"pmid":"27939447","url":"https://pubmed.ncbi.nlm.nih.gov/27939447/","claims":"MYOC secretion up, MMP2 secretion/activity down after 100nM DEX in HTM cultures"},
       {"pmid":"11726619","url":"https://pubmed.ncbi.nlm.nih.gov/11726619/","claims":"MYOC glucocorticoid responsiveness can be measured by qPCR and immunoblot"},
       {"pmid":"41636427","url":"https://pubmed.ncbi.nlm.nih.gov/41636427/","claims":"independent TM 100nM DEX RNAseq at 1h,6h,48h"}]}
    write_json(a.out/"P9c_prospective_donor_randomized_TM_validation_design.json",plan)
    conditions=[{"study":"PRIMARY_RANDOMIZED_48H","timepoint_h":48,"dose_nM":dose,
       "vehicle_matched":dose==0,"donor_assigned_all_three_doses":True,
       "planned_not_executed":True} for dose in (0,100,500)]
    write_tsv(a.out/"P9c_minimal_prospective_48h_condition_matrix.tsv",conditions)
    print("P9_PROTOCOL",json.dumps({"status":plan["status"],"conditions":len(conditions),"positive":PRIMARY_REFERENCE},ensure_ascii=False),flush=True)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--mode",required=True,choices=("panel","audit","protocol"))
    for name in ("p7-panel","p7-emt","p8-both","primary"):
        parser.add_argument("--"+name,required=True,type=Path)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    if a.mode=="panel":collect_panel(a)
    elif a.mode=="audit":audit(a)
    else:protocol(a)
if __name__=="__main__":main()
