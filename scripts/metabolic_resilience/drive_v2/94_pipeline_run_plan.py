#!/usr/bin/env python3
from __future__ import annotations
import json,os
from pathlib import Path
ROOT=Path(os.environ.get('IS_ANALYSIS_ROOT','/srv/is-analysis'))
REPO=Path(os.environ.get('IS_ANALYSIS_REPO',str(ROOT/'IS_Analysis_V3')))
CODE=REPO/'scripts/metabolic_resilience/drive_v2'
steps=[
 ('31','bash', '31_stage3c0r_build_1000g_eur_references.sh','results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C0R_LD_MAPPING_SUMMARY.json'),
 ('32','python3','32_stage3c1_ld_clump.py','results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C1_LD_CLUMP_SUMMARY.json'),
 ('33','python3','33_stage3c2_freeze_instruments.py','results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C2_INSTRUMENT_FREEZE.json'),
 ('44','bash','44_stage3d4_run_confirmatory_mr.sh','results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D3_MR_RESULTS.tsv'),
 ('46','python3','46_stage3d6_multiple_testing_overlap.py','results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D6_MULTIPLE_TESTING_AND_OVERLAP.json'),
 ('47','python3','47_stage3d7_steiger_sensitivity.py','results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D7_STEIGER_SENSITIVITY.json'),
 ('53','python3','53_stage3e3_summarize_coloc.py','results/metabolic_resilience/stage3_coloc/STAGE3E3_COLOC_SUMMARY.tsv'),
 ('54','python3','54_stage3e4_multisignal_coloc_gate.py','results/metabolic_resilience/stage3_coloc/STAGE3E4_MULTISIGNAL_DECISION.json'),
 ('55','python3','55_stage3f0_variant_artifact_audit.py','results/metabolic_resilience/stage3_variant_artifact/STAGE3F0_VARIANT_ARTIFACT_SUMMARY.json'),
 ('65','python3','65_stage4f_eas_ld_reference_plan.py','results/metabolic_resilience/stage4_eas/STAGE4F_EAS_LD_REFERENCE_PLAN.json'),
 ('73','python3','73_stage5d_write_analysis_contract.py','results/metabolic_resilience/stage5_koges/STAGE5D_ANALYSIS_CONTRACT.json'),
 ('80','python3','80_stage6a_functional_annotation_scaffold.py','results/metabolic_resilience/stage6_functional/STAGE6A_FUNCTIONAL_ANNOTATION_REGISTRY.json'),
 ('91','python3','91_master_analysis_gate.py','results/metabolic_resilience/MASTER_ANALYSIS_GATE.json'),
 ('92','bash','92_environment_snapshot.sh','results/metabolic_resilience/reproducibility'),
]
rows=[]
for order,exe,script,output in steps:
    p=ROOT/output; exists=p.exists() and (p.is_dir() or p.stat().st_size>0)
    rows.append({'order':order,'script':script,'command':f'{exe} {CODE/script}','expected_output':str(p),'state':'DONE' if exists else 'PENDING'})
payload={'steps':rows,'next':[r for r in rows if r['state']=='PENDING'][:1]}; print(json.dumps(payload,indent=2))
