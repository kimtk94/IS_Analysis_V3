#!/usr/bin/env python3
from __future__ import annotations
import json,os
from pathlib import Path
ROOT=Path(os.environ.get('IS_ANALYSIS_ROOT','/srv/is-analysis'))
checks=[
 ('stage2_shortlist',ROOT/'results/metabolic_resilience/stage2_gwas/final_shortlist/STAGE2_FINAL_SCREENING_SHORTLIST.tsv',True),
 ('stage3_ld_reference',ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C0R_LD_MAPPING_SUMMARY.json',True),
 ('stage3_clump',ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C1_LD_CLUMP_SUMMARY.json',True),
 ('stage3_frozen_iv',ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C2_INSTRUMENT_FREEZE.json',True),
 ('stage3_mr',ROOT/'results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D3_MR_RESULTS.tsv',True),
 ('stage3_multiplicity_overlap',ROOT/'results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D6_MULTIPLE_TESTING_AND_OVERLAP.json',True),
 ('stage3_coloc',ROOT/'results/metabolic_resilience/stage3_coloc/STAGE3E3_COLOC_SUMMARY.tsv',True),
 ('stage3_multisignal_gate',ROOT/'results/metabolic_resilience/stage3_coloc/STAGE3E4_MULTISIGNAL_DECISION.json',True),
 ('variant_artifact_audit',ROOT/'results/metabolic_resilience/stage3_variant_artifact/STAGE3F0_VARIANT_ARTIFACT_SUMMARY.json',True),
 ('eas_registry',ROOT/'results/metabolic_resilience/stage4_eas/STAGE4D_EAS_RESOURCE_REGISTRY_SUMMARY.json',False),
 ('koges_contract',ROOT/'results/metabolic_resilience/stage5_koges/STAGE5D_ANALYSIS_CONTRACT.json',False),
 ('functional_registry',ROOT/'results/metabolic_resilience/stage6_functional/STAGE6A_FUNCTIONAL_ANNOTATION_REGISTRY.json',False),
]
rows=[]; blocking=0
for name,path,required in checks:
    exists=path.exists() and path.stat().st_size>0
    status='EXISTS' if exists else 'MISSING'
    if required and not exists: blocking+=1
    rows.append({'name':name,'path':str(path),'required_for_core_claims':required,'status':status})
payload={'status':'PASS_CORE_READY' if blocking==0 else 'HOLD_CORE_INCOMPLETE','blocking_missing':blocking,'checks':rows,'rule':'Final manuscript-level causal claims remain HOLD until every required core gate exists and any REVIEW/HOLD statuses inside those outputs are resolved.'}
out=ROOT/'results/metabolic_resilience/MASTER_ANALYSIS_GATE.json'; out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(payload,indent=2)+'\n'); print(json.dumps(payload,indent=2)); raise SystemExit(0 if blocking==0 else 2)
