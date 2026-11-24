#!/usr/bin/env python3
"""Inventory core output files without certifying scientific readiness.

This command does not parse or validate the output schemas or their internal
REVIEW/HOLD decisions. Even a complete inventory remains scientifically HOLD.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

CHECKS = [
 ('stage2_shortlist','results/metabolic_resilience/stage2_gwas/final_shortlist/STAGE2_FINAL_SCREENING_SHORTLIST.tsv',True),
 ('stage3_ld_reference','results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C0R_LD_MAPPING_SUMMARY.json',True),
 ('stage3_clump','results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C1_LD_CLUMP_SUMMARY.json',True),
 ('stage3_frozen_iv','results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C2_INSTRUMENT_FREEZE.json',True),
 ('stage3_mr','results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D3_MR_RESULTS.tsv',True),
 ('stage3_multiplicity_overlap','results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D6_MULTIPLE_TESTING_AND_OVERLAP.json',True),
 ('stage3_coloc','results/metabolic_resilience/stage3_coloc/STAGE3E3_COLOC_SUMMARY.tsv',True),
 ('stage3_multisignal_gate','results/metabolic_resilience/stage3_coloc/STAGE3E4_MULTISIGNAL_DECISION.json',True),
 ('variant_artifact_audit','results/metabolic_resilience/stage3_variant_artifact/STAGE3F0_VARIANT_ARTIFACT_SUMMARY.json',True),
 ('eas_registry','results/metabolic_resilience/stage4_eas/STAGE4D_EAS_RESOURCE_REGISTRY_SUMMARY.json',False),
 ('koges_contract','results/metabolic_resilience/stage5_koges/STAGE5D_ANALYSIS_CONTRACT.json',False),
 ('functional_registry','results/metabolic_resilience/stage6_functional/STAGE6A_FUNCTIONAL_ANNOTATION_REGISTRY.json',False),
]


def inventory(root: Path) -> dict:
    rows = []
    blocking = 0
    for name, relative_path, required in CHECKS:
        path = root / relative_path
        exists = path.is_file() and path.stat().st_size > 0
        if required and not exists:
            blocking += 1
        rows.append({
            'name': name,
            'path': str(path),
            'required_for_core_claims': required,
            'status': 'EXISTS' if exists else 'MISSING',
        })
    return {
        'status': ('FILES_PRESENT_NOT_SCIENTIFICALLY_VALIDATED' if blocking == 0
                   else 'HOLD_CORE_INCOMPLETE'),
        'scientific_status': 'HOLD',
        'scientific_validation_performed': False,
        'validation_scope': 'nonempty_regular_files_only',
        'blocking_missing': blocking,
        'checks': rows,
        'rule': (
            'File presence does not establish scientific readiness. Final '
            'manuscript-level causal claims remain HOLD until required output '
            'schemas, scientific results, and internal REVIEW/HOLD decisions '
            'have been validated and resolved.'
        ),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--root', type=Path,
        default=Path(os.environ.get('IS_ANALYSIS_ROOT', '/srv/is-analysis')),
        help='Analysis root (default: IS_ANALYSIS_ROOT or /srv/is-analysis).',
    )
    args = parser.parse_args(argv)
    payload = inventory(args.root)
    out = args.root / 'results/metabolic_resilience/MASTER_ANALYSIS_GATE.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, indent=2) + '\n'
    out.write_text(text, encoding='utf-8')
    print(text, end='')
    return 2  # File inventory alone never clears the scientific HOLD.


if __name__ == '__main__':
    raise SystemExit(main())
