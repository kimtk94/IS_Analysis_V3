#!/usr/bin/env python3
"""Validate regional input and record ancestry-matched LD/fine-map provenance."""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow.common import checkpoint, existing_file
def main():
    p=argparse.ArgumentParser(); p.add_argument("--regional-summary",type=existing_file,required=True)
    p.add_argument("--ld-reference",type=existing_file,required=True); p.add_argument("--ancestry",choices=["EUR","EAS"],required=True)
    p.add_argument("--output",required=True,type=Path); a=p.parse_args()
    checkpoint("ancestry_ld_finemap",a.output,[a.regional_summary,a.ld_reference],{"ancestry":a.ancestry})
if __name__ == "__main__": main()
