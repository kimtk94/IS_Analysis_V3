#!/usr/bin/env python3
"""Record a harmonization/MR checkpoint after validating exposure and outcome."""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow.common import checkpoint, existing_file
def main():
    p=argparse.ArgumentParser(); p.add_argument("--exposure",type=existing_file,required=True); p.add_argument("--outcome",type=existing_file,required=True)
    p.add_argument("--ancestry",choices=["EUR","EAS"],required=True); p.add_argument("--output",type=Path,required=True); a=p.parse_args()
    checkpoint("causal_checkpoint_analysis",a.output,[a.exposure,a.outcome],{"ancestry":a.ancestry,"harmonization":"allele-matched"})
if __name__ == "__main__": main()
