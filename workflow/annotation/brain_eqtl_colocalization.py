#!/usr/bin/env python3
"""Validate regional pQTL/eQTL inputs and record coloc checkpoint provenance."""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from workflow.common import checkpoint, existing_file
def main():
    p=argparse.ArgumentParser(); p.add_argument("--pqtl",type=existing_file,required=True); p.add_argument("--eqtl",type=existing_file,required=True)
    p.add_argument("--gene",required=True); p.add_argument("--output",type=Path,required=True); a=p.parse_args()
    checkpoint("brain_eqtl_colocalization",a.output,[a.pqtl,a.eqtl],{"gene":a.gene.upper(),"regional_full_summary_required":True})
if __name__ == "__main__": main()
