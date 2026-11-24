#!/usr/bin/env python3
"""Variant-allele-specific VEP audit for fixed seven discovery sentinels."""
import argparse
import hashlib
import json
from pathlib import Path
from ckd_korean_validation import anchors, get, write

def run(root,out):
    panel=anchors(root)
    results=[]
    for rs,rr in panel.items():
        a=rr[0]
        variant=a["sentinel"].split(":")
        chrom,pos,ref,alt=variant[0].removeprefix("chr"),int(variant[1]),variant[2],variant[3]
        url="https://rest.ensembl.org/vep/human/id/"+rs+"?canonical=1&hgvs=1&mane=1"
        r={"rsid":rs,"source_sentinel":a["sentinel"],"url":url,"status":"HOLD"}
        try:
            raw=get(url,4000000)
            data=json.loads(raw)
            r["sha256"]=hashlib.sha256(raw).hexdigest()
            if not isinstance(data,list):raise ValueError("Unexpected VEP response")
            hits=[v for v in data if str(v.get("seq_region_name"))==chrom and v.get("start")==v.get("end")==pos and v.get("assembly_name")=="GRCh38" and v.get("allele_string","").split("/")[0]==ref and alt in v.get("allele_string","").split("/")[1:]]
            if len(hits)!=1:raise ValueError("No unique GRCh38 matching allele record")
            v=hits[0]
            consequences=[]
            for x in v.get("transcript_consequences",[]):
                if x.get("variant_allele")!=alt:continue
                keys=("gene_symbol","gene_id","transcript_id","biotype","canonical","mane_select","mane_plus_clinical","consequence_terms","hgvsc","hgvsp","amino_acids","protein_start","variant_allele")
                consequences.append({k:x[k] for k in keys if k in x})
            flags={"missense_variant","stop_gained","stop_lost","start_lost","frameshift_variant","inframe_insertion","inframe_deletion","splice_acceptor_variant","splice_donor_variant"}
            coding=[x for x in consequences if flags.intersection(x.get("consequence_terms",[]))]
            r.update(status="MATCHED_GRCH38_ALLELE",all_alleles=v["allele_string"],query_record_most_severe=v.get("most_severe_consequence"),
                     alt_allele_terms=sorted({t for x in consequences for t in x.get("consequence_terms",[])}),
                     canonical_or_mane=[x for x in consequences if x.get("canonical")==1 or x.get("mane_select") or x.get("mane_plus_clinical")],
                     protein_or_splice_altering_transcripts=coding,
                     note="VEP annotation, not a linkage-disequilibrium PAV screen or proof of assay interference.")
        except (OSError,ValueError,KeyError,TypeError) as e:r["reason"]=str(e)[:300]
        results.append(r)
    write(out,{"status":"ALLELE_SPECIFIC_VEP_ANNOTATION","results":results,"limit":"Current Ensembl annotation can include multiple alternate alleles/transcripts. Only exact study ALT is assessed. No functional or causal inference."})

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.root,a.out)
