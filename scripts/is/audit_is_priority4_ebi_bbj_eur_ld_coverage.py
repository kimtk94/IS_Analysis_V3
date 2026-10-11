#!/usr/bin/env python3
import csv,subprocess,json,collections,hashlib,argparse
from pathlib import Path
def run(inp,vcfs,out):
    groups=collections.defaultdict(set)
    with inp.open() as f:
        for r in csv.DictReader(f,delimiter="\t"):
            groups[r["locus"]].add(r["variant_GRCh37_REF_ALT"])
    expected={"BBJ_IS_L001":2836,"BBJ_IS_L002":2830}
    if {k:len(v) for k,v in groups.items()}!=expected:raise ValueError("SNP denominator drift")
    records=[];detail=[]
    for locus,ch in (("BBJ_IS_L001","4"),("BBJ_IS_L002","10")):
        p=vcfs/f"chr{ch}.EUR503.fullcis.b37.vcf.gz"
        result=subprocess.run(["bcftools","query","-f","%CHROM:%POS:%REF:%ALT\n",str(p)],capture_output=True,text=True,check=True)
        alleles=set()
        for line in result.stdout.splitlines():
            spl=line.split(":")
            if len(spl)!=4:continue # symbolic/structural alternate alleles are not usable SNP keys
            for alt in spl[3].split(","):alleles.add(":".join((*spl[:3],alt)))
        ok=groups[locus]&alleles
        for v in sorted(groups[locus]):
            detail.append({"locus":locus,"variant_GRCh37":v,"in_EUR503_regional_VCF":v in ok,"state":"REFERENCE_GENOTYPES_NOT_GTEX_COHORT"})
        records.append(dict(locus=locus,chromosome=ch,ABF_shared_SNPs=len(groups[locus]),
            EUR503_exact_variant_SNPs=len(ok),EUR503_missing=len(groups[locus])-len(ok),
            EUR503_exact_coverage=len(ok)/len(groups[locus]),
            eur_reference_not_qtl_donor_LD=True))
    out.mkdir(parents=True,exist_ok=True)
    for name,rr in (("IS_PRIORITY4_ABF_EUR_REFERENCE_LD_COVERAGE.tsv",records),
                    ("IS_PRIORITY4_ABF_EUR_REFERENCE_LD_SNP_KEYS.tsv",detail)):
        with (out/name).open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=rr[0].keys(),delimiter="\t")
            w.writeheader();w.writerows(rr)
    (out/"IS_PRIORITY4_ABF_EUR_REFERENCE_LD_COVERAGE.json").write_text(json.dumps({"status":"REFERENCE_PANEL_COVERAGE_ONLY","results":records,"QTL_study_matched_LD":False,"SuSiE_executed":False,"candidate_genes_unchanged":2225},indent=2))
    return records
if __name__=="__main__":
    a=argparse.ArgumentParser()
    for k in ("inp","vcfs","out"):a.add_argument("--"+k,type=Path,required=True)
    x=a.parse_args();print(json.dumps(run(x.inp,x.vcfs,x.out),indent=2))
