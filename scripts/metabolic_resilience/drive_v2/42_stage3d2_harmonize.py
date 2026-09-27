#!/usr/bin/env python3
from pathlib import Path
import csv,gzip,json,os,sys

ROOT=Path(os.environ.get("IS_ANALYSIS_ROOT","/srv/is-analysis"))
INST=ROOT/"results/metabolic_resilience/stage3_full_pgwas/instruments/STAGE3C4_FROZEN_INSTRUMENTS_ANNOTATED.tsv"
OUTCOMES=ROOT/"results/metabolic_resilience/stage3_confirmatory_mr/standardized_outcomes"
OUT=ROOT/"results/metabolic_resilience/stage3_confirmatory_mr/harmonized"
AUDIT=ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit"
OUT.mkdir(parents=True,exist_ok=True)

COMP=str.maketrans("ACGT","TGCA")
def comp_snp(a):
    return a.translate(COMP) if len(a)==1 else None
def pal(a,b):
    return len(a)==1 and len(b)==1 and {a,b} in ({"A","T"},{"C","G"})

if not INST.exists():
    print("[FAIL] C4 annotated instruments missing")
    sys.exit(2)

with INST.open("r",encoding="utf-8",newline="") as f:
    inst=list(csv.DictReader(f,delimiter="\t"))

traits={}
for p in OUTCOMES.glob("*.instrument_outcomes.tsv.gz"):
    with gzip.open(p,"rt",encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    if rows:
        traits[rows[0]["trait"]]=rows

combined=[]
audit=[]

for trait,outcomes in sorted(traits.items()):
    by_vid={}
    for o in outcomes:
        vid=o.get("matched_variant_id","")
        if vid:
            by_vid.setdefault(vid,[]).append(o)

    for x in inst:
        candidates=by_vid.get(x["variant_id"],[])
        status="NO_OUTCOME"; selected=None; aligned_beta=None; orientation=""; reason=""

        ex_ea=x["exposure_effect_allele"].upper()
        ex_oa=x["exposure_other_allele"].upper()
        try: ex_eaf=float(x["exposure_eaf"])
        except Exception: ex_eaf=None

        for o in candidates:
            oe=o["effect_allele"].upper(); oo=o["other_allele"].upper()

            direct=(oe==ex_ea and oo==ex_oa)
            swapped=(oe==ex_oa and oo==ex_ea)

            # Strand-flip matching is restricted to SNPs. Indels require exact allele representation.
            strand_direct=False; strand_swapped=False
            if len(oe)==len(oo)==len(ex_ea)==len(ex_oa)==1:
                strand_direct=(comp_snp(oe)==ex_ea and comp_snp(oo)==ex_oa)
                strand_swapped=(comp_snp(oe)==ex_oa and comp_snp(oo)==ex_ea)

            if pal(ex_ea,ex_oa):
                try:
                    out_eaf=float(o["eaf"]) if o.get("eaf") not in ("",None,"None") else None
                except Exception:
                    out_eaf=None
                if ex_eaf is None or out_eaf is None:
                    reason="PALINDROMIC_NO_EAF"; continue
                d_direct=abs(ex_eaf-out_eaf)
                d_flip=abs(ex_eaf-(1-out_eaf))
                if min(d_direct,d_flip)>0.10 or abs(d_direct-d_flip)<0.05:
                    reason="PALINDROMIC_AMBIGUOUS"; continue
                if d_direct<d_flip:
                    orientation="PAL_FREQ_DIRECT"; aligned_beta=float(o["beta"])
                else:
                    orientation="PAL_FREQ_SWAPPED"; aligned_beta=-float(o["beta"])
                selected=o; status="KEEP"; break

            if direct:
                orientation="DIRECT"; aligned_beta=float(o["beta"])
            elif swapped:
                orientation="SWAPPED"; aligned_beta=-float(o["beta"])
            elif strand_direct:
                orientation="STRAND_DIRECT"; aligned_beta=float(o["beta"])
            elif strand_swapped:
                orientation="STRAND_SWAPPED"; aligned_beta=-float(o["beta"])
            else:
                reason="ALLELE_MISMATCH"; continue

            selected=o; status="KEEP"; break

        audit.append({
            "gene":x["gene_symbol"],"mode":x["mode"],"trait":trait,
            "variant_id":x["variant_id"],"rsid":x.get("rsid",""),
            "status":status,"orientation":orientation,"reason":reason,
            "match_method":selected.get("match_method","") if selected else ""
        })

        if selected is None:
            continue

        combined.append({
            "gene":x["gene_symbol"],"protein_id":x["protein_id"],"mode":x["mode"],
            "trait":trait,"domain":selected["domain"],"outcome_type":selected["outcome_type"],
            "variant_id":x["variant_id"],"rsid":x.get("rsid",""),
            "chrom_hg19":x["chrom_hg19"],"pos_hg19":x["pos_hg19"],
            "effect_allele":ex_ea,"other_allele":ex_oa,"exposure_eaf":x["exposure_eaf"],
            "beta_exposure":float(x["exposure_beta"]),"se_exposure":float(x["exposure_se"]),
            "F":float(x["exposure_F"]),"beta_outcome":aligned_beta,
            "se_outcome":float(selected["se"]),"outcome_eaf_raw":selected.get("eaf",""),
            "outcome_n":selected.get("n",""),"orientation":orientation,
            "match_method":selected.get("match_method","")
        })

combined_path=OUT/"STAGE3D2_HARMONIZED.tsv.gz"
fields=[
    "gene","protein_id","mode","trait","domain","outcome_type","variant_id","rsid",
    "chrom_hg19","pos_hg19","effect_allele","other_allele","exposure_eaf",
    "beta_exposure","se_exposure","F","beta_outcome","se_outcome",
    "outcome_eaf_raw","outcome_n","orientation","match_method"
]
with gzip.open(combined_path,"wt",encoding="utf-8",newline="") as f:
    wr=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
    wr.writeheader(); wr.writerows(combined)

audit_path=AUDIT/"STAGE3D2_HARMONIZATION_AUDIT.tsv"
with audit_path.open("w",encoding="utf-8",newline="") as f:
    afields=list(audit[0].keys()) if audit else ["status"]
    wr=csv.DictWriter(f,fieldnames=afields,delimiter="\t",lineterminator="\n")
    wr.writeheader(); wr.writerows(audit)

trait_status={}
for trait in traits:
    a=[r for r in audit if r["trait"]==trait]
    trait_status[trait]={
        "keep":sum(r["status"]=="KEEP" for r in a),
        "excluded":sum(r["status"]!="KEEP" for r in a),
    }

summary={
    "harmonized_rows":len(combined),
    "keep":sum(r["status"]=="KEEP" for r in audit),
    "excluded":sum(r["status"]!="KEEP" for r in audit),
    "pal_freq_resolved":sum(r["orientation"].startswith("PAL_FREQ") for r in audit),
    "traits":trait_status,
    "all_traits_have_kept_variants":all(v["keep"]>0 for v in trait_status.values()) if trait_status else False,
}
(AUDIT/"STAGE3D2_HARMONIZATION_SUMMARY.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
print(json.dumps(summary,indent=2))

sys.exit(0 if summary["all_traits_have_kept_variants"] else 2)
