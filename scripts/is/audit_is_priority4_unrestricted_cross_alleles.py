#!/usr/bin/env python3
"""Cross match four unrestricted leads by position/allele-pair, never raw ID."""
import gzip,csv,json,sys
from pathlib import Path
SITES={"4":{81182554,81181072},"10":{105459834,105599770}}
FILES={"BBJ":"/srv/is-analysis/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz",
"GIGA":"/srv/is-analysis/data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz"}
found={}
for name,fp in FILES.items():
 with gzip.open(fp,"rt") as f:
  for r in csv.DictReader(f,delimiter="\t"):
   c=r["chr"]
   if c in SITES and int(r["pos"]) in SITES[c]:
    k=c+":"+r["pos"]
    if (name,k) in found:raise ValueError("Duplicated coordinate")
    found[name,k]=r
  print("SCANNED",name,flush=True)
rows=[]
for c,poss in SITES.items():
 for pos in sorted(poss):
  k=c+":"+str(pos)
  a=found.get(("BBJ",k));b=found.get(("GIGA",k))
  v=dict(site=k,BBJ_variant=a["variant_id"] if a else "",GIGA_variant=b["variant_id"] if b else "",
   BBJ_p=a["p"] if a else "",GIGA_p=b["p"] if b else "",
   BBJ_effect=a["effect_allele"] if a else "",GIGA_effect=b["effect_allele"] if b else "",
   BBJ_beta=a["beta"] if a else "",GIGA_beta=b["beta"] if b else "",
   pair_match="",palindromic="",aligned_beta_sign_agreement="")
  if a and b:
   bp=set((a["ref"],a["alt"]));gp=set(b["variant_id"].split(":")[-2:])
   v["pair_match"]=bp==gp
   v["palindromic"]=bp in ({"A","T"},{"C","G"})
   if bp==gp and a["effect_allele"] in bp and b["effect_allele"] in bp:
    same=(float(a["beta"])*(1 if a["effect_allele"]==a["alt"] else -1))*(float(b["beta"])*(1 if b["effect_allele"]==a["alt"] else -1))
    v["aligned_beta_sign_agreement"]=same>0 if not v["palindromic"] else "PALINDROME_STRAND_UNRESOLVED"
  rows.append(v)
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
with (out/"LEAD_POS_ALLELE_CROSS.tsv").open("w") as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
(out/"LEAD_POS_ALLELE_CROSS.json").write_text(json.dumps(rows,indent=2)+"\n")
print(json.dumps(rows,indent=2))
