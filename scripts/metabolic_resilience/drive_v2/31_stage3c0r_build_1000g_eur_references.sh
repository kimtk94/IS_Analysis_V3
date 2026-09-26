#!/usr/bin/env bash
ROOT="${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
cd "$ROOT" || exit 1
set +e
set +u
set +o pipefail 2>/dev/null || true

LDROOT="$ROOT/data/metabolic_resilience/stage2_gwas/ld_reference_1kg_eur"
WINDOWS="$ROOT/results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2_GENE_CIS_WINDOWS.tsv"
PRECLUMP="$ROOT/results/metabolic_resilience/stage3_full_pgwas/ld_clump"
EUR="$LDROOT/EUR503.samples"
VCFDIR="$LDROOT/stage3_regional_vcf"
PGENDIR="$LDROOT/stage3_regional_pgen"
CACHE="$LDROOT/stage3_chromosome_cache"
AUDIT="$ROOT/results/metabolic_resilience/stage3_full_pgwas/audit"
REGIONS="$AUDIT/STAGE3C0R_REGIONS.tsv"
REFSUMMARY="$AUDIT/STAGE3C0R_REFERENCE_BUILD.tsv"
MAPSUMMARY="$AUDIT/STAGE3C0R_LD_MAPPING_SUMMARY.tsv"
MAPJSON="$AUDIT/STAGE3C0R_LD_MAPPING_SUMMARY.json"
mkdir -p "$VCFDIR" "$PGENDIR" "$CACHE" "$AUDIT"

BCFTOOLS="$(command -v bcftools)"
PLINK2="$ROOT/tools/plink2-a6.39-20260919/plink2"
[ -x "$PLINK2" ] || PLINK2="$(command -v plink2)"

echo "===================================================="
echo "STAGE 3-C0R : REBUILD 8 CANDIDATE 1000G EUR REFERENCES"
echo "===================================================="

READY=1
for F in "$WINDOWS" "$EUR"; do
  if [ -s "$F" ]; then echo "[PASS] $F"; else echo "[FAIL] $F"; READY=0; fi
done
if [ -x "$BCFTOOLS" ]; then echo "[PASS] bcftools=$BCFTOOLS"; else echo "[FAIL] bcftools"; READY=0; fi
if [ -x "$PLINK2" ]; then echo "[PASS] plink2=$PLINK2"; else echo "[FAIL] plink2"; READY=0; fi
echo "EUR sample count: $(grep -cve '^[[:space:]]*$' "$EUR" 2>/dev/null)"

# clean stale outputs from failed v5a run
find "$CACHE" -maxdepth 1 -type f -name '*integrated_v5a*' -print -delete 2>/dev/null
rm -f "$VCFDIR"/*.vcf.gz "$VCFDIR"/*.vcf.gz.tbi "$PGENDIR"/*.pgen "$PGENDIR"/*.pvar "$PGENDIR"/*.psam "$PGENDIR"/*.log 2>/dev/null

python3 - <<'PY'
from pathlib import Path
import csv
ROOT=Path('/srv/is-analysis')
WINDOWS=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2_GENE_CIS_WINDOWS.tsv'
OUT=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C0R_REGIONS.tsv'
wanted={'AOC1','OGN','TNFRSF6B','TFPI','SULT1A1','CTRL','IDUA','COMT'}
rows=[]
with WINDOWS.open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['gene_symbol'] in wanted and r['assembly']=='GRCh37':
            rows.append({'protein_id':r['protein_id'],'gene':r['gene_symbol'],'chrom':r['chrom'],'start':r['cis_start'],'end':r['cis_end']})
rows.sort(key=lambda x:(int(x['chrom']),int(x['start'])))
with OUT.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print('targets=',len(rows))
for r in rows: print(r['gene'],f"chr{r['chrom']}:{r['start']}-{r['end']}")
PY
REGION_RC=$?

BASEURL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
CHROMS="$(awk -F'\t' 'NR>1 {print $3}' "$REGIONS" | sort -n | uniq)"
URL_READY=1
for CHR in $CHROMS; do
  BASE="ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  VCF_HTTP="$(curl -L -s -o /dev/null -w '%{http_code}' -I "$BASEURL/$BASE")"
  TBI_HTTP="$(curl -L -s -o /dev/null -w '%{http_code}' -I "$BASEURL/$BASE.tbi")"
  echo "chr$CHR VCF=$VCF_HTTP TBI=$TBI_HTTP"
  if [ "$VCF_HTTP" != "200" ] || [ "$TBI_HTTP" != "200" ]; then URL_READY=0; fi
done
[ "$URL_READY" = "1" ] && echo "[PASS] all required v5b URLs available" || echo "[FAIL] one or more v5b URLs unavailable"

printf 'gene\tchrom\tstart\tend\tvcf_variants\tpgen_variants\tpgen_samples\tstatus\n' > "$REFSUMMARY"

if [ "$READY" = "1" ] && [ "$REGION_RC" = "0" ] && [ "$URL_READY" = "1" ]; then
for CHR in $CHROMS; do
  BASE="ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  URL="$BASEURL/$BASE"; VCF="$CACHE/$BASE"; TBI="$VCF.tbi"
  if [ ! -s "$VCF" ]; then rm -f "$VCF"; wget -c -O "$VCF" "$URL"; VCF_RC=$?; else VCF_RC=0; fi
  if [ ! -s "$TBI" ]; then rm -f "$TBI"; wget -c -O "$TBI" "$URL.tbi"; TBI_RC=$?; else TBI_RC=0; fi
  echo "chr$CHR vcf_rc=$VCF_RC tbi_rc=$TBI_RC"
  if [ "$VCF_RC" != "0" ] || [ "$TBI_RC" != "0" ] || [ ! -s "$VCF" ] || [ ! -s "$TBI" ]; then
    echo "[FAIL] chr$CHR source download"; rm -f "$VCF" "$TBI"; continue
  fi
  "$BCFTOOLS" index -n "$VCF" >/dev/null 2>&1; INDEX_RC=$?
  if [ "$INDEX_RC" != "0" ]; then echo "[FAIL] chr$CHR source/index"; rm -f "$VCF" "$TBI"; continue; fi

  while IFS=$'\t' read -r PID GENE GCHR START END; do
    [ "$GCHR" != "$CHR" ] && continue
    OUTVCF="$VCFDIR/$GENE.1000G_EUR.b37.vcf.gz"; TMP="$VCFDIR/$GENE.tmp.vcf.gz"; PREFIX="$PGENDIR/$GENE.1000G_EUR.b37"
    rm -f "$TMP" "$TMP.tbi" "$OUTVCF" "$OUTVCF.tbi"
    "$BCFTOOLS" view --threads 2 -r "$CHR:$START-$END" -S "$EUR" -Oz -o "$TMP" "$VCF"; VIEW_RC=$?
    if [ "$VIEW_RC" != "0" ] || [ ! -s "$TMP" ]; then
      printf '%s\t%s\t%s\t%s\t0\t0\t0\tFAIL_VCF\n' "$GENE" "$CHR" "$START" "$END" >> "$REFSUMMARY"; continue
    fi
    "$BCFTOOLS" index -t -f "$TMP"; IDX_RC=$?
    if [ "$IDX_RC" != "0" ]; then
      printf '%s\t%s\t%s\t%s\t0\t0\t0\tFAIL_INDEX\n' "$GENE" "$CHR" "$START" "$END" >> "$REFSUMMARY"; rm -f "$TMP" "$TMP.tbi"; continue
    fi
    mv "$TMP" "$OUTVCF"; mv "$TMP.tbi" "$OUTVCF.tbi"
    N_VCF="$("$BCFTOOLS" view -H "$OUTVCF" 2>/dev/null | wc -l)"
    rm -f "$PREFIX.pgen" "$PREFIX.pvar" "$PREFIX.psam" "$PREFIX.log"
    "$PLINK2" --vcf "$OUTVCF" --double-id --maf 0.01 --min-alleles 2 --max-alleles 2 --set-all-var-ids '@:#:$r:$a' --make-pgen --threads 2 --memory 2500 --out "$PREFIX"
    PLINK_RC=$?
    if [ "$PLINK_RC" = "0" ] && [ -s "$PREFIX.pgen" ] && [ -s "$PREFIX.pvar" ] && [ -s "$PREFIX.psam" ]; then
      N_PVAR="$(grep -vc '^#' "$PREFIX.pvar")"; N_PSAM="$(grep -vc '^#' "$PREFIX.psam")"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' "$GENE" "$CHR" "$START" "$END" "$N_VCF" "$N_PVAR" "$N_PSAM" >> "$REFSUMMARY"
    else
      printf '%s\t%s\t%s\t%s\t%s\t0\t0\tFAIL_PGEN\n' "$GENE" "$CHR" "$START" "$END" "$N_VCF" >> "$REFSUMMARY"
    fi
  done < <(tail -n +2 "$REGIONS")

  rm -f "$VCF" "$TBI"
done
fi

column -t -s $'\t' "$REFSUMMARY"
REF_PASS="$(awk -F'\t' 'NR>1 && $NF=="PASS" {n++} END {print n+0}' "$REFSUMMARY")"
echo "reference_pass=$REF_PASS/8"

python3 - <<'PY'
from pathlib import Path
import csv,json
ROOT=Path('/srv/is-analysis')
PRE=ROOT/'results/metabolic_resilience/stage3_full_pgwas/ld_clump'
REF=ROOT/'data/metabolic_resilience/stage2_gwas/ld_reference_1kg_eur/stage3_regional_pgen'
AUDIT=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit'
OUT=AUDIT/'STAGE3C0R_LD_MAPPING_SUMMARY.tsv'; JSONOUT=AUDIT/'STAGE3C0R_LD_MAPPING_SUMMARY.json'
genes=['AOC1','OGN','TNFRSF6B','TFPI','SULT1A1','CTRL','IDUA','COMT']; rows=[]
for gene in genes:
    cf=PRE/f'{gene}_primary_preclump.tsv'; pvar=REF/f'{gene}.1000G_EUR.b37.pvar'
    if not cf.exists() or not pvar.exists():
        rows.append({'gene':gene,'candidate_n':0,'position_match':0,'allele_match':0,'position_pct':0,'allele_pct':0,'status':'FAIL_MISSING_INPUT'});continue
    candidates={}
    with cf.open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            key=(int(r['chrom_hg19']),int(r['pos_hg19']))
            candidates.setdefault(key,[]).append({'variant_id':r['variant_id'],'a0':r['allele0'].upper(),'a1':r['allele1'].upper()})
    pos_ids=set(); allele_ids=set(); header=None
    with pvar.open(errors='replace') as f:
        for line in f:
            if line.startswith('##'): continue
            if header is None:
                header=line.lstrip('#').split(); H={x:i for i,x in enumerate(header)}; continue
            p=line.split()
            if len(p)<len(header): continue
            try: key=(int(p[H['CHROM']]),int(p[H['POS']]))
            except: continue
            if key not in candidates: continue
            ref=p[H['REF']].upper(); alt=p[H['ALT']].upper()
            for c in candidates[key]:
                pos_ids.add(c['variant_id'])
                if {ref,alt}=={c['a0'],c['a1']}: allele_ids.add(c['variant_id'])
    all_ids={c['variant_id'] for arr in candidates.values() for c in arr}; n=len(all_ids); np=len(pos_ids); na=len(allele_ids)
    pp=100*np/n if n else 0; ap=100*na/n if n else 0
    status='PASS' if ap>=90 else ('WARN' if ap>=70 else 'FAIL_LOW_COVERAGE')
    rows.append({'gene':gene,'candidate_n':n,'position_match':np,'allele_match':na,'position_pct':round(pp,3),'allele_pct':round(ap,3),'status':status})
with OUT.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
total=sum(r['candidate_n'] for r in rows); matched=sum(r['allele_match'] for r in rows)
s={'targets':len(rows),'pass_ge90pct':sum(r['status']=='PASS' for r in rows),'warn_70_90pct':sum(r['status']=='WARN' for r in rows),
   'fail_lt70pct':sum(r['status'].startswith('FAIL') for r in rows),'primary_candidates_total':total,'allele_matched_total':matched,
   'overall_allele_match_pct':round(100*matched/total,3) if total else 0,'reference':'1000 Genomes Phase 3 EUR503 GRCh37 v5b',
   'reference_maf':0.01,'reference_variant_type':'biallelic SNPs and indels','next_step_if_pass':'PLINK2 LD clumping r2<0.01 using candidate-specific reference'}
JSONOUT.write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2))
PY
MAP_RC=$?

echo "mapping_rc=$MAP_RC"
column -t -s $'\t' "$MAPSUMMARY" 2>/dev/null
cat "$MAPJSON" 2>/dev/null
if [ "$REF_PASS" = "8" ] && [ "$MAP_RC" = "0" ]; then
  echo "[PASS] STAGE 3-C0R REFERENCE BUILD COMPLETE"
else
  echo "[PARTIAL] Reference build requires review."
  echo "Shell remains active."
fi
