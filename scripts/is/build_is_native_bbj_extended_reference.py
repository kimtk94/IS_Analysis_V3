#!/usr/bin/env python3
"""Freeze original BBJ GWAS and chain-lifted four-locus reference WITHOUT QTL input.

The reference is independent from GTEx source input. It is not an independent
GRCh38 FASTA base proof, molecular cohort LD, or causal inference result.
"""
import csv,gzip,hashlib,json,sys,zipfile,collections,io
from decimal import Decimal, InvalidOperation
from pathlib import Path

WINDOWS={
 "BBJ_IS_L001":("4",79681087,82684150),
 "BBJ_IS_L002":("10",103839152,106959690),
 "BBJ_IS_L003":("12",110849271,114615871),
 "BBJ_IS_L004":("13",109454304,112459229),
}
ZIP=Path("/srv/is-analysis/data/is/east_asia/japan/bbj/hum0197.v3.BBJ.IS.v1.zip")
CHAIN=Path("/srv/is-analysis/data/is/reference/hg19ToHg38.over.chain.gz")
VENDOR=Path("/srv/is-analysis/data/is/reference/vendor_pyliftover")
OUT=Path("/srv/is-analysis/IS_G2_SDY_20261011/data/is_reference/native_bbj_is_4loci_extended_20261011_v2.fingerprints.tsv.gz")
REPORT=Path("/srv/is-analysis/IS_G2_SDY_20261011/data/is_reference/native_bbj_is_4loci_extended_20261011_v2.manifest.json")
if OUT.exists() or REPORT.exists():raise FileExistsError("NO_REFERENCE_OVERWRITE")
sys.path.insert(0,str(VENDOR))
from pyliftover import LiftOver
lo=LiftOver(str(CHAIN))
def revcomp(s):
    return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
total=0
hits=[]
dup=collections.Counter()
reject=collections.Counter()
with zipfile.ZipFile(ZIP) as z:
    names=[x for x in z.namelist() if x.endswith(".auto.txt.gz")]
    if len(names)!=1:raise ValueError("SINGLE_AUTOSOME_MEMBER_EXPECTED")
    with z.open(names[0]) as zf:
        with gzip.open(zf,"rt") as f:
            reader=csv.DictReader(f,delimiter="\t")
            required=("v","CHR","POS","Allele1","Allele2","BETA","AF_Allele2")
            if not set(required).issubset(reader.fieldnames):raise ValueError("NATIVE_SCHEMA")
            for r in reader:
                total+=1
                for locus,(ch,start,stop) in WINDOWS.items():
                    if r["CHR"]!=ch:continue
                    pos_numeric=Decimal(r["POS"])
                    if not pos_numeric.is_finite() or pos_numeric!=int(pos_numeric):
                        reject[(locus,"NONINTEGER_NATIVE_POSITION")]+=1
                        break
                    pos=int(pos_numeric)
                    if pos<start or pos>stop:break
                    try:
                        k=r["v"]
                        parts=k.split(":")
                        if len(parts)!=4:raise ValueError("ID_PARTS")
                        ch1,p1,ref,alt=parts
                        if ch1!=ch or int(p1)!=pos:raise ValueError("CHR_POS")
                        ref=ref.upper();alt=alt.upper()
                        if not ref or not alt or any(x not in "ACGT" for x in ref+alt):raise ValueError("NON_SNV")
                        a1=r["Allele1"].upper();a2=r["Allele2"].upper()
                        if {a1,a2}!={ref,alt}:raise ValueError("NATIVE_ALLELES_NOT_ID")
                        if a2!=alt:raise ValueError("NATIVE_EFFECT_ALLELE_NOT_ID_ALT")
                        mappings=lo.convert_coordinate("chr"+ch,pos-1)
                        if not mappings or len(mappings)!=1:raise ValueError("NOT_UNIQUE_CHAIN")
                        c,p,strand,_=mappings[0]
                        if strand not in ("+","-"):raise ValueError("CHAIN_STRAND")
                        er,ea=(ref,alt) if strand=="+" else (revcomp(ref),revcomp(alt))
                        mk=f"{c.removeprefix('chr')}:{int(p)+1}:{er}:{ea}"
                        hits.append((locus,k,mk,r["BETA"],r["AF_Allele2"],strand))
                    except Exception as exc:
                        reject[(locus,str(exc))]+=1
                    break
def source_fingerprint(locus,variant_id,match_key,beta,eaf):
    normalized="IS_BBJ_NATIVE_4LOCI_V1\t"+"\t".join((
        locus,variant_id,match_key,format(float(beta),".12f"),format(float(eaf),".12f")))
    return hashlib.sha256(normalized.encode()).hexdigest()[:40]
fields=["locus","variant_id","native_variant_fingerprint_sha256_160bits"]
hits.sort(key=lambda t:(t[0],t[1]))
if len(hits)!=len(set((x[0],x[1]) for x in hits)):
    raise ValueError("ORIGINAL_NATIVE_DUPLICATE_CHR37_SNP_IN_WINDOW")
OUT.parent.mkdir(parents=True,exist_ok=True)
with OUT.open("wb") as raw, gzip.GzipFile(fileobj=raw,filename="",mode="wb",mtime=0,compresslevel=9) as gz, io.TextIOWrapper(gz,encoding="utf8",newline="") as f:
    w=csv.writer(f,delimiter="\t",lineterminator="\n")
    w.writerow(fields)
    for locus,variant,mk,beta,eaf,strand in hits:
        w.writerow([locus,variant,source_fingerprint(locus,variant,mk,beta,eaf)])
info=dict(schema="IS_NATIVE_BBJ_ORIGINAL_ZIP_LIFTOVER_FOUR_LOCI_PADDED_V1",
    original_native_GWAS_stream_rows=total,reference_variants=len(hits),
    locus_counts=dict(collections.Counter(x[0] for x in hits)),
    native_excluded_records=dict((f"{k[0]}::{k[1]}",v) for k,v in reject.items()),
    GRCh37_windows=WINDOWS,original_native_zip_size=ZIP.stat().st_size,
    native_zip_SHA256=hashlib.sha256(ZIP.read_bytes()).hexdigest(),
    chain_SHA256=hashlib.sha256(CHAIN.read_bytes()).hexdigest(),
    gz_sha256=hashlib.sha256(OUT.read_bytes()).hexdigest(),
    required_effect_allele="BBJ_NATIVE_ALLELE2_EQUALS_GRCH37_VARIANT_ID_ALT",
    interpretation="Source original native GWAS effect allele + chain verified fingerprints only (no redistributed native GWAS beta/EAF). Not GRCh38 FASTA, not matched cohort LD")
REPORT.write_text(json.dumps(info,ensure_ascii=False,indent=2)+"\n")
print(json.dumps(info,ensure_ascii=False))
