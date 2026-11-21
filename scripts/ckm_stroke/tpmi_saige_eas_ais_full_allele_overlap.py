#!/usr/bin/env python3
"""Separate TPMI 433.21 SAIGE vs GIGASTROKE EAS AIS GRCh37 allele intersection.

Never use TPMI PLINK2 OBS_CT, effect values or its 1,507,261 SNP result here.
Source hashes/chain checksum are pinned, source QC=PASS is mandatory, and
palindromic/ambiguous records are held out. SQLite prevents multi-million-row
RAM accumulation and de-duplicates on genomic-coordinate keys. Export only
aggregate QC; no per-SNP or individual-level GWAS data.
"""
from __future__ import annotations
import argparse
from collections import Counter
import gzip
import json
import math
from pathlib import Path
import sqlite3
import tempfile
from urllib.parse import urlsplit

from scripts.ckm_stroke.qc_tpmi_43321_saige_full_casecontrol import (
    URL as SAIGE_URL, SOURCE_BYTES as SAIGE_BYTES,
    EXPECTED_HEADER as SAIGE_HEADER,
)
from scripts.ckm_stroke.qc_gigastroke_eas_ais_original import (
    URL as AIS_URL, SCHEMA as AIS_HEADER,
)
from scripts.ckm_stroke.tpmi_gigastroke_eas_genomewide_overlap import (
    download_pinned, download_chain, CHAIN_URL, CHAIN_MD5,
    AIS_SHA256, AIS_ROWS, unique_lift, translate_snv, allele_align,
)

SAIGE_SHA256="9f235bfec6265442b7699accecbbd2a69c5acb369f7ee1bcc9e79ec94285a33e"
SAIGE_ROWS=15_588_353
SAIGE_QC_PASS=8_042_797
SAIGE_CASES=9249
SAIGE_CONTROLS=304660
DNA=set("ACGT")
AUTOSOMES={str(i) for i in range(1,23)}
PALINDROMES=({"A","T"},{"C","G"})
CHUNK_SIZE=12500

def saige_candidate(row):
    """One source SAIGE row -> CHR38 POS38 REF ALT EFFECT_AF and effect stats.

    MarkerID stores ref/alt allele strings; order of Allele1/Allele2 is not
    assumed to equal ref/alt. BETA and AF_Allele2 are for Allele2.
    Return reason and tuple. Fail closed on any source-QC other than PASS.
    """
    if len(row)!=len(SAIGE_HEADER):
        raise ValueError("SAIGE original header-width drift")
    if row[19]!="PASS":
        return "SOURCE_QC_"+row[19][:24].upper(),None
    chrom,pos,marker,allele1,allele2=row[:5]
    if chrom not in AUTOSOMES or not pos.isdecimal() or int(pos)<=0:
        return "NONAUTOSOMAL_OR_BAD_COORDINATE",None
    pieces=marker.removeprefix("chr").replace("_",":").split(":")
    if len(pieces)!=4 or pieces[0]!=chrom or pieces[1]!=pos:
        return "MALFORMED_MARKER_ID",None
    ref,alt=pieces[2].upper(),pieces[3].upper()
    a1,a2=allele1.upper(),allele2.upper()
    if not all(len(a)==1 and a in DNA for a in (ref,alt,a1,a2)) or ref==alt:
        return "NOT_BIALLELIC_SNV",None
    if {ref,alt}!={a1,a2}:return "MARKER_EFFECT_ALLELE_MISMATCH",None
    if {ref,alt} in PALINDROMES:return "PALINDROMIC_UNRESOLVED",None
    try:
        af=float(row[6]);info=float(row[7]);beta=float(row[8])
        se=float(row[9]);pval=float(row[12])
        ncase=int(row[17]);nctrl=int(row[18])
    except (TypeError,ValueError,OverflowError):
        return "BAD_NUMERIC",None
    if (not all(math.isfinite(v) for v in (af,info,beta,se,pval)) or
            not 0<=af<=1 or not 0<=info<=1.001 or se<=0 or
            not 0<=pval<=1 or
            (ncase,nctrl)!=(SAIGE_CASES,SAIGE_CONTROLS)):
        return "BAD_NUMERIC_OR_N",None
    return "CANDIDATE",(chrom,int(pos),ref,alt,a2,af,info,beta,se,pval)

def make_sqlite(path):
    conn=sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute("PRAGMA cache_size=-180000")
    conn.execute("PRAGMA locking_mode=EXCLUSIVE")
    conn.execute("""CREATE TABLE outcome(
      chrom TEXT NOT NULL,pos INTEGER NOT NULL,ea TEXT,oa TEXT,eaf REAL,
      beta REAL,se REAL,pval REAL,ambiguous INTEGER NOT NULL DEFAULT 0,
      PRIMARY KEY(chrom,pos)) WITHOUT ROWID""")
    conn.execute("""CREATE TABLE saige(
      chrom TEXT NOT NULL,pos INTEGER NOT NULL,r TEXT,a TEXT,eff TEXT,
      freq REAL,info REAL,beta REAL,se REAL,pval REAL,
      ambiguous INTEGER NOT NULL DEFAULT 0,
      PRIMARY KEY(chrom,pos)) WITHOUT ROWID""")
    return conn

def stage_outcome(conn,path):
    c=Counter(); batch=[]
    sql="""INSERT INTO outcome(chrom,pos,ea,oa,eaf,beta,se,pval)
             VALUES(?,?,?,?,?,?,?,?)
             ON CONFLICT(chrom,pos) DO UPDATE SET ambiguous=1"""
    with gzip.open(path,"rt",encoding="utf-8-sig") as f:
        if tuple(f.readline().rstrip("\r\n").split("\t"))!=AIS_HEADER:
            raise ValueError("GIGASTROKE original EAS AIS header drift")
        for line in f:
            c["all_original_rows"]+=1
            row=line.rstrip("\r\n").split("\t")
            if len(row)!=len(AIS_HEADER):raise ValueError("EAS original row schema drift")
            if row[0] not in AUTOSOMES or not row[1].isdecimal():
                c["invalid_chr_or_pos"]+=1;continue
            ea,oa=row[9].upper(),row[10].upper()
            if (len(ea)!=1 or len(oa)!=1 or
                ea not in DNA or oa not in DNA or ea==oa):
                c["non_biallelic_SNV"]+=1;continue
            try:
                af,beta,se,p=float(row[2]),float(row[3]),float(row[4]),float(row[5])
            except ValueError:
                c["invalid_numeric"]+=1;continue
            if not all(map(math.isfinite,(af,beta,se,p))) or not 0<=af<=1 or se<=0 or not 0<=p<=1:
                c["invalid_numeric"]+=1;continue
            batch.append((row[0],int(row[1]),ea,oa,af,beta,se,p))
            c["candidate_outcome_rows"]+=1
            if len(batch)>=CHUNK_SIZE:
                conn.executemany(sql,batch);conn.commit();batch=[]
        if batch:conn.executemany(sql,batch);conn.commit()
    if c["all_original_rows"]!=AIS_ROWS:
        raise ValueError("GIGASTROKE source rows no longer pinned")
    c["unique_outcome_coordinates"]=conn.execute("SELECT count(*) FROM outcome").fetchone()[0]
    c["outcome_duplicate_coordinate_positions"]=conn.execute(
       "SELECT count(*) FROM outcome WHERE ambiguous=1").fetchone()[0]
    return dict(c)

def stage_saige(conn,path,lifter):
    c=Counter(); batch=[]
    sql="""INSERT INTO saige(chrom,pos,r,a,eff,freq,info,beta,se,pval)
             VALUES(?,?,?,?,?,?,?,?,?,?)
             ON CONFLICT(chrom,pos) DO UPDATE SET ambiguous=1"""
    with gzip.open(path,"rt",encoding="utf-8-sig") as f:
        if tuple(f.readline().rstrip("\r\n").split("\t"))!=SAIGE_HEADER:
            raise ValueError("TPMI original SAIGE header drift")
        for line in f:
            c["all_original_rows"]+=1
            row=line.rstrip("\r\n").split("\t")
            status,x=saige_candidate(row)
            c[status]+=1
            if status!="CANDIDATE":continue
            chrom,p38,ref,alt,eff,af,info,beta,se,p=x
            outcome,coord=unique_lift(lifter,chrom,p38)
            if outcome!="UNIQUE":
                c["LIFT_"+outcome]+=1;continue
            c37,p37,strand=coord
            r,a,eff=translate_snv(ref,alt,eff,strand)
            batch.append((c37,p37,r,a,eff,af,info,beta,se,p))
            if len(batch)>=CHUNK_SIZE:
                conn.executemany(sql,batch);conn.commit();batch=[]
            c["LIFT_UNIQUE"]+=1
        if batch:conn.executemany(sql,batch);conn.commit()
    if c["all_original_rows"]!=SAIGE_ROWS or c["CANDIDATE"]>SAIGE_QC_PASS:
        raise ValueError("TPMI SAIGE source rows/QC source count drift")
    if c["SOURCE_QC_PASS"]:
        raise ValueError("source QC=PASS should enter CANDIDATE rather than filtered count")
    if c["all_original_rows"]-c["SOURCE_QC_FAIL"]!=SAIGE_QC_PASS:
        raise ValueError("TPMI SAIGE pinned QC=PASS count drift")
    c["unique_saige_lifted_coords"]=conn.execute("SELECT count(*) FROM saige").fetchone()[0]
    c["saige_ambiguous_lifted_coords"]=conn.execute(
      "SELECT count(*) FROM saige WHERE ambiguous=1").fetchone()[0]
    return dict(c)

def compare(conn):
    c=Counter();freq=Counter();chr_match=Counter()
    q="""SELECT s.chrom,s.r,s.a,s.eff,s.freq,s.info,s.beta,s.se,s.pval,
               o.ea,o.oa,o.eaf,o.beta,o.se,o.pval
        FROM saige AS s INNER JOIN outcome AS o
        ON s.chrom=o.chrom AND s.pos=o.pos
        WHERE s.ambiguous=0 AND o.ambiguous=0"""
    for v in conn.execute(q):
        chrom,r,a,eff,af,info,beta,se,pval,ea,oa,eaf,ob,ose,op=v
        c["COORDINATE_INTERSECTIONS"]+=1
        status,sign,delta=allele_align(r,a,eff,af,ea,oa,eaf)
        c[status]+=1
        if status!="MATCHED":continue
        chr_match[chrom]+=1
        if delta<=0.05:freq["DELTA_LE_0.05"]+=1
        elif delta<=0.1:freq["DELTA_GT_0.05_LE_0.1"]+=1
        elif delta<=0.2:freq["DELTA_GT_0.1_LE_0.2"]+=1
        else:freq["DELTA_GT_0.2"]+=1
        if info>=0.8:c["INFO_GE_0.8"]+=1
        if min(af,1-af)>=0.01:c["SAIGE_MAF_GE_0.01"]+=1
        if info>=0.8 and min(af,1-af)>=0.01:
            c["INFO_GE_0.8_AND_SAIGE_MAF_GE_0.01"]+=1
            if delta<=0.05:
                c["INFO_GE_0.8_MAF_GE_0.01_AF_DELTA_LE_0.05"]+=1
            if delta<=0.10:
                c["INFO_GE_0.8_MAF_GE_0.01_AF_DELTA_LE_0.10"]+=1
        if delta>0.05:c["AF_DELTA_GT_0.05"]+=1
        if delta>0.1:c["AF_DELTA_GT_0.1"]+=1
    return dict(c),dict(freq),dict(chr_match)

def execute():
    from pyliftover import LiftOver
    with tempfile.TemporaryDirectory(prefix="is_saige_eas_overlap_") as d:
        root=Path(d)
        saige=root/"saige_source.gz"
        eas=root/"eas_source.gz"
        chain_zip=root/"chain.gz"
        chain_text=root/"chain.over"
        src=download_pinned(SAIGE_URL,urlsplit(SAIGE_URL).hostname,
            saige,SAIGE_SHA256,900000000,SAIGE_BYTES)
        outcome=download_pinned(AIS_URL,"ftp.ebi.ac.uk",eas,AIS_SHA256,1500000000)
        ucsc=download_chain(chain_zip)
        with gzip.open(chain_zip,"rb") as fi,chain_text.open("wb") as fo:
            while buf:=fi.read(1048576):fo.write(buf)
        lifter=LiftOver(str(chain_text))
        # Independent original TPMI PLINK2 known locus validates the chain
        # direction, but NOT SAIGE locus-specific effects or variant identity.
        code,x=unique_lift(lifter,"19",15462568)
        if code!="UNIQUE" or x!=("19",15573379,"+"):
            raise ValueError("GRCh38->37 UCSC chain control failed")
        con=make_sqlite(root/"join.sqlite")
        eas_qc=stage_outcome(con,eas)
        saige_qc=stage_saige(con,saige,lifter)
        results,af_bins,bychr=compare(con)
        con.close()
        status=("PASS_TECHNICAL_SAIGE_EAS_ALLELE_INTERSECTION" if
                results.get("MATCHED",0)>0 else "HOLD_NO_ALLELE_MATCHES")
        return {
           "schema":"TPMI-SAIGE-43321-vs-GIGASTROKE-EAS-AIS-allele-overlap-aggregate-v1",
           "technical_status":status,
           "TPMI_source_model":"SAIGE",
           "NOT_TPMI_PLINK2":True,
           "SAIGE_original_sha256_pinned":src["sha256"],
           "SAIGE_original_rows":SAIGE_ROWS,
           "SAIGE_full_original_QC_PASS":SAIGE_QC_PASS,
           "SAIGE_source_GRCh38_official":"https://pheweb.ibms.sinica.edu.tw/about",
           "SAIGE_effect_allele":"Allele2",
           "SAIGE_effect_allele_source":"https://saigegit.github.io/SAIGE-doc/docs/single_step2.html",
           "GIGASTROKE_EAS_AIS_source_sha256_pinned":outcome["sha256"],
           "GIGASTROKE_EAS_AIS_source_build":"GRCh37",
           "UCSC_hg38ToHg19_chain_md5":ucsc["md5"],
           "SAIGE_source_filter_liftover_QC":saige_qc,
           "GIGASTROKE_EAS_source_index_QC":eas_qc,
           "allele_overlap_QC":results,
           "AF_delta_bins":af_bins,
           "matched_by_chromosome":bychr,
           "effective_N_SAIGE_case":SAIGE_CASES,
           "effective_N_SAIGE_control":SAIGE_CONTROLS,
           "unadjusted_binary_harmonic_N_nominal_not_LDSC_Neff":
               4*SAIGE_CASES*SAIGE_CONTROLS/(SAIGE_CASES+SAIGE_CONTROLS),
           "unadjusted_harmonic_N_not_approved_for_LDSC":True,
           "effective_N_is_model_specific":True,
           "LDSC_MR_coloc_performed":False,
           "EAS_sample_overlap_certified":False,
           "biological_effect_concordance_tested":False,
           "individual_SNP_effects_exported":False,
           "scientific_status":"HOLD",
        }

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    try:res=execute()
    except Exception as e:
        res={"schema":"TPMI-SAIGE-43321-vs-GIGASTROKE-EAS-AIS-allele-overlap-aggregate-v1",
             "technical_status":"HOLD_SOURCE_BUILD_ALLELE_OR_LIFTOVER",
             "error_type":type(e).__name__,
             "error":str(e)[:340],"scientific_status":"HOLD"}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(res,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "technical_status":res["technical_status"],
        "SAIGE_total_rows":res.get("SAIGE_original_rows"),
        "SAIGE_source_filter_liftover_QC":res.get("SAIGE_source_filter_liftover_QC"),
        "GIGASTROKE_EAS_source_index_QC":res.get("GIGASTROKE_EAS_source_index_QC"),
        "allele_overlap_QC":res.get("allele_overlap_QC"),
        "AF_delta_bins":res.get("AF_delta_bins"),
        "error":res.get("error"),
        "scientific_status":res["scientific_status"]},ensure_ascii=False))
    if not res["technical_status"].startswith("PASS"):raise SystemExit(1)
if __name__=="__main__":main()
