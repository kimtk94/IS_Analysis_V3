#!/usr/bin/env python3
"""Read-only original ImmuNexUT E-GEAD-420 nominal all-association source.

Verifies publisher filelist MD5 for small IDF metadata and probes only
ZIP64 central directory and <=2MiB first original data bytes of the huge
43.5GB public processed archive. ZIP publisher MD5 CANNOT be verified without
full source download: never claim it has passed. No participant values, donor
IDs, genotypes, or QTL associations are exported.

Study: Japanese immune cells n=416, 28 subsets, mixed 10 immune diseases and
healthy donors. Full nominal cis-pairs including nonsignificant tests are
described in official NBDC documentation, NOT independently row-complete yet.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import re
import struct
from urllib.request import Request,urlopen
from urllib.parse import urlparse
import urllib.error
import zlib

BASE="https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-420/"
FILELIST=BASE+"E-GEAD-420.filelist.txt"
IDF=BASE+"E-GEAD-420.idf.txt"
ZIP=BASE+"E-GEAD-420.processed.zip"
MAX_TAIL=262144
MAX_CD=2*1024*1024
MAX_PREFIX=2*1024*1024
EXPECTED_ZIP_SIZE=43507433191
EXPECTED_PUBLISHED_ZIP_MD5="0231fc7ca80b4a604f31661d72ae62d8"


def check_url(url):
    p=urlparse(url)
    if p.scheme!="https" or p.hostname!="ddbj.nig.ac.jp":
        raise ValueError("Unapproved host")
    if not p.path.startswith("/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-420/"):
        raise ValueError("Source is not the exact official dataset")


def get_small(url,limit=512000):
    check_url(url)
    req=Request(url,headers={"User-Agent":"IS_Analysis_V3_ImmuNexUT_public_source_audit/1"})
    with urlopen(req,timeout=45) as r:
        b=r.read(limit+1)
    if len(b)>limit:
        raise ValueError("Small public metadata unexpected size")
    return b


def read_range(url,spec,limit=MAX_PREFIX):
    check_url(url)
    req=Request(url,headers={
        "User-Agent":"IS_Analysis_V3_ImmuNexUT_public_source_audit/1",
        "Range":"bytes="+spec,"Accept-Encoding":"identity",
    })
    with urlopen(req,timeout=80) as r:
        data=r.read(limit+1)
        content_range=r.headers.get("Content-Range","")
        status=r.status
    if status!=206 or len(data)>limit:
        raise ValueError("Official archive does not support required bounded byte ranges")
    m=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",content_range)
    if not m:
        raise ValueError("Source HTTP Content-Range absent or malformed")
    a,b,total=map(int,m.groups())
    if b-a+1!=len(data) or total!=EXPECTED_ZIP_SIZE:
        raise ValueError("Published source size / range mismatch")
    if spec.startswith("-"):
        if b!=total-1 or b-a+1!=min(total,int(spec[1:])):
            raise ValueError("Source suffix bytes mismatch")
    else:
        st,end=map(int,spec.split("-"))
        if (a,b)!=(st,end):
            raise ValueError("Source explicit byte range mismatch")
    return data,{"start":a,"end":b,"total":total}


def source_filelist(raw):
    lines=raw.decode("utf-8","strict").splitlines()
    meta={}
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        cols=line.split("\t")
        if len(cols)!=5:
            raise ValueError("Unexpected DDBJ source manifest TSV columns")
        kind,name,modified,size,md5=cols
        if not re.fullmatch("[0-9a-f]{32}",md5) or not size.isdecimal():
            raise ValueError("Published source metadata malformed")
        meta[name]={"type":kind,"size_bytes":int(size),"publisher_MD5":md5}
    if meta.get("E-GEAD-420.processed.zip",{}).get("size_bytes")!=EXPECTED_ZIP_SIZE:
        raise ValueError("Publisher reported ZIP bytes changed")
    if meta["E-GEAD-420.processed.zip"]["publisher_MD5"]!=EXPECTED_PUBLISHED_ZIP_MD5:
        raise ValueError("Publisher ZIP MD5 manifest drift")
    if meta.get("nominal_eQTL_2.tar.gz",{}).get("size_bytes")!=43495814577:
        raise ValueError("Publisher nominal nested TAR size changed")
    return meta


def zip64_eocd(tail,meta):
    e=tail.rfind(b"PK\x05\x06")
    if e<0 or e+22>len(tail):
        raise ValueError("No ZIP end-of-central-directory record")
    signature,disk,cd_disk,disc_entries,entries,cd_bytes,cd_offset,comment=struct.unpack_from("<IHHHHIIH",tail,e)
    if (disk!=0 or cd_disk!=0 or e+22+comment>len(tail)):
        raise ValueError("Multi-disk or malformed ZIP EOCD")
    needs64=(entries==0xFFFF or cd_bytes==0xFFFFFFFF or cd_offset==0xFFFFFFFF)
    if needs64:
        l=e-20
        if l<0 or tail[l:l+4]!=b"PK\x06\x07":
            raise ValueError("ZIP64 locator missing")
        marker,zip64_disk,record_offset,ndisks=struct.unpack_from("<IIQI",tail,l)
        if zip64_disk!=0 or ndisks!=1:
            raise ValueError("Multi-disk ZIP64 archive unsupported")
        ix=record_offset-meta["start"]
        if not 0<=ix<=len(tail)-56:
            raise ValueError("ZIP64 record outside expected bounded tail")
        vals=struct.unpack_from("<IQHHIIQQQQ",tail,ix)
        if vals[0]!=0x06064B50 or vals[1]<44 or vals[4]!=0 or vals[5]!=0 or vals[6]!=vals[7]:
            raise ValueError("Invalid ZIP64 eocd")
        entries,cd_bytes,cd_offset=vals[7],vals[8],vals[9]
    if not 1<=entries<=5000 or not 0<cd_bytes<=MAX_CD or not 0<=cd_offset<meta["total"]:
        raise ValueError("Implausible source ZIP central directory")
    if cd_offset+cd_bytes>meta["total"]:
        raise ValueError("Central directory exceeds source size")
    return {
        "central_directory_start":cd_offset,
        "central_directory_end_inclusive":cd_offset+cd_bytes-1,
        "central_directory_byte_count":cd_bytes,
        "members_expected":entries,
        "ZIP64_footer_required":needs64,
    }


def zip_directory_members(blob,n):
    items=[]
    pos=0
    for _ in range(n):
        if pos+46>len(blob):
            raise ValueError("Truncated original ZIP central directory")
        f=struct.unpack_from("<IHHHHHHIIIHHHHHII",blob,pos)
        if f[0]!=0x02014B50:
            raise ValueError("Malformed central directory entry")
        namesz,extrasz,commentsz=f[10],f[11],f[12]
        end=pos+46+namesz+extrasz+commentsz
        if end>len(blob):
            raise ValueError("Central member name exceeds ZIP size")
        name=blob[pos+46:pos+46+namesz].decode("utf-8","replace")
        extra=blob[pos+46+namesz:pos+46+namesz+extrasz]
        comp,uncomp,offset=f[8],f[9],f[16]
        # ZIP64 extra field 0x0001 values occur in order of 32-bit sentinel flags
        j=0
        while j+4<=len(extra):
            eid,esize=struct.unpack_from("<HH",extra,j)
            field=extra[j+4:j+4+esize]
            if j+4+esize>len(extra):
                raise ValueError("Corrupted ZIP extra field")
            if eid==1:
                k=0
                def next64():
                    nonlocal k
                    if k+8>len(field):
                        raise ValueError("Incomplete ZIP64 member offset/size")
                    result=struct.unpack_from("<Q",field,k)[0]
                    k+=8
                    return result
                if uncomp==0xFFFFFFFF:uncomp=next64()
                if comp==0xFFFFFFFF:comp=next64()
                if offset==0xFFFFFFFF:offset=next64()
            j+=4+esize
        if comp==0xFFFFFFFF or uncomp==0xFFFFFFFF or offset==0xFFFFFFFF:
            raise ValueError("Unresolved ZIP64 member size")
        items.append({"name":name,"method":f[4],"compressed_size":comp,
                      "uncompressed_size":uncomp,"local_header_offset":offset,
                      "original_crc32":f[7],"source_zip_content_not_fully_verified":True})
        pos=end
    if pos!=len(blob):
        raise ValueError("Unaccounted bytes after ZIP directory")
    return items


def source_member_prefix(entry,cap=MAX_PREFIX):
    h,info=read_range(ZIP,f"{entry['local_header_offset']}-{entry['local_header_offset']+29}",30)
    x=struct.unpack_from("<IHHHHHIIIHH",h,0)
    if x[0]!=0x04034B50 or x[3]!=entry["method"]:
        raise ValueError("Original ZIP local header invalid")
    namelen,exlen=x[-2:]
    nh,_=read_range(ZIP,f"{entry['local_header_offset']+30}-{entry['local_header_offset']+29+namelen+exlen}",16384)
    if nh[:namelen].decode("utf-8","replace")!=entry["name"]:
        raise ValueError("Central and local member name mismatch")
    start=entry["local_header_offset"]+30+namelen+exlen
    size=min(entry["compressed_size"],cap)
    raw,_=read_range(ZIP,f"{start}-{start+size-1}",cap)
    if entry["method"]==0:
        first=raw
    elif entry["method"]==8:
        first=zlib.decompressobj(-15).decompress(raw,2*MAX_PREFIX)
    else:
        raise ValueError("Unsupported original ZIP member compression")
    if not first or len(first)>2*MAX_PREFIX:
        raise ValueError("Unable to decode bounded member prefix")
    return first,{"source_compressed_member_prefix_bytes":size,
                  "first_decompressed_zip_member_bytes":len(first),
                  "full_43GB_member_downloaded":False}


def nested_gzip_tar_head(raw):
    if not raw.startswith(b"\x1f\x8b"):
        raise ValueError("Nested original nominal eQTL is not gzip")
    nested=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(raw,2*MAX_PREFIX)
    if len(nested)<512:
        raise ValueError("Cannot decode nested TAR entry")
    fields=[]
    offset=0
    while offset+512<=len(nested) and len(fields)<12:
        row=nested[offset:offset+512]
        if row==b"\0"*512:break
        name=row[:100].split(b"\0",1)[0].decode("utf-8","replace")
        size=row[124:136].strip(b"\0 ")
        if not re.fullmatch(rb"[0-7]+",size):
            raise ValueError("Nested TAR member has bad size")
        n=int(size,8)
        fields.append({"member":name,"uncompressed_bytes":n,"contents_not_validated":True})
        next_offset=offset+512+((n+511)//512)*512
        if next_offset>len(nested):
            break
        offset=next_offset
    if not fields:
        raise ValueError("No TAR members inspected")
    return fields



IMMUNEXUT_13_COL_HEADER=(
    "Gene_id","Gene_name","CHR","TSS_position","Number_of_variants_cis",
    "Variant_ID","REF","ALT","Variant_CHR","Variant_position_start",
    "Variant_position_end","nominal_P_value","slope(ALT)",
)


def validate_real_immuNexUT_13col_rows(pieces):
    """Validate genuinely observed 13-column ImmuNexUT source, not QTLtools guess.

    ImmuNexUT publishes a postprocessed table with ALT-oriented slope,
    nominal p, cis-test count, REF/ALT and coordinate. It does not publish
    SE/MAF/per-celltype sample n in this original observed table header.
    Multiple rsIDs and multi-base indel alleles are permitted, not collapsed.
    """
    import math
    if not pieces or tuple(pieces[0])!=IMMUNEXUT_13_COL_HEADER:
        raise ValueError("ImmuNexUT actual source 13-column header mismatch")
    if len(pieces)<4:
        raise ValueError("Less than 3 actual association rows")
    summary={
        "source_13col_header_verified":True,
        "n_actual_complete_source_association_rows":len(pieces)-1,
        "source_nominal_non_significant_p_gt_0_05":0,
        "source_nominal_significant_p_le_0_05":0,
        "source_reported_cis_variants_for_first_gene":None,
        "source_first_gene_id":None,
        "distinct_genes_in_bounded_prefix":None,
        "variant_id_semicolon_multi_ids_rows":0,
        "SNVs_and_indels_supported":True,
        "beta_effect_allele":"ALT",
        "slope_beta_present":True,
        "nominal_p_value_present":True,
        "ref_alt_alleles_present":True,
        "variant_coordinates_present":True,
        "slope_se_column_present":False,
        "allele_frequency_MAF_EAF_present":False,
        "per_celltype_sample_n_column_present":False,
        "no_association_SE_reconstruction_performed":True,
        "full_study_CRC_MD5_checked":False,
        "full_cis_association_completeness":"NOT_VERIFIED",
        "GWAS_liftover_and_harmonized_variants":"NOT_RUN",
        "coloc_ABF_or_SuSiE":"BLOCKED",
    }
    unique_genes=set()
    for i,cols in enumerate(pieces[1:],start=1):
        if len(cols)!=len(IMMUNEXUT_13_COL_HEADER):
            raise ValueError("Original association row column count mismatch")
        (gene,gname,gchr,tss,n_cis,rsids,ref,alt,vchr,vstart,vend,pvalue,beta)=cols
        if not gene.startswith("ENSG") or not gchr.startswith("chr") or not vchr.startswith("chr"):
            raise ValueError("Unexpected original gene/chr provenance")
        if gchr!=vchr:
            raise ValueError("Original cis variant chromosome differs from gene")
        if not gname or not rsids:
            raise ValueError("Source gene or variant identifier missing")
        if not set(ref)<=set("ACGT") or not set(alt)<=set("ACGT") or not ref or not alt or ref==alt:
            raise ValueError("Invalid or identical REF ALT alleles")
        try:
            pos=int(vstart); end=int(vend); transcription=int(tss); n=int(n_cis)
            pv=float(pvalue); b=float(beta)
        except ValueError as exc:
            raise ValueError("Invalid original ImmuNexUT numeric association") from exc
        if not (n>0 and pos>0 and end>=pos and transcription>0 and
                math.isfinite(pv) and 0<=pv<=1 and math.isfinite(b)):
            raise ValueError("Original ImmuNexUT association values out of bounds")
        unique_genes.add(gene)
        summary["source_nominal_non_significant_p_gt_0_05"]+=(pv>0.05)
        summary["source_nominal_significant_p_le_0_05"]+=(pv<=0.05)
        summary["variant_id_semicolon_multi_ids_rows"]+=(";" in rsids)
        if i==1:
            summary["source_reported_cis_variants_for_first_gene"]=n
            summary["source_first_gene_id"]=gene
    summary["distinct_genes_in_bounded_prefix"]=len(unique_genes)
    return summary


def first_gene_full_cis_block_probe(gzip_tar_member_prefix):
    """Look for the first COMPLETE gene block in a bounded nested gzip prefix.

    Requires seeing a following gene (or full end-of-member) to establish that
    the first gene's 835 reported rows are bounded and not just 80 records.
    Count matching confirms *first gene block row completeness* against its
    own source metadata, not reference-genome full cis variant coverage.
    """
    if not gzip_tar_member_prefix.startswith(b"\x1f\x8b"):
        raise ValueError("Nested gene association stream not gzip")
    decompressor=zlib.decompressobj(16+zlib.MAX_WBITS)
    buffer=decompressor.decompress(gzip_tar_member_prefix,2*MAX_PREFIX)
    if len(buffer)<1536:
        raise ValueError("No complete gene association block prefix")
    hdr=buffer[:512]
    filename=hdr[:100].split(b"\0",1)[0].decode("utf-8","replace")
    if not filename.endswith("_nominal.txt"):
        raise ValueError("First nested file is not nominal SNP-gene associations")
    body=buffer[512:]
    lines=body.split(b"\n")
    # Last decoded line may be partial; do not treat it as a valid record.
    complete_lines=lines[:-1]
    if not complete_lines:
        raise ValueError("No complete source eQTL gene rows in prefix")
    header=complete_lines[0].decode("utf-8","strict").rstrip("\r").split()
    if tuple(header)!=IMMUNEXUT_13_COL_HEADER:
        raise ValueError("Original 13-column file schema changed")
    first_gene=None
    expected=None
    nrows=0
    n_nonsig=0
    loci=set()
    saw_next_gene=False
    for raw_line in complete_lines[1:]:
        if not raw_line.strip():
            continue
        row=raw_line.decode("utf-8","strict").rstrip("\r").split()
        if len(row)!=13:
            raise ValueError("Malformed original cis gene block row")
        gene=row[0]
        if first_gene is None:
            first_gene=gene
            expected=int(row[4])
            if expected<1 or expected>2000000:
                raise ValueError("Unreasonable declared original first-gene cis variant count")
        if gene!=first_gene:
            saw_next_gene=True
            break
        if int(row[4])!=expected:
            raise ValueError("First gene cis tested count changes across rows")
        nrows+=1
        try:
            p=float(row[11])
        except ValueError as exc:
            raise ValueError("Source nominal pvalue nonnumeric") from exc
        if not 0<=p<=1:
            raise ValueError("Source nominal pvalue out of range")
        n_nonsig+=(p>0.05)
        locus=(row[8],row[9],row[10],row[6],row[7])
        if locus in loci:
            raise ValueError("Repeated original variant coordinate/REF/ALT in first gene")
        loci.add(locus)
    if first_gene is None:
        raise ValueError("No first source gene found")
    result={
        "source_gene_id":first_gene,
        "source_declared_cis_tested_variants":expected,
        "original_gene_association_rows_seen":nrows,
        "distinct_original_variant_chrom_pos_allele_keys":len(loci),
        "source_nonsignificant_p_gt_0_05_rows":n_nonsig,
        "observed_start_of_second_gene":saw_next_gene,
        "fully_observed_first_gene_block":saw_next_gene,
        "number_of_cis_variants_matches_first_gene_association_rows":nrows==expected if saw_next_gene else None,
        "original_source_metadata_vs_row_count":(
            "COUNT_MATCHED" if saw_next_gene and nrows==expected
            else "COUNT_MISMATCH" if saw_next_gene
            else "UNFINISHED_GENE_PREFIX_HOLD"),
        "raw_variant_gwas_alleles_harmonized":False,
        "full_genome_cis_catalogue_covered":False,
        "colocalization":"NOT_RUN",
        "source_publication_claim_of_unfiltered_cis":"PUBLISHED_NOT_ALL_GENES_INDEPENDENTLY_VERIFIED",
    }
    return result


def first_nominal_tsv_preview(gzip_tar_member_prefix,max_rows=80):
    """Inspect only start of the first ~8GB nominal QTLtools TAR member.

    This is source-byte schema evidence, not full member CRC or sample N.
    QTLtools original nominal output can be headerless, so no guessed
    column names are hard-coded as established evidence.
    """
    if not gzip_tar_member_prefix.startswith(b"\x1f\x8b"):
        raise ValueError("Not an original gzip member")
    inflated=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(
        gzip_tar_member_prefix,2*MAX_PREFIX)
    if len(inflated)<1536:
        raise ValueError("Too few bytes to inspect original nominal rows")
    header=inflated[:512]
    name=header[:100].split(b"\0",1)[0].decode("utf-8","replace")
    if not name.endswith("_nominal.txt"):
        raise ValueError("First nested nominal TAR member not association TSV")
    declared=header[124:136].strip(b"\0 ")
    if not re.fullmatch(rb"[0-7]+",declared):
        raise ValueError("Original nominal TAR member size malformed")
    declared_size=int(declared,8)
    if declared_size<100:
        raise ValueError("Nominal TSV member unexpectedly small")
    raw=inflated[512:512+min(declared_size,len(inflated)-512)]
    decoded=raw.decode("utf-8","strict")
    rows=decoded.split("\n")
    if len(rows)<4:
        raise ValueError("Original association source has fewer than 3 full lines")
    first=rows[:max_rows] if decoded.endswith("\n") else rows[:max_rows]
    complete=[x.rstrip("\r") for x in first if x.strip()]
    # Exclude final partial line if it lies at the decompression boundary.
    if not decoded.endswith("\n") and len(rows)<=max_rows:
        complete=complete[:-1]
    if len(complete)<3:
        raise ValueError("Too few complete original nominal associations")
    pieces=[x.split() for x in complete]
    num_cols=[len(p) for p in pieces]
    if not 7<=num_cols[0]<=40 or any(count!=num_cols[0] for count in num_cols):
        raise ValueError("Inconsistent association TSV field counts in source")
    source_13col_stats=None
    if pieces[0] and pieces[0][0]=="Gene_id":
        source_13col_stats=validate_real_immuNexUT_13col_rows(pieces)
    sample=pieces[:8]
    # Upstream QTLtools cis --nominal format is WHITESPACE-separated.
    # slope_se exists only when --std-err was supplied. Optional group
    # columns move positional indices: do not guess full schema if unexpected.
    header_fields=[v.lower().lstrip("#") for v in sample[0]]
    source_has_header=("phe_id" in header_fields or "grp_id" in header_fields)
    data_rows=pieces[1:] if source_has_header else pieces
    qtltools_schema=None
    if not source_has_header and len(data_rows[0]) in (15,16):
        first_data=data_rows[0]
        # Standard ungrouped positional QTLtools fields; prior to declaring
        # a genuine format, require chromosome and nominal p-value semantics.
        try:
            n_cis=int(first_data[5])
            qtl_p=float(first_data[11])
            qtl_slope=float(first_data[13])
            chr1=first_data[1]
            chr2=first_data[8]
            if n_cis>0 and 0<=qtl_p<=1 and chr1==chr2 and qtl_slope==qtl_slope and first_data[-1] in ("0","1"):
                qtltools_schema={
                    "format":"QTLtools nominal ungrouped positional candidate",
                    "nominal_pvalue":qtl_p,
                    "n_variants_tested_cis_for_first_gene":n_cis,
                    "first_phenotype_id":first_data[0],
                    "first_variant_id":first_data[7],
                    "slope_beta":qtl_slope,
                    "slope_se_present":len(first_data)==16,
                    "slope_se_column_option_required":"--std-err",
                    "first_record_without_confirmed_REF_ALT":True,
                    "first_record_without_proven_per_celltype_N":True,
                }
                if len(first_data)==16:
                    se=float(first_data[14])
                    if not 0<se<float("inf"):
                        raise ValueError("Nonpositive QTLtools slope_se")
                    qtltools_schema["slope_se"]=se
        except (ValueError,IndexError):
            pass
    header_guess=any(word in [v.strip().lower() for v in sample[0]]
                     for word in ("variant_id","slope","slope_se","pval_nominal","phenotype_id","gene_id"))
    numeric_by_col=[]
    for j in range(num_cols[0]):
        inspected=pieces[1:40] if header_guess else pieces[:40]
        numeric=0
        for row in inspected:
            try:
                value=float(row[j])
                if value==value and value not in (float("inf"),float("-inf")):
                    numeric+=1
            except ValueError:
                pass
        numeric_by_col.append({"column_index_zero_based":j,
                               "numeric_fraction_first_rows":round(numeric/max(1,len(inspected)),3)})
    return {
        "source_nominal_file":name,
        "declared_first_nominal_uncompressed_bytes":declared_size,
        "sampled_complete_association_lines":len(complete),
        "whitespace_delimited_columns_per_line":num_cols[0],
        "tab_delimiter_required":False,
        "qtltools_nominal_ungrouped_schema_probe":qtltools_schema,
        "actual_ImmuNexUT_13_column_source_data_QC":source_13col_stats,
        "looks_like_header":header_guess,
        "first_line_fields_verbatim":sample[0][:40],
        "second_line_fields_verbatim":sample[1][:40],
        "numeric_field_evidence":numeric_by_col,
        "no_source_per_person_data_extracted":True,
        "full_nominal_member_size_or_crc_checked":False,
        "entire_cis_distribution_examined":False,
        "coloc_allele_sample_size_compatibility":"NOT_YET_VERIFIED",
    }


def run(out):
    report={
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "dataset":"ImmuNexUT E-GEAD-420 nominal QTLtools cis eQTL all data",
        "source_publication":"NBDC hum0214, 2021, Japanese 416 participants across disease/healthy cohorts",
        "association_source":"https://humandbs-production.ddbj.nig.ac.jp/en/dataset/E-GEAD-420",
        "source_GWAS_build":"GIGASTROKE GRCh37",
        "nominal_source_build":"GRCh38",
        "source_status":"SOURCE_OR_QC_HOLD",
        "published_whole_ZIP_MD5_verified_by_full_download":False,
        "source_filelist_publisher_MD5_only":True,
        "full_genome_cis_beta_SE_EAF_rows_read":False,
        "GWAS_liftover":"NOT_RUN",
        "coloc":"NOT_RUN",
        "disease_covar_and_per_celltype_N_checked":False,
    }
    try:
        filelist=get_small(FILELIST,6000)
        manifest=source_filelist(filelist)
        report["published_manifest"]=manifest
        idf=get_small(IDF,12000)
        if len(idf)!=manifest["E-GEAD-420.idf.txt"]["size_bytes"]:
            raise ValueError("Original IDF size differs")
        if hashlib.md5(idf).hexdigest()!=manifest["E-GEAD-420.idf.txt"]["publisher_MD5"]:
            raise ValueError("Original IDF publisher MD5 mismatch")
        report["small_IDF_original_MD5"]="PASS"
        report["IDF_field_names"]=[row.split("\t",1)[0] for row in idf.decode("utf-8","replace").splitlines() if row]
        tail,window=read_range(ZIP,f"-{MAX_TAIL}",MAX_TAIL)
        footer=zip64_eocd(tail,window)
        report["ZIP_directory_footer"]=footer
        start,end=footer["central_directory_start"],footer["central_directory_end_inclusive"]
        if start>=window["start"] and end<=window["end"]:
            directory=tail[start-window["start"]:start-window["start"]+footer["central_directory_byte_count"]]
        else:
            directory,_=read_range(ZIP,f"{start}-{end}",MAX_CD)
        members=zip_directory_members(directory,footer["members_expected"])
        report["ZIP_member_inventory"]=members[:100]
        matching=[x for x in members if x["name"].endswith("nominal_eQTL_2.tar.gz")]
        if len(matching)!=1 or matching[0]["uncompressed_size"]!=manifest["nominal_eQTL_2.tar.gz"]["size_bytes"]:
            raise ValueError("Expected original 43GB nominal TAR member not found or size conflict")
        first,source_meta=source_member_prefix(matching[0])
        report["original_nominal_member_prefix"]=source_meta
        report["nested_nominal_TAR_first_members"]=nested_gzip_tar_head(first)
        report["first_actual_nominal_association_records"]=first_nominal_tsv_preview(first)
        report["first_gene_full_cis_block_probe"]=first_gene_full_cis_block_probe(first)
        report["source_status"]="OFFICIAL_PUBLIC_43GB_ZIP64_CENTRAL_AND_NOMINAL_TAR_PREFIX_VERIFIED"
    except (OSError,ValueError,TimeoutError,urllib.error.URLError,struct.error,zlib.error) as e:
        report["source_status"]="SOURCE_OR_QC_HOLD"
        report["source_issue"]=f"{type(e).__name__}: {str(e)[:230]}"
    out.mkdir(parents=True,exist_ok=True)
    (out/"E_GEAD_420_source_zip64_evidence.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=["# ImmuNexUT E-GEAD-420 real 43GB nominal cis eQTL source","",
           "**"+report["source_status"]+"**","",
           "Source says all nominal tested cis associations, including nonsignificant, but actual full rows not downloaded.",
           ""]
    if y:=report.get("ZIP_member_inventory"):
        lines+=["| Original ZIP member | Declared member size | ZIP compression |","|---|---:|---:|"]
        for v in y[:15]:
            lines.append(f"| {v['name']} | {v['uncompressed_size']:,} | {v['method']} |")
    if y:=report.get("nested_nominal_TAR_first_members"):
        lines+=["","| First nested nominal TAR member | Declared bytes |","|---|---:|"]
        for v in y[:12]:
            lines.append(f"| {v['member']} | {v['uncompressed_bytes']:,} |")
    if v:=report.get("first_actual_nominal_association_records"):
        lines.extend(["", "## Actual original source nominal SNP-gene table preview",
                      f"- First file: {v['source_nominal_file']}",
                      f"- Complete source lines inspected: {v['sampled_complete_association_lines']}",
                      f"- Whitespace-delimited fields per row: {v['whitespace_delimited_columns_per_line']}",
                      f"- QTLtools positional schema probe: {v.get('qtltools_nominal_ungrouped_schema_probe')}",
                      f"- Actual ImmuNexUT 13-column data audit: {v.get('actual_ImmuNexUT_13_column_source_data_QC')}",
                      f"- Looks like a header: {v['looks_like_header']}",
                      "- First row: `" + " | ".join(v["first_line_fields_verbatim"])[:500] + "`",
                      "- Second row: `" + " | ".join(v["second_line_fields_verbatim"])[:500] + "`"])
    if first_gene:=report.get("first_gene_full_cis_block_probe"):
        lines.extend(["", "## First-gene complete block count cross-check",
                      f"- Gene: {first_gene['source_gene_id']}",
                      f"- Original declared cis variants: {first_gene['source_declared_cis_tested_variants']}",
                      f"- Source rows observed until next gene: {first_gene['original_gene_association_rows_seen']}",
                      f"- Next gene observed: {first_gene['observed_start_of_second_gene']}",
                      f"- Result: {first_gene['original_source_metadata_vs_row_count']}",
                      f"- Non-significant tested variants in this gene block: {first_gene['source_nonsignificant_p_gt_0_05_rows']}"])
    if err:=report.get("source_issue"):
        lines.append("\nHOLD: "+err)
    lines+=["","GRCh38 public immune cell QTL vs GRCh37 original stroke GWAS; no positional join or causal claims.",
            "Original 43GB publisher MD5 not verified (full transfer intentionally avoided).",
            "N=416 overall not guaranteed per cell type; disease covariates and donor mix need source review."]
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md,encoding="utf-8")
    print(md)
    if f:=os.getenv("GITHUB_STEP_SUMMARY"):
        with open(f,"a",encoding="utf-8") as stream:stream.write(md)
    return report


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--out-dir",default="artifacts/is_EGEAD420_pilot")
    args=a.parse_args()
    audit=run(Path(args.out_dir))
    if audit.get("status",audit.get("source_status"))!="OFFICIAL_PUBLIC_43GB_ZIP64_CENTRAL_AND_NOMINAL_TAR_PREFIX_VERIFIED":
        raise SystemExit(1)
