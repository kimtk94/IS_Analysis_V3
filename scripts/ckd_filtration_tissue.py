import argparse
import collections
import gzip
import hashlib
import io
import json
import math
import re
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

PIN="ad0d9181d3d1b4d598b7c9e5a3d5be11825db8681a877069a4410c96834328c5"
EXPECTED={"rs16846391","rs2070600","rs4390169","rs17050272","rs77924615","rs10185736","rs11162405"}
RECORD="16976585"
COMP=str.maketrans("ACGT","TGCA")
PHEWEB="https://koges.leelabsg.org"

def digest(p,alg="sha256"):
    with p.open("rb") as f: return hashlib.file_digest(f,alg).hexdigest()

def get(url,cap=2000000):
    with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"CKD-public-summary-research/1.0","Accept":"application/json, */*"}),timeout=25) as f:
        raw=f.read(cap+1)
    if len(raw)>cap: raise ValueError("Response exceeds bound")
    return raw

def anchors(root):
    p=root/"egfr.zip"
    if digest(p)!=PIN: raise ValueError("Prior eGFR artifact digest changed")
    with zipfile.ZipFile(p) as z:
        if z.testzip(): raise ValueError("Prior archive CRC")
        name=[n for n in z.namelist() if n.rsplit("/",1)[-1]=="assay_results.json"]
        if len(name)!=1: raise ValueError("Prior assay table missing")
        rows=json.loads(z.read(name[0]))
    out=collections.defaultdict(list)
    for r in rows:
        if r["status"]=="EXPLORATORY_WALD" and r["bonferroni_p"]<.05:
            out[r["rsid"]].append(r)
    if set(out)!=EXPECTED or sum(map(len,out.values()))!=10: raise ValueError("Frozen family changed")
    for rs,rr in out.items():
        a=rr[0]
        if any((x["chr37"],x["pos37"],x["ref37"],x["alt37"])!=(a["chr37"],a["pos37"],a["ref37"],a["alt37"]) for x in rr):
            raise ValueError("Anchor discordance")
    return dict(sorted(out.items()))

def orient(ea,oa,a,b):
    if {ea,oa} in ({"A","T"},{"C","G"}):return None
    if (ea,oa)==(a,b):return 1
    if (ea,oa)==(b,a):return -1
    return None

def wald(bx,sx,by,sy):
    if not all(math.isfinite(x) for x in (bx,sx,by,sy)) or bx==0 or min(sx,sy)<=0:raise ValueError("Nonfinite estimate")
    b=by/bx
    se=math.hypot(sy,b*sx)/abs(bx)
    return {"mr_beta":b,"mr_se":se,"mr_ci95":[b-1.959964*se,b+1.959964*se],
            "mr_p":math.erfc(abs(b/se)/math.sqrt(2))}

def summarize(rs,rr,by,sy,ea,oa,source,expected_egfr=True):
    a=rr[0]
    s=orient(a["ea37"],a["oa37"],ea,oa)
    if s is None:raise ValueError("Allele mismatch or palindrome")
    p=math.erfc(abs(by/sy)/math.sqrt(2))
    result={"rsid":rs,"outcome_source":source,"outcome_ea":ea,"outcome_oa":oa,
            "outcome_beta":by,"outcome_se":sy,"two_sided_p":p,
            "bonferroni_p_7":min(1,7*p),"egfr_direction_agrees":bool(by*s*a["outcome_beta_aligned"]>0) if expected_egfr else None,
            "assays":[]}
    for r in rr:
        sign=orient(r["ea37"],r["oa37"],ea,oa)
        if sign is None:raise ValueError("Assay allele mismatch")
        m=wald(r["beta"],r["se"],by*sign,sy)
        result["assays"].append({"platform":r["platform"],"assay_id":r["assay_id"],
            "exposure_beta":r["beta"],"exposure_se":r["se"],"previous_egfr_mr_beta":r["mr_beta"],
            "new_outcome_beta_aligned":by*sign,**m})
    return result

def write(out,obj):
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")
    print(json.dumps(obj,indent=2,allow_nan=False))


# Fixed seven variants and seven genes from the prior eGFR discovery union.
GENES=("FAP","INHBB","EFNA1","AGER","UMOD","CGREF1","GIPC2")
def renal(root,out,trait):
    panel=anchors(root)
    path=root/"outcome.gz"
    # Official CKDGen source schemas, build 37. No proxy or complement rescue.
    schema=({"chr":"chr","pos":"pos","id":"rsid","ea":"allele1","oa":"allele2",
             "beta":"effect","se":"stderr.gc","p":"p.value.gc","freq":"freq1","n":"n"}
            if trait=="cys" else
            {"chr":"chr","pos":"pos_b37","id":"rsid","ea":"allele1","oa":"allele2",
             "beta":"effect","se":"stderr","p":"p-value","freq":"freq1","n":"n_total_sum"})
    hits=collections.defaultdict(list)
    total=bad=0
    with gzip.open(path,"rt",encoding="utf-8") as f:
        h=f.readline().strip().split()
        hh=[x.lower() for x in h]
        if len(set(hh))!=len(hh) or not set(schema.values())<=set(hh):
            raise ValueError("Source schema mismatch "+repr(h))
        ix={k:hh.index(v) for k,v in schema.items()}
        for line in f:
            if not line.strip():continue
            total+=1
            q=line.split()
            if len(q)!=len(h):
                bad+=1
                continue
            rs=q[ix["id"]]
            if rs in panel: hits[rs].append({k:q[i] for k,i in ix.items()})
    if total<1000000 or bad>0:raise ValueError("Source row/schema failure "+str((total,bad)))
    rows=[]
    for rs,rr in panel.items():
        a=rr[0]
        valid=[]
        for r in hits[rs]:
            try:
                if (int(r["chr"].removeprefix("chr")),int(r["pos"]))!=(a["chr37"],a["pos37"]):continue
                sign=orient(a["ea37"],a["oa37"],r["ea"].upper(),r["oa"].upper())
                if sign is not None:valid.append((r,sign))
            except (ValueError,TypeError):continue
        if len(valid)!=1:
            rows.append(dict(rsid=rs,status="HOLD_UNMATCHED_OR_AMBIGUOUS",rsid_rows=len(hits[rs])))
            continue
        r,sign=valid[0]
        by,sy,p,freq=map(float,(r["beta"],r["se"],r["p"],r["freq"]))
        if not (all(math.isfinite(v) for v in (by,sy,p,freq)) and sy>0 and 0<=p<=1 and 0<freq<1):
            rows.append(dict(rsid=rs,status="HOLD_NUMERIC"));continue
        aligned=by*sign
        consistent=aligned*a["outcome_beta_aligned"]*(1 if trait=="cys" else -1)>0
        rows.append(dict(rsid=rs,status="MATCHED_SOURCE_ASSOCIATION",assays=[x["assay_id"] for x in rr],
            effect_allele=a["ea37"],other_allele=a["oa37"],beta_aligned=aligned,se=sy,
            ci95=[aligned-1.959964*sy,aligned+1.959964*sy],source_p=p,bonferroni_p_7=min(1,p*7),
            frequency_aligned=freq if sign==1 else 1-freq,source_n=r["n"],
            expected_direction_agrees=consistent,passes=consistent and p*7<.05))
    write(out,dict(status="SUPPORTING_CROSS_ANCESTRY_SNP_ASSOCIATIONS_NOT_MR",trait=trait,
         source_bytes=path.stat().st_size,source_sha256=digest(path),header=h,rows_scanned=total,
         malformed=bad,family=sorted(panel),results=rows,
         limit="Predominantly European renal outcomes, Chinese exposure not re-estimated. SNP associations, not ancestry-matched MR or independent tests across overlapping endpoints. No publisher checksum available."))

def korean_probe(out):
    # Request only the first compressed megabyte, never treat partial content as a verified GWAS.
    url="https://zenodo.org/api/records/16976585/files/GWAS_summary_statistics.txt.gz/content"
    import zlib
    req=urllib.request.Request(url,headers={"Range":"bytes=0-1048575","User-Agent":"CKD-public-summary-research/1.0","Accept-Encoding":"identity"})
    with urllib.request.urlopen(req,timeout=30) as f:
        b=f.read(1048576)
        status=f.status
        content_range=f.headers.get("Content-Range")
    raw=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(b,250000)
    lines=raw.decode("utf-8",errors="replace").splitlines()[:8]
    # Only aggregate GWAS header and first seven source records, no individual data.
    write(out,dict(status="PARTIAL_SOURCE_SCHEMA_ONLY",http_status=status,content_range=content_range,
         compressed_bytes_sampled=len(b),first_lines=lines,whitespace_fields=[len(x.split()) for x in lines],
         tab_fields=[len(x.split("\t")) for x in lines],comma_fields=[len(x.split(",")) for x in lines]))

def table(path,start,end,gz=True):
    op=gzip.open if gz else open
    with op(path,"rt",encoding="utf-8-sig",errors="strict") as f:
        pre=[]
        for line in f:
            if line.strip().lower()==start:break
            pre.append(line)
        else:raise ValueError("Missing table start "+start)
        body=[]
        for line in f:
            if line.strip().lower()==end:break
            body.append(line)
        else:raise ValueError("Missing table end "+end)
    import csv
    return pre,csv.reader(body,delimiter="\t")

def biopsy(root,out):
    from scipy.stats import rankdata,mannwhitneyu
    import statistics
    from collections import Counter
    pre,a=table(root/"platform.txt","!platform_table_begin","!platform_table_end",False)
    if not any(x.strip()=="^PLATFORM = GPL14951" for x in pre):raise ValueError("Wrong platform")
    h=next(a)
    if not {"ID","Symbol","ILMN_Gene","Species"}<=set(h):raise ValueError("Platform schema "+repr(h))
    mapping={}
    seen=set()
    for row in a:
        if len(row)!=len(h):raise ValueError("Malformed platform row")
        r=dict(zip(h,row));pid=r["ID"]
        if pid in seen:raise ValueError("Duplicate probe")
        seen.add(pid)
        g=r["Symbol"].strip()
        if g in GENES and r["ILMN_Gene"].strip()==g and r["Species"].strip().casefold()=="homo sapiens":
            mapping[pid]=g
    pre,a=table(root/"matrix.gz","!series_matrix_table_begin","!series_matrix_table_end")
    import csv
    h=next(a);ids=h[1:]
    expected={}
    for start,n,prefix,group in [(4301,22,"minimal_lesion_IgAN_","minimal"),(4323,8,"active_lesion_IgAN_","active"),
                                  (4331,12,"chronic_lesion_IgAN_","chronic"),(4343,10,"MIXED_lesions_IgAN_","mixed"),
                                  (4353,22,"non IgAN GN_","non_IgAN"),(4375,7,"LD_","LD")]:
        for i in range(n):expected["GSM324"+str(start+i)]=(prefix+str(i+1),group)
    if h[0]!="ID_REF" or set(ids)!=set(expected) or len(ids)!=81:raise ValueError("Wrong series sample set")
    titles=[next(csv.reader([x],delimiter="\t"))[1:] for x in pre if x.startswith("!Sample_title\t")]
    if len(titles)!=1 or titles[0]!=[expected[i][0] for i in ids]:raise ValueError("GEO title mismatch")
    groups=[expected[i][1] for i in ids]
    probes=collections.defaultdict(dict);seen=set();n=0;incomplete=[]
    for row in a:
        if len(row)!=len(h) or row[0] in seen:raise ValueError("Matrix row schema")
        pid=row[0];seen.add(pid);n+=1
        if pid not in mapping:continue
        try:
            v=list(map(float,row[1:]))
            if not all(math.isfinite(x) for x in v):raise ValueError("Nonfinite")
        except ValueError:incomplete.append(pid);continue
        # Relative within-probe rank across all 81 biopsies; no assumption of log2 scale.
        probes[mapping[pid]][pid]=list(rankdata(v,method="average")/len(v))
    if n!=29366:raise ValueError("Matrix row count differs "+str(n))
    rows=[]
    for gene in GENES:
        if not probes[gene]:
            rows.append(dict(gene=gene,status="HOLD_NO_COMPLETE_MAPPED_PROBE"));continue
        values=[statistics.median(xs) for xs in zip(*probes[gene].values())]
        x=[v for v,g in zip(values,groups) if g=="chronic"]
        y=[v for v,g in zip(values,groups) if g=="minimal"]
        u,p=mannwhitneyu(x,y,alternative="two-sided",method="asymptotic")
        rows.append(dict(gene=gene,status="UNADJUSTED_BIOPSY_ASSOCIATION",probes=sorted(probes[gene]),
             medians={g:statistics.median(v for v,k in zip(values,groups) if k==g) for g in sorted(set(groups))},
             chronic_n=len(x),minimal_n=len(y),rank_biserial=float(2*u/(len(x)*len(y))-1),
             p=float(p),bonferroni_p_7=min(1,float(p)*7),
             probe_rank_biserial={pid:float(2*mannwhitneyu([v for v,g in zip(vs,groups) if g=="chronic"],[v for v,g in zip(vs,groups) if g=="minimal"],method="asymptotic").statistic/(len(x)*len(y))-1) for pid,vs in probes[gene].items()}))
    write(out,dict(status="EXPLORATORY_IGAN_BIOPSY_CHRONIC_VS_MINIMAL",source="GSE116626 GPL14951",
         platform_sha256=digest(root/"platform.txt"),matrix_sha256=digest(root/"matrix.gz"),matrix_rows=n,
         sample_group_counts=dict(Counter(groups)),incomplete_candidate_probes=incomplete,results=rows,
         limit="One reported biopsy per participant, 12 chronic vs 22 minimal; unadjusted ordinal associations, potential age/sex/center/eGFR/cell mixture confounding. Published processed values, not fold change or tissue protein; no genetic mediation."))

def hpa(root,out):
    import csv
    with zipfile.ZipFile(root/"hpa.zip") as z:
        if z.testzip():raise ValueError("HPA ZIP CRC")
        names=[n for n in z.namelist() if n.rsplit("/",1)[-1]=="rna_single_cell_type_tissue.tsv"]
        if len(names)!=1:raise ValueError("HPA source table missing")
        with z.open(names[0]) as f:
            rd=csv.DictReader(io.TextIOWrapper(f,encoding="utf-8-sig"),delimiter="\t")
            h=rd.fieldnames
            req={"Gene","Gene name","Tissue","Cell type","Cluster","nTPM"}
            if not req<=set(h or []):raise ValueError("HPA schema "+repr(h))
            rows=[];total=0;groups=set();tissues=set()
            for r in rd:
                total+=1;tissues.add(r["Tissue"])
                if r["Tissue"].casefold()!="kidney":continue
                groups.add((r["Cell type"],r["Cluster"]))
                if r["Gene name"] not in GENES:continue
                v=float(r["nTPM"])
                if not math.isfinite(v) or v<0:raise ValueError("Invalid nTPM")
                rows.append({"gene":r["Gene name"],"ensembl":r["Gene"],"cell_type":r["Cell type"],"cluster":r["Cluster"],"nTPM":v})
    if not groups or not rows:raise ValueError("No exact kidney dataset "+repr(sorted(tissues)))
    for g in GENES:
        rr=[r for r in rows if r["gene"]==g]
        if len(rr)!=len(groups) or len({(r["cell_type"],r["cluster"]) for r in rr})!=len(groups) or len({r["ensembl"] for r in rr})!=1:
            raise ValueError("Incomplete gene-cluster mapping "+g)
    write(out,dict(status="NORMAL_REFERENCE_CLUSTER_EXPRESSION_DESCRIPTIVE",version="HPA 22.0",
        source_sha256=digest(root/"hpa.zip"),source_rows=total,kidney_clusters=len(groups),
        top_clusters={g:sorted([r for r in rows if r["gene"]==g],key=lambda r:-r["nTPM"])[:5] for g in GENES},
        rows=rows,limit="Cluster-aggregated normal reference RNA, not donor-level disease differential expression or protein abundance."))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=["cys","ckd","korean-probe","biopsy","hpa"])
    p.add_argument("--root",type=Path,default=Path("."))
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if a.mode in ("cys","ckd"):renal(a.root,a.out,a.mode)
    elif a.mode=="korean-probe":korean_probe(a.out)
    elif a.mode=="biopsy":biopsy(a.root,a.out)
    else:hpa(a.root,a.out)
