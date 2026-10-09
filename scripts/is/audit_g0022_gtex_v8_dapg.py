#!/usr/bin/env python3
"""GTEx v8 DAP-G fine-mapping SNP overlap with GRCh37 EAS AIS G0022.
The GTEx v8 endpoint is GRCh38. STRICT liftover+alleles required.
Variant-level PIP overlap is NOT molecular-QTL colocalization.
"""
import argparse,csv,gzip,hashlib,json,os,re,sys,time
from collections import Counter,defaultdict
from pathlib import Path
from urllib.parse import urlparse

VENDOR=Path("/srv/is-analysis/data/is/reference/vendor_pyliftover")
if VENDOR.exists():sys.path.insert(0,str(VENDOR))
try:
    from pyliftover import LiftOver
except ImportError:
    LiftOver=None  # execution requires isolated pyliftover install; unit tests use mock

ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
CHAIN=Path("/srv/is-analysis/data/reference/hg38ToHg19.over.chain")
GENES=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/IS_ALL_GENE_WINDOW_UNIVERSE.tsv")
OUT=ROOT/"gtex_v8_dapg_liftover"
SYMBOLS=("ACAD10","ALDH2","NAA25","HECTD4","BRAP","PTPN11","IFT81","ATP2A2")
BASE="https://gtexportal.org/api/v2"
MAX_TOTAL_PER_GENE=25000
def stable(x):return str(x).split(".")[0]
def read_tsv(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def complement(a):
    return a.translate(str.maketrans("ACGT","TGCA"))
def variant_liftover(variant,converter):
    m=re.fullmatch(r"chr(\d+)_([1-9]\d*)_([ACGT]+)_([ACGT]+)_b38",str(variant))
    if not m:return None,"BAD_GTEX_VARIANT_ID"
    ch,pos,ref,alt=m.groups()
    if ch!="12" or len(ref)!=1 or len(alt)!=1:return None,"CHROM_OR_NON_SNP"
    mapped=converter.convert_coordinate("chr"+ch,int(pos)-1)
    if len(mapped)!=1:return None,"LIFTOVER_AMBIGUOUS_OR_UNMAPPED"
    c,point,strand,score=mapped[0]
    if c!="chr12" or int(point)<0:return None,"LIFTOVER_OTHER_CHROMOSOME"
    if point!=int(point):return None,"NON_INTEGER_MAPPING"
    if strand=="-":ref,alt=complement(ref),complement(alt)
    elif strand!="+":return None,"BAD_MAPPING_STRAND"
    return (f"12:{int(point)+1}:{ref}:{alt}",str(strand)),"LIFTOVER_ONE_TO_ONE"
def make_json_file(path,obj):
    tmp=Path(str(path)+".part")
    tmp.write_text(json.dumps(obj,ensure_ascii=False,separators=(",",":")))
    tmp.replace(path)
def request_json(session,endpoint,params):
    import requests
    url=BASE+endpoint
    if urlparse(url).netloc!="gtexportal.org" or endpoint not in (
        "/reference/gene","/association/fineMapping"):
        raise RuntimeError("Source endpoint not allowlisted")
    err=None
    for attempt in range(3):
        try:
            r=session.get(url,params=params,timeout=(10,30))
            r.raise_for_status()
            obj=r.json()
            if "data" not in obj or not isinstance(obj["data"],list):
                raise ValueError("GTEx API data schema invalid")
            return obj
        except (requests.RequestException,ValueError) as e:
            err=e
            if attempt<2:time.sleep(attempt+1)
    raise RuntimeError(f"GTEx endpoint unavailable {endpoint}: {err}")
def fetch_gene(session,gene,cache):
    lookup=cache/(gene+".v26.gene.json")
    if lookup.exists():
        raw=json.loads(lookup.read_text())
    else:
        raw=request_json(session,"/reference/gene",{"geneId":gene,"gencodeVersion":"v26",
             "genomeBuild":"GRCh38/hg38","itemsPerPage":5})
        make_json_file(lookup,raw)
    hits=[x for x in raw["data"] if x["geneSymbol"].upper()==gene and
          x["genomeBuild"]=="GRCh38/hg38"]
    if len(hits)!=1:raise ValueError("GTEx gene not uniquely resolved "+gene)
    gencode=hits[0]["gencodeId"]
    entries=[];page=0
    while True:
        path=cache/(gene+f".v8.dapg.page{page:03d}.json")
        if path.exists():
            data=json.loads(path.read_text())
        else:
            data=request_json(session,"/association/fineMapping",{
                "gencodeId":gencode,"datasetId":"gtex_v8",
                "page":page,"itemsPerPage":2500})
            make_json_file(path,data)
        paging=data.get("paging_info",{})
        n=int(paging.get("totalNumberOfItems",len(data["data"])))
        pages=int(paging.get("numberOfPages",1))
        if n>MAX_TOTAL_PER_GENE or pages>30:raise ValueError("GTEx page cap exceeded; fail closed")
        entries.extend(data["data"])
        if page+1>=pages:break
        if not data["data"]:raise ValueError("Unexpected truncated GTEx pagination")
        page+=1
    if len(entries)!=n:raise ValueError(f"GTEx pagination mismatch {gene}: {len(entries)} vs {n}")
    return gencode,entries
def fetch(args):
    if not args.chain.is_file():raise FileNotFoundError("Need official hg38ToHg19 chain")
    gout=read_tsv(args.genes)
    src={r["gene_symbol"]:r for r in gout if r["group_id"]=="IS_XDATA_G0022"}
    if not all(x in src for x in SYMBOLS):raise ValueError("GENCODE v19 positional symbols missing")
    gw=read_tsv(args.root/"variants.tsv")
    pips={x["variant_id"]:x for x in gw}
    if len(gw)<4500:raise ValueError("Full-locus reference genotype data unavailable")
    v8=json.loads((args.root/"G0022_FULL_AIS_SUSIE_SUMMARY.json").read_text())
    target={s for sets in v8["credible_set_variant_ids"].values() for s in sets}
    if len(target)!=4 or not target.issubset(pips):raise ValueError("Expected exactly four valid AIS full CS variants")
    plan={"status":"PLAN_ONLY","genes":list(SYMBOLS),"gwas_grch37_snps":len(gw),
       "credible_set_snps":len(target),"source":"GTEx_Portal_v2_GTex_v8_DAPG_GRCh38",
       "assembly_conversion":"hg38ToHg19.over.chain_1to1","scope":"QTL_FINEMAP_VARIANT_OVERLAP_ONLY"}
    if not args.execute:
        print(json.dumps(plan,indent=2));return plan
    if LiftOver is None:
        raise RuntimeError("pyliftover required for GTEx execution; install in isolated reference/vendor_pyliftover")
    import requests
    converter=LiftOver(str(args.chain))
    args.out.mkdir(parents=True,exist_ok=True)
    cache=args.out/"raw_gtex_v8";cache.mkdir(exist_ok=True)
    session=requests.Session()
    session.headers.update({"User-Agent":"IS-Master-Degree-research/1.0 (public GTEx API research metadata)"})
    allrows=[];qc=Counter();by_gene={}
    for symbol in SYMBOLS:
        gencode,eqtl=fetch_gene(session,symbol,cache)
        if stable(gencode)!=stable(src[symbol]["gene_id"]):
            raise ValueError("Ensembl stable gene mismatch between v19/v26: "+symbol)
        by_gene[symbol]={"gtEx_v8_gencode_id":gencode,"gtex_dapg_rows":len(eqtl)}
        for e in eqtl:
            if e.get("gencodeId")!=gencode or e.get("datasetId")!="gtex_v8":
                raise ValueError("Source row gene/dataset misidentified")
            variant=e.get("variantId")
            mapped,state=variant_liftover(variant,converter)
            qc[state]+=1
            if mapped is None:continue
            snpid,strand=mapped
            gw_hit=pips.get(snpid)
            # A reverse REF/ALT pair is a separate orientation conflict. Do not
            # assume a swapped pair represents a valid GTEx effect allele.
            if gw_hit is None:
                pos=snpid.split(":")[1]
                bypos=[k for k in pips if k.startswith("12:"+pos+":")]
                category="POSITION_PRESENT_ALLELES_INCOMPATIBLE" if bypos else "NO_GWAS_GENOTYPE_SNP_MATCH"
                qc[category]+=1
                continue
            pp=float(e["pip"])
            if not 0<=pp<=1:raise ValueError("Invalid GTEx PIP")
            tissue=e.get("tissueSiteDetailId") or ""
            if not tissue:raise ValueError("Missing GTEx tissue context")
            allrows.append({
                "gtex_gene_symbol":symbol,"gtex_gencode_v26":gencode,
                "gencode_v19":src[symbol]["gene_id"],
                "gtex_variant_grch38":variant,
                "variant_grch37_ref_alt":snpid,
                "liftover_strand":strand,"tissue":tissue,
                "gtex_method":e.get("method") or "",
                "gtex_finemap_set_id":e.get("setId"),
                "gtex_set_size":e.get("setSize"),
                "gtex_finemap_pip":pp,
                "gwas_ais_full_locus_pip":float(gw_hit.get("z") or 0), # corrected later
                "gwas_is_95pct_cs":int(snpid in target),
                "source_dataset":"GTEX_V8_DAPG",
                "interpretation":"FINEMAP_VARIANT_OVERLAP_NOT_COLOCALIZATION"})
        print(json.dumps({"gene":symbol,"dapg_rows":len(eqtl),
           "mapped_variant_overlaps_so_far":len(allrows)}),flush=True)
    # Replace placeholder with full-sample SuSiE output PIP, not GWAS Z score.
    assoc={x["variant_id"]:x for x in read_tsv(args.root/"G0022_FULL_AIS_PIP.tsv")}
    for row in allrows:
        row["gwas_ais_full_locus_pip"]=float(assoc[row["variant_grch37_ref_alt"]]["pip"])
    dest=args.out/"G0022_GTEX_V8_DAPG_GWAS_VARIANT_OVERLAPS.tsv"
    header=["gtex_gene_symbol","gtex_gencode_v26","gencode_v19","gtex_variant_grch38",
       "variant_grch37_ref_alt","liftover_strand","tissue","gtex_method",
       "gtex_finemap_set_id","gtex_set_size","gtex_finemap_pip",
       "gwas_ais_full_locus_pip","gwas_is_95pct_cs","source_dataset","interpretation"]
    with dest.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=header,delimiter="\t")
        w.writeheader();w.writerows(allrows)
    per_gene=[]
    for gene in SYMBOLS:
        r=[x for x in allrows if x["gtex_gene_symbol"]==gene]
        cs=[x for x in r if x["gwas_is_95pct_cs"]]
        per_gene.append({
            "gene":gene,"prior_annotation":"GENCODE_V19_G0022_POSITIONAL",
            "gtex_v8_dapg_rows":by_gene[gene]["gtex_dapg_rows"],
            "gwas_genotype_exact_allele_overlap_records":len(r),
            "gwas_cs_exact_allele_overlap_records":len(cs),
            "gwas_cs_distinct_overlap_snps":len({x["variant_grch37_ref_alt"] for x in cs}),
            "gtex_tissues_with_cs_overlap":";".join(sorted({x["tissue"] for x in cs})),
            "maximum_gtex_pip_among_cs":max((x["gtex_finemap_pip"] for x in cs),default=""),
            "qtl_colocalization_status":"NOT_TESTED_FULL_SUMSTATS_NOT_YET_IMPORTED"})
    with (args.out/"G0022_GTEX_V8_DAPG_PER_GENE_AUDIT.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(per_gene[0]),delimiter="\t")
        w.writeheader();w.writerows(per_gene)
    report={**plan,"status":"GTEX_V8_DAPG_FINEMAP_OVERLAP_AUDITED",
       "gene_assayed":len(SYMBOLS),"total_gtex_dapg_rows":sum(x["gtex_dapg_rows"] for x in by_gene.values()),
       "gwas_allele_exact_finemap_overlap_records":len(allrows),
       "gwas_credible_set_overlap_records":sum(x["gwas_is_95pct_cs"] for x in allrows),
       "gwas_credible_set_unique_qtl_overlap_snps":len({x["variant_grch37_ref_alt"] for x in allrows if x["gwas_is_95pct_cs"]}),
       "liftover_and_variant_qc":dict(qc),"per_gene":by_gene,
       "prior_BBJ_coloc_not_transferred":True,"QTL_colocalization_done":False,
       "per_variant_effective_n_verified":False,
       "interpretation":"GTEx DAP-G QTL fine-map SNP overlap only; cannot infer shared causal signal or independent replication"}
    (args.out/"G0022_GTEX_V8_DAPG_OVERLAP_SUMMARY.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return report
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--genes",type=Path,default=GENES)
    p.add_argument("--chain",type=Path,default=CHAIN)
    p.add_argument("--out",type=Path,default=OUT)
    p.add_argument("--execute",action="store_true")
    args=p.parse_args()
    fetch(args)
