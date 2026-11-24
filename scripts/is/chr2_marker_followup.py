#!/usr/bin/env python3
"""Fixed chr2 marker follow-up, not a novel locus or gene-causality claim."""
import argparse,csv,hashlib,itertools,json,math,re,urllib.request,urllib.parse,zipfile
from pathlib import Path
from statistics import NormalDist

POS=2645097
KEY="2:2645097:C:T"
REGION="IS_XDATA_G0003"
PIN="541c31269241e63e1995baa4101e4d1c63029f5975045bcdea1df3b27d530ef5"
HOST37="https://grch37.rest.ensembl.org"
HOST38="https://rest.ensembl.org"

def need(ok,msg):
    if not ok:raise ValueError(msg)
def get(base,args=None):
    url=base+("?" + urllib.parse.urlencode(args) if args else "")
    req=urllib.request.Request(url,headers={"Accept":"application/json","Content-Type":"application/json","User-Agent":"IS-public-locus-research/1.0"})
    with urllib.request.urlopen(req,timeout=18) as f:raw=f.read(5000001)
    need(len(raw)<=5000000,"Response cap")
    return json.loads(raw),hashlib.sha256(raw).hexdigest(),url
def safe(fn,*args):
    try:return fn(*args)
    except (OSError,ValueError,KeyError,TypeError,OverflowError) as e:return {"status":"HOLD","reason":str(e)[:350]}
def write(out,r):
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
    print(json.dumps(r,indent=2,allow_nan=False),flush=True)
def reconstructed_se(beta,p):
    need(math.isfinite(beta) and beta!=0 and 0<p<1,"Invalid beta/p for normal approximation")
    return abs(beta)/(-NormalDist().inv_cdf(p/2))
def contrast(root,legacy,out):
    path=root/"prior.zip"
    with path.open("rb") as f:need(hashlib.file_digest(f,"sha256").hexdigest()==PIN,"Previous artifact digest")
    with zipfile.ZipFile(path) as z:
        need(z.testzip() is None,"CRC")
        names=[n for n in z.namelist() if n.rsplit("/",1)[-1]=="contrast.json"]
        need(len(names)==1,"Contrast JSON member")
        prior=json.loads(z.read(names[0]))
    rows=[r for r in prior["results"] if r["region_id"]==REGION and r["variant"]==KEY]
    need(len(rows)==1 and prior["family"]==24,"Original region family")
    a=rows[0]
    need(a["status"]=="CONDITIONAL_NUMERIC_CONTRAST" and a["EAS"]["alt"]==a["EUR"]["alt"]=="T","Original allele alignment")
    with legacy.open() as f:
        matches=[r for r in csv.DictReader(f,delimiter="\t") if r["region_id"]==REGION]
    labels={"BBJ_IS","GWS_EAS_AS","GWS_EAS_AIS","GWS_EAS_CES","GWS_EAS_LAS","GWS_EAS_SVS"}
    need(len(matches)==6 and {r["source_dataset"] for r in matches}==labels,"Frozen archive rows")
    d={}
    for r in matches:
        need(r["group_source_lead_variant"]==KEY and r["group_source_lead_allele_status"]=="SOURCE_VARIANT_UNORDERED_ALLELE_PAIR_MATCH","Legacy allele identity")
        beta,p,eaf=map(float,(r["group_source_lead_beta_ALT"],r["group_source_lead_p"],r["group_source_lead_ALT_AF"]))
        se=reconstructed_se(beta,p)
        d[r["source_dataset"]]={"beta_alt":beta,"source_p":p,"alt_eaf":eaf,"se_inferred_normal":se,"approx_ci95":[beta-1.959964*se,beta+1.959964*se],"regional_gws_count":int(r["n_genomewide_p_lt_5e8"])}
    need(math.isclose(d["GWS_EAS_AIS"]["beta_alt"],a["EAS"]["beta_alt"],abs_tol=1e-9),"EAS beta mismatch")
    c=[]
    for x,y in itertools.combinations(("LAS","CES","SVS"),2):
        r,s=d["GWS_EAS_"+x],d["GWS_EAS_"+y]
        delta=r["beta_alt"]-s["beta_alt"];s1,s2=r["se_inferred_normal"],s["se_inferred_normal"]
        p0=math.erfc(abs(delta)/math.hypot(s1,s2)/math.sqrt(2));pw=math.erfc(abs(delta)/(s1+s2)/math.sqrt(2))
        c.append(dict(pair=x+"-"+y,delta=delta,approx_p_zero_cov=p0,approx_p_maxvar=pw,approx_bonferroni_p_zero_cov=min(1,3*p0),approx_bonferroni_p_maxvar=min(1,3*pw)))
    write(out,{"status":"ARCHIVAL_SUBTYPE_NORMAL_APPROXIMATION_NOT_CONFIRMATORY","variant":KEY,"verified_ancestry_contrast":a,"legacy_source_rows":d,"subtype_comparisons":c,
        "limits":["SE inferred from rounded archived beta and source P, conditional on two-sided normal Wald statistics; not actual subtype SE.",
                  "Unknown shared controls and model differences; maximum variance is a covariance bound conditional on these inferred SEs.",
                  "Three pairwise tests exploratory after marker selection; does not address genome-wide or prior selection multiplicity.",
                  "BBJ contributes to GIGASTROKE EAS, not independent replication."]})
def identity():
    data,sha,u=get(HOST37+f"/overlap/region/human/2:{POS}-{POS}",{"feature":"variation"})
    need(isinstance(data,list),"Unexpected overlap schema")
    ids=sorted({x["id"] for x in data if x.get("assembly_name")=="GRCh37" and str(x.get("seq_region_name"))=="2" and x.get("start")==x.get("end")==POS and {"C","T"}<=set(x.get("alleles",[])) and re.fullmatch(r"rs\d+",str(x.get("id","")))})
    need(len(ids)==1,"No unique exact GRCh37 rsID: "+repr(ids))
    rs=ids[0]
    seq,seqsha,_=get(HOST37+f"/sequence/region/human/2:{POS}..{POS}:1")
    need(seq["seq"].upper()=="C","GRCh37 reference differs")
    d,sha2,u2=get(HOST38+"/variation/human/"+rs)
    maps={(int(m["start"]),m["allele_string"]) for m in d["mappings"] if m.get("assembly_name")=="GRCh38" and str(m.get("seq_region_name"))=="2" and m.get("start")==m.get("end") and m.get("strand")==1 and m.get("allele_string") in ("C/T","G/A")}
    need(len(maps)==1,"No unique compatible biallelic GRCh38 mapping")
    pos38,als=next(iter(maps));ref,alt=als.split("/")
    lift,lsha,lu=get(HOST38+f"/map/human/GRCh37/2:{POS}..{POS}:1/GRCh38")
    ms=lift["mappings"];need(len(ms)==1,"Ambiguous build lift")
    m=ms[0]["mapped"]
    need(str(m["seq_region_name"])=="2" and m["start"]==m["end"]==pos38 and m["strand"]==(1 if als=="C/T" else -1),"Build/allele discordance")
    seq38,seq38sha,_=get(HOST38+f"/sequence/region/human/2:{pos38}..{pos38}:1")
    need(seq38["seq"].upper()==ref,"GRCh38 reference differs")
    return dict(status="EXACT_BUILD_AND_RSID",rsid=rs,variant37=KEY,pos38=pos38,ref38=ref,alt38=alt,variant_id=f"chr2_{pos38}_{ref}_{alt}_b38",source_hashes=[sha,sha2,lsha,seqsha,seq38sha],urls=[u,u2,lu])
def koges():
    data,sha,u=get("https://koges.leelabsg.org/api/variant/2-2645097-C-T")
    need((str(data.get("chrom")),int(data.get("pos",0)),data.get("ref"),data.get("alt"))==("2",POS,"C","T"),"KoGES exact variant mismatch")
    rsids=sorted(set(re.findall(r"rs\d+",str(data.get("rsids",data.get("rsid",""))))))
    need(isinstance(data.get("phenos"),list),"KoGES phenotype schema")
    hits=[]
    for p in data["phenos"]:
        label=str(p.get("phenostring",""));code=str(p.get("phenocode",""))
        if not re.search(r"stroke|cerebral infarct|cerebrovascular|cerebral arter|brain infarct|ischemic|ischaemic|carotid",label+" "+code,re.I):continue
        row=dict(label=label,code=code,available_fields=sorted(p))
        for k in ("beta","sebeta","pval","num_cases","num_controls","af","maf"):
            if k in p:row[k]=p[k]
        hits.append(row)
    return {"status":"EXACT_COORDINATE_ALLELE_PHEWAS_LOOKUP","variant":KEY,"rsids":rsids,"n_all_phenotypes":len(data["phenos"]),"stroke_related_rows":hits,"sha256":sha,"url":u,
        "limits":"Only documented KoGES-only phenotype can be a supportive lookup; do not treat coronary disease, stroke without subtype or BBJ meta as independent AIS replication."}
def vep(v):
    data,sha,u=get(HOST38+"/vep/human/id/"+v["rsid"],{"canonical":1})
    hits=[r for r in data if r.get("assembly_name")=="GRCh38" and str(r.get("seq_region_name"))=="2" and r.get("start")==r.get("end")==v["pos38"] and r.get("allele_string","").split("/")[0]==v["ref38"] and v["alt38"] in r.get("allele_string","").split("/")[1:]]
    need(len(hits)==1,"VEP exact allele mapping")
    r=hits[0];tc=[q for q in r.get("transcript_consequences",[]) if q.get("variant_allele")==v["alt38"]]
    return dict(status="EXACT_ALT_VEP",all_alt_terms=sorted({x for t in tc for x in t.get("consequence_terms",[])}),transcripts=[{k:t[k] for k in ("gene_symbol","gene_id","transcript_id","canonical","consequence_terms","distance","hgvsp") if k in t} for t in tc],sha256=sha,url=u)
def gtex(v):
    allrows=[];sources=[]
    for page in range(5):
        data,sha,u=get("https://gtexportal.org/api/v2/association/singleTissueEqtl",{"datasetId":"gtex_v10","variantId":v["variant_id"],"itemsPerPage":1000,"page":page})
        need(isinstance(data.get("data"),list),"GTEx schema")
        rows=data["data"]
        for r in rows:
            need(r["variantId"]==v["variant_id"] and r.get("datasetId","gtex_v10")=="gtex_v10","GTEx identity")
            allrows.append({k:r[k] for k in ("gencodeId","geneSymbol","tissueSiteDetailId","nes","pValue","variantId") if k in r})
        sources.append({"url":u,"sha256":sha})
        if len(rows)<1000:break
    else:raise ValueError("GTEx bounded pagination")
    need(len({(r["gencodeId"],r["tissueSiteDetailId"]) for r in allrows})==len(allrows),"Duplicate GTEx rows")
    return dict(status="SIGNIFICANT_ONLY",n=len(allrows),rows=allrows,sources=sources)
def catalog(v):
    # Exact rsID curated associations, no negative novelty claim from absence.
    data,sha,u=get("https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/"+v["rsid"]+"/associations",{"size":100})
    rows=data.get("_embedded",{}).get("associations",[])
    need(isinstance(rows,list),"Catalog schema")
    return dict(status="CURATED_ASSOCIATIONS_FIRST_PAGE",total=data.get("page",{}).get("totalElements"),rows=[{k:r[k] for k in ("pvalue","pvalueMantissa","pvalueExponent","riskFrequency","orPerCopyNum","betaNum","betaUnit","range","description","efoTraits","loci","_links") if k in r} for r in rows],sha256=sha,url=u)
def ref(out):
    k=safe(koges);v=safe(identity)
    if v.get("rsid") and k.get("rsids") and v["rsid"] not in k["rsids"]:
        k={"status":"HOLD","reason":"KoGES rsID discordance","source":k}
    r={"status":"REFERENCE_AND_PUBLIC_VARIANT_CONTEXT_NOT_REPLICATION","reference":v,"koges":k}
    if v.get("status")=="EXACT_BUILD_AND_RSID":
        r.update(vep=safe(vep,v),gtex=safe(gtex,v),catalog=safe(catalog,v))
    r["limits"]=["A single selected variant lookup, no independent cohort overlap verification or phenotype harmonization yet.",
                 "Significant-only QTL and curated GWAS associations not full regional evidence; no causal gene, fine-mapping, or coloc.",
                 "Missing rows or failed queries are not null associations."]
    write(out,r)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["contrast","reference"]);p.add_argument("--root",type=Path,default=Path("."));p.add_argument("--legacy",type=Path);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    if a.mode=="contrast":contrast(a.root,a.legacy,a.out)
    else:ref(a.out)
