#!/usr/bin/env python3
"""Frozen CKM common-instrument recovery; Skin official credible-set inventory (not MR)."""
import argparse,collections,hashlib,json,math,re,sys,tempfile,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
for folder in ("prior","ckm","indexed"):
    for sub in ("scripts","scripts/is","scripts/publication"):
        sys.path.insert(0,str(ROOT/folder/sub))
import early_track_acceleration as base
import publication_common_instrument_pilot as old
import fixed_subtype_indexed as idx
base.HOSTS.add("api.platform.opentargets.org")
REQUIRED_IDS=("rs646776","rs1367117","rs9875338","rs6818397","rs12916","rs1564348","rs2073547","rs2980875","rs579459","rs2419604","rs174583","rs1169288","rs4942486","rs8017377","rs2000999","rs6504872","rs7412","rs6065311","rs5763662")
OUTCOMES={"LAA":("EUR_LAS","Large artery stroke"),"CE":("EUR_CES","Cardioembolic stroke"),"SVD":("EUR_SVS","Small vessel stroke")}
SKIN=("GCST90094903","GCST90094904")
require=base.require
def mapping(obj,r,assembly):
    if not isinstance(obj,dict):raise ValueError("missing exact variant object")
    require(obj.get("name")==r["rsid"],"variant identity differs")
    q=[m for m in obj.get("mappings",[]) if m.get("assembly_name")==assembly and str(m.get("seq_region_name"))==str(r["chrom"]) and m.get("start")==m.get("end") and m.get("strand")==1]
    require(len(q)==1,"missing/nonunique primary coordinate")
    x=q[0];a=x["allele_string"].split("/")
    require(len(a)==2 and len(set(a))==2 and all(re.fullmatch("[ACGT]",v) for v in a),"not exact biallelic SNV")
    require(set(a) not in (set("AT"),set("CG")),"palindromic reference")
    require(old.orientation(r["ea"],r["oa"],*a) is not None,"exposure allele mismatch")
    return int(x["start"]),a
def resolve(r,a,b):
    p37,v37=mapping(a,r,"GRCh37");p38,v38=mapping(b,r,"GRCh38")
    require(p37==r["pos"],"source GRCh37 coordinate differs")
    require(set(v37)==set(v38),"allele set differs across assemblies")
    ref,alt=v38
    sign=old.orientation(r["ea"],r["oa"],alt,ref)[0]
    return {"status":"MATCHED_REFERENCE","chr":str(r["chrom"]),"pos":p38,"ref":ref,"alt":alt,"exposure_sign":sign,
            "grch37_position":p37,"grch37_alleles":v37,"reference_order_changed":v37!=v38}
def reference():
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/"ldl.gz";src=old.download(old.LDL_URL,p,60_000_000,180)
        require(src["sha256"]=="ddf75ed591ffcce8c3e81ae2f4c77f4428cd7dbf4628b24467aa8cf6e85e3dce","frozen exposure hash")
        panel,selection=old.load_ldl(p)
    require(tuple(r["rsid"] for r in panel)==REQUIRED_IDS and selection["source_rows"]==2437751 and selection["eligible_genomewide_significant_rows"]==2511,"frozen panel drift")
    refs={}
    for assembly,host in (("GRCh37","grch37.rest.ensembl.org"),("GRCh38","rest.ensembl.org")):
        raw,prov=base.get("https://"+host+"/variation/homo_sapiens",10_000_000,body=json.dumps({"ids":list(REQUIRED_IDS)}).encode())
        refs[assembly]=(json.loads(raw),prov)
    res={}
    for r in panel:
        k=r["rsid"]
        try:res[k]=resolve(r,refs["GRCh37"][0].get(k),refs["GRCh38"][0].get(k))
        except (ValueError,KeyError,TypeError,IndexError) as e:res[k]={"status":"HOLD_REFERENCE","reason":str(e)[:150]}
    return {"status":"FROZEN_LDL_REFERENCE_RECOVERY","panel":panel,"selection":selection,"exposure_source":src,"reference_sources":{k:v[1] for k,v in refs.items()},"reference":res,
            "matched":sum(v["status"]=="MATCHED_REFERENCE" for v in res.values()),
            "limits":["Same 19 frozen exposure SNPs; no replacements, proxy or multiallelic rescue. GRCh37 exact source position and both unique primary assemblies required. Reference allele ordering may legitimately change between assemblies; set and effect alignment must still agree. No LD certification."]}
def read_reference(root):
    x=json.loads((root/"reference.json").read_text())
    require(x.get("status")=="FROZEN_LDL_REFERENCE_RECOVERY" and tuple(v["rsid"] for v in x["panel"])==REQUIRED_IDS and set(x["reference"])==set(REQUIRED_IDS),"reference artifact contract")
    return x
def extract(name,root):
    import urllib.error
    x=read_reference(root);refs=x["reference"];idx.PANEL=REQUIRED_IDS
    get=idx.get
    def retry(url,cap,span=None):
        for attempt in range(2):
            try:return get(url,cap,span)
            except urllib.error.HTTPError as e:
                if e.code not in (429,500,502,503,504) or attempt:raise
            except OSError:
                if attempt:raise
            time.sleep(4)
    idx.get=retry
    study,trait=OUTCOMES[name];z=idx.extract(study,refs);md=z.get("metadata",{})
    verified=(z.get("accession")==idx.STUDIES[study] and md.get("trait_description")==[trait] and {v for row in md.get("samples",[]) for v in row.get("sample_ancestry_category",[])}=={"European"} and bool(z.get("fields",{}).get("rsid")))
    matched={}
    if verified:
        for rs,row in z["rows"].items():
            if row["status"]=="MATCHED_HARMONISED_ROW":
                require(refs[rs]["status"]=="MATCHED_REFERENCE","reference gate lost")
                matched[rs]={"by":row["beta_alt"]*refs[rs]["exposure_sign"],"sy":row["se"]}
    return {"status":"CKM_INDEPENDENT_INDEXED_OUTCOME","outcome":name,"reference_exposure_sha":x["exposure_source"]["sha256"],"source_metadata_qualified":verified,"matched":len(matched),"effects":matched,"readback":z,
            "limits":["Exact rsID, coordinate, resolved nonpalindromic allele, metadata and bounded BGZF ranges only. Not original whole-file checksum. Network errors are HOLD, not absent associations."]}
def combine(root):
    x=read_reference(root);ans={};coverage={}
    for k in OUTCOMES:
        p=root/(k+".json")
        v=json.loads(p.read_text()) if p.exists() else {}
        require(not v.get("effects") or v.get("outcome")==k and v.get("reference_exposure_sha")==x["exposure_source"]["sha256"],"outcome identity")
        ans[k]=v.get("effects",{}) if v.get("source_metadata_qualified") is True else {}
        require(set(ans[k])<=set(REQUIRED_IDS),"unexpected outcome instrument")
        coverage[k]={"matched":len(ans[k]),"status":v.get("status","HOLD_MISSING_ARTIFACT")}
    return {"status":"CKM_FROZEN_COMMON_PANEL_READOUT","reference_matched":x["matched"],"coverage":coverage,"analysis":old.analyze(x["panel"],ans),
            "limits":["Original #334 minimum common instrument count, fixed six contrast family, worst-case per-SNP cross-outcome covariance, leave-one-out. No claims of causal or ancestry-specific mechanism; exposure standardization, sample overlap and between-SNP LD still unresolved. Missing not zero."]}
def gql(query):
    body=json.dumps({"query":query}).encode()
    raw,src=base.get("https://api.platform.opentargets.org/api/v4/graphql",4_000_000,body=body)
    obj=json.loads(raw);require(not obj.get("errors"),"GraphQL schema/response error: "+str(obj.get("errors",""))[:180])
    require(isinstance(obj.get("data"),dict),"missing GraphQL data")
    return obj["data"],src
def skin(acc):
    require(acc in SKIN,"unexpected skin accession")
    folder=base.SKIN_BASE+acc+"/";catalog={}
    for label,url in (("original_directory",folder),("harmonised_directory",folder+"harmonised/")):
        try:
            raw,src=base.get(url,2_000_000)
            files=[v for v in base.links(raw) if re.fullmatch(r"[A-Za-z0-9_.-]+\.(?:tsv(?:\.gz)?|txt|yaml|tbi|md5|log)",v)]
            catalog[label]={"source":src,"files":files}
        except (OSError,ValueError) as e:catalog[label]={"status":"HOLD_DIRECTORY","reason":str(e)[:150]}
    # Open Targets API is an alternative processed derivative of the SAME Korean study.
    # No whole sumstats or MR proxy; retrieve all credible sets within a strict page bound.
    md,ms=gql('query { study(studyId:"'+acc+'") { id traitFromSource pubmedId publicationFirstAuthor nSamples hasSumstats summarystatsLocation qualityControls } }')
    s=md.get("study");require(isinstance(s,dict) and s.get("id")==acc and str(s.get("pubmedId"))=="35421404" and "wrinkle" in str(s.get("traitFromSource","")).lower(),"Korean wrinkle study identity")
    q='query { credibleSets(studyIds:["'+acc+'"],page:{index:0,size:100}) { count rows { studyLocusId studyId chromosome position beta standardError pValueMantissa pValueExponent finemappingMethod confidence purityMinR2 qualityControls variant { id } } } }'
    data,ps=gql(q);p=data["credibleSets"];rows=p["rows"]
    require(isinstance(p["count"],int) and p["count"]==len(rows) and len(rows)<=100,"incomplete credible-set page")
    require(len({x["studyLocusId"] for x in rows})==len(rows) and all(x["studyId"]==acc for x in rows),"study locus identity")
    for z in rows:
        v=z.get("variant") or {};require(re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT)_\d+_[ACGT]+_[ACGT]+",v.get("id","")) is not None,"unresolved lead variant")
        require(str(v["id"].split("_")[0])==str(z["chromosome"]) and int(v["id"].split("_")[1])==int(z["position"]),"lead coordinate disagreement")
        for k in ("beta","standardError","pValueMantissa","purityMinR2"):
            require(z.get(k) is None or isinstance(z[k],(int,float)) and math.isfinite(z[k]),"invalid lead statistic")
    return {"status":"SKIN_DERIVED_CREDIBLE_SET_INVENTORY","accession":acc,"catalog":catalog,"study":s,"study_source":ms,"credible_set_source":ps,"count":p["count"],"credible_sets":rows,
            "limits":["Open Targets postprocessed fine-mapping of the same original Korean GWAS (Lee et al 2022), not independent replication or full original summary statistics. Lead effects may be missing; no new colocalization, instrument selection, MR or merging distinct wrinkle phenotypes. Check ancestry and LD source before mechanistic use."]}
def selftest():
    r={"rsid":"rs123","chrom":1,"pos":10,"ea":"A","oa":"C"}
    def obj(build,pos,alleles):
        return {"name":"rs123","mappings":[{"assembly_name":build,"seq_region_name":"1","start":pos,"end":pos,"strand":1,"allele_string":alleles}]}
    z=resolve(r,obj("GRCh37",10,"A/C"),obj("GRCh38",20,"C/A"))
    assert z["reference_order_changed"] and z["exposure_sign"]==1 and z["pos"]==20
    for a,b in [(obj("GRCh37",11,"A/C"),obj("GRCh38",20,"C/A")),(obj("GRCh37",10,"A/C"),obj("GRCh38",20,"C/A/G"))]:
        try:resolve(r,a,b)
        except ValueError:pass
        else:raise AssertionError("invalid mapping accepted")
    assert len(REQUIRED_IDS)==19 and len(set(REQUIRED_IDS))==19 and len(OUTCOMES)==3 and len(SKIN)==2
    print("Reference-order flip, exact coordinate, multiallelic rejection, and fixed panel checks passed.")
def main():
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["reference","combine","selftest",*OUTCOMES,*SKIN]);p.add_argument("--root",type=Path,default=Path("inputs"));p.add_argument("--out",type=Path);a=p.parse_args()
    if a.mode=="selftest":selftest();return
    try:
        if a.mode=="reference":x=reference()
        elif a.mode=="combine":x=combine(a.root)
        elif a.mode in OUTCOMES:x=extract(a.mode,a.root)
        else:x=skin(a.mode)
    except (OSError,ValueError,TypeError,KeyError,IndexError) as e:x={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:220]}
    x["mode"]=a.mode;x["scientific_status"]="HOLD_INDEPENDENT_VALIDATION"
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
