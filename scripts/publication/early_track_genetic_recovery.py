#!/usr/bin/env python3
"""Frozen CKM common-instrument recovery; Skin official credible-set inventory (not MR)."""
import argparse,collections,concurrent.futures,csv,gzip,hashlib,json,math,re,subprocess,sys,tempfile,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
for folder in ("prior","ckm","indexed"):
    for sub in ("scripts","scripts/is","scripts/publication"):
        sys.path.insert(0,str(ROOT/folder/sub))
import early_track_acceleration as base
import publication_common_instrument_pilot as old
import fixed_subtype_indexed as idx
base.HOSTS.update({"api.platform.opentargets.org","api.ncbi.nlm.nih.gov"})
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
def ref_query(host,ids):
    from urllib.request import Request,urlopen
    from urllib.parse import urlsplit
    require(host in {"grch37.rest.ensembl.org","rest.ensembl.org"},"reference host")
    if len(ids)>1:
        url="https://"+host+"/variation/homo_sapiens"
        data=json.dumps({"ids":ids}).encode()
    else:
        url="https://"+host+"/variation/human/"+ids[0]+"?content-type=application/json"
        data=None
    req=Request(url,data=data,headers={"User-Agent":"IS-frozen-reference-recovery/1.1","Accept":"application/json","Content-Type":"application/json","Accept-Encoding":"identity"})
    last=None
    for attempt in range(2):
        try:
            with urlopen(req,timeout=18) as f:
                require(f.status==200 and urlsplit(f.geturl()).hostname==host,"reference HTTP/redirect")
                raw=f.read(2_000_001)
            require(0<len(raw)<=2_000_000,"reference transfer cap")
            obj=json.loads(raw)
            if len(ids)==1:obj={ids[0]:obj}
            require(isinstance(obj,dict),"reference JSON object")
            return obj,{"host":host,"ids":ids,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()}
        except (OSError,ValueError) as e:
            last=e
            if attempt==0:time.sleep(2)
    raise last
def reference_build(host,ids):
    obs={};provenance=[];failures={}
    try:
        obj,src=ref_query(host,ids);provenance.append(src)
        obs.update({k:v for k,v in obj.items() if k in ids and isinstance(v,dict) and v.get("name")==k})
    except (OSError,ValueError) as e:failures["batch"]=str(e)[:150]
    missing=[k for k in ids if k not in obs]
    def one(k):
        try:
            o,src=ref_query(host,[k]);v=o.get(k)
            require(isinstance(v,dict) and v.get("name")==k,"individual reference ID")
            return k,v,src,None
        except (OSError,ValueError,TypeError) as e:return k,None,None,str(e)[:150]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for k,v,src,err in pool.map(one,missing):
            if v is not None:obs[k]=v;provenance.append(src)
            else:failures[k]=err
    return obs,{"host":host,"received":len(obs),"batch_or_individual_sources":provenance,"transport_failures":failures}
def refsnp_placement(obj,r,build):
    # RefSNP SPDI position is 0-based; input and output genomic positions are 1-based.
    require(str(obj.get("refsnp_id"))==r["rsid"][2:],"exact dbSNP identity")
    snap=obj.get("primary_snapshot_data",{})
    require(snap.get("variant_type")=="snv","no primary SNV snapshot")
    chrom_prefix="NC_%06d."%int(r["chrom"])
    found=[]
    for p in snap.get("placements_with_allele",[]):
        ann=p.get("placement_annot",{});seq=p.get("seq_id","")
        if not seq.startswith(chrom_prefix) or ann.get("seq_type")!="refseq_chromosome" or ann.get("mol_type")!="genomic":continue
        traits=[a for a in ann.get("seq_id_traits_by_assembly",[]) if
            re.fullmatch(build+r"(?:\.p\d+)?",a.get("assembly_name","")) and
            re.fullmatch(r"GCF_000001405\.\d+",a.get("assembly_accession","")) and
            a.get("is_top_level") is True and a.get("is_chromosome") is True and
            a.get("is_alt") is False and a.get("is_patch") is False]
        if not traits:continue
        alleles=[a.get("allele",{}).get("spdi",{}) for a in p.get("alleles",[])]
        require(len(alleles)==2,"multiallelic or incomplete placement")
        require(all(v.get("seq_id")==seq and isinstance(v.get("position"),int) and
            re.fullmatch("[ACGT]",v.get("deleted_sequence","")) and
            re.fullmatch("[ACGT]",v.get("inserted_sequence","")) for v in alleles),"noncanonical SPDI allele")
        pos={v["position"] for v in alleles};ref={v["deleted_sequence"] for v in alleles};alt={v["inserted_sequence"] for v in alleles}
        require(len(pos)==len(ref)==1 and len(alt)==2 and ref<=alt,"SPDI site/reference mismatch")
        p0=next(iter(pos));a=next(iter(ref));b=next(iter(alt-ref))
        require(p0>=0 and {a,b} not in ({"A","T"},{"C","G"}),"palindromic or invalid reference")
        found.append((p0+1,[a,b],seq,sorted({t["assembly_name"] for t in traits})))
    require(len(found)==1,"no unique primary chromosome placement "+build)
    return found[0]
def resolve_dbsnp(obj,r):
    p37,a37,s37,b37=refsnp_placement(obj,r,"GRCh37")
    p38,a38,s38,b38=refsnp_placement(obj,r,"GRCh38")
    require(p37==r["pos"],"source GRCh37 coordinate differs")
    require(set(a37)==set(a38),"two-build allele sets differ")
    ref,alt=a38;orient=old.orientation(r["ea"],r["oa"],alt,ref)
    require(orient is not None,"exposure alleles differ")
    return {"status":"MATCHED_REFERENCE","chr":str(r["chrom"]),"pos":p38,"ref":ref,"alt":alt,"exposure_sign":orient[0],
            "grch37_position":p37,"grch37_alleles":a37,"grch37_refseq":s37,"grch38_refseq":s38,"assembly_labels":[b37,b38],
            "reference_order_changed":a37!=a38}
def dbsnp_one(r):
    rs=r["rsid"];url="https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/"+rs[2:]
    error=None;source=None
    for attempt in range(2):
        try:
            raw,source=base.get(url,16_000_000)
            break
        except OSError as e:
            error=str(e)[:150]
            if attempt==0:time.sleep(2)
    if source is None:return rs,{"status":"HOLD_REFERENCE","reason":"dbSNP transport: "+str(error)},None
    try:result=resolve_dbsnp(json.loads(raw),r)
    except (ValueError,KeyError,TypeError,IndexError) as e:result={"status":"HOLD_REFERENCE","reason":str(e)[:150]}
    return rs,result,source

def reference():
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/"ldl.gz";exposure_source=old.download(old.LDL_URL,p,60_000_000,180,retry=1)
        require(exposure_source["sha256"]=="ddf75ed591ffcce8c3e81ae2f4c77f4428cd7dbf4628b24467aa8cf6e85e3dce","frozen exposure hash")
        panel,selection=old.load_ldl(p)
    require(tuple(r["rsid"] for r in panel)==REQUIRED_IDS and selection["source_rows"]==2437751 and selection["eligible_genomewide_significant_rows"]==2511,"frozen panel drift")
    res={};sources={}
    # Independent official source, not repeated requests to the failing Ensembl GRCh37 server.
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for rs,z,src in pool.map(dbsnp_one,panel):
            res[rs]=z
            if src:sources[rs]=src
    return {"status":"FROZEN_LDL_REFERENCE_RECOVERY","panel":panel,"selection":selection,"exposure_source":exposure_source,"reference_sources":{"dbSNP_RefSNP":sources},"reference":res,
            "matched":sum(v["status"]=="MATCHED_REFERENCE" for v in res.values()),
            "limits":["Same frozen 19 SNPs. Independent official RefSNP exact primary chromosome SPDI on both builds; original GRCh37 position and nonpalindromic biallelic allele set must match. No proxy, merged IDs, multiallelic rescue or LD certification. If <10 reference matches, stop outcome queries."]}

def read_reference(root):
    x=json.loads((root/"reference.json").read_text())
    require(x.get("status")=="FROZEN_LDL_REFERENCE_RECOVERY" and tuple(v["rsid"] for v in x["panel"])==REQUIRED_IDS and set(x["reference"])==set(REQUIRED_IDS),"reference artifact contract")
    return x
def extract(name,root):
    import urllib.error
    x=read_reference(root);refs=x["reference"];idx.PANEL=REQUIRED_IDS
    if x["matched"]<old.MIN_N:
        return {"status":"HOLD_REFERENCE_BELOW_COMMON_PANEL_MINIMUM","outcome":name,"reference_exposure_sha":x["exposure_source"]["sha256"],"reference_matched":x["matched"],"matched":0,"effects":{},"source_metadata_qualified":False,"limits":["Fewer than the fixed minimum of ten eligible references; no heavy outcome requests."]}
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
def skin_source_contract(acc,catalog,leads):
    """Older GWAS files are not indexed; exact complete harmonised file only, no reconstructed SE."""
    import yaml
    from urllib.request import Request,urlopen
    from urllib.parse import urlsplit
    folder=base.SKIN_BASE+acc+"/harmonised/"
    names=catalog.get("harmonised_directory",{}).get("files",[])
    matches=[n for n in names if re.fullmatch(r"35421404-"+acc+r"-EFO_\d+\.h\.tsv\.gz",n)]
    out={"status":"HOLD_FULL_SOURCE","lead_total":len(leads),"matched":0}
    stage="filename"
    try:
        require(len(matches)==1,"unique publication-prefixed harmonised file absent")
        name=matches[0];url=folder+name;out["filename"]=name
        stage="metadata"
        md_raw,ms=base.get(url+"-meta.yaml",500_000);md=yaml.safe_load(md_raw)
        require(isinstance(md,dict),"metadata object")
        out.update(metadata_source=ms,metadata={k:md.get(k) for k in ("gwas_id","data_file_name","genome_assembly","is_harmonised","trait_description","samples","data_file_md5sum")})
        require(md.get("gwas_id")==acc and md.get("data_file_name")==name,"metadata accession/filename mismatch")
        require(str(md.get("genome_assembly","")).upper()=="GRCH38" and md.get("is_harmonised") is True,"not verified harmonised GRCh38")
        stage="checksum"
        checks,cs=base.get(folder+"md5sum.txt",100_000)
        found=[]
        for line in checks.decode("ascii").splitlines():
            m=re.fullmatch(r"([a-fA-F0-9]{32})\s+\*?(?:\./)?(.+)",line.strip())
            if m and m[2]==name:found.append(m[1].lower())
        require(len(found)<=1,"duplicate exact named manifest MD5")
        md5_from_yaml=str(md.get("data_file_md5sum","")).lower()
        require(re.fullmatch("[a-f0-9]{32}",md5_from_yaml),"canonical filename-bound YAML MD5 absent")
        if found:
            require(found[0]==md5_from_yaml,"published checksums differ")
            checksum_basis="exact_manifest_and_exact_yaml"
        else:
            # YAML itself is official and names this exact accession/file. No inferred checksum
            # from a neighboring phenotype or another file. Retain manifest digest and absence.
            checksum_basis="exact_filename_bound_yaml_only_manifest_entry_absent"
        expected=md5_from_yaml
        out.update(expected_md5=expected,checksum_source=cs,checksum_basis=checksum_basis,manifest_exact_matches=len(found))
        stage="size"
        _,probe=base.get(url,1,(0,0))
        m=re.fullmatch(r"bytes 0-0/(\d+)",probe.get("content_range",""))
        require(m is not None,"missing exact remote size")
        size=int(m[1]);out.update(remote_size=size,size_probe=probe)
        require(0<size<=600_000_000,"compressed source over 600 MB cap")
        stage="complete_download"
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"source.gz";digest=hashlib.md5();sha=hashlib.sha256();n=0
            req=Request(url,headers={"User-Agent":"IS-fixed-skin-source/1.0","Accept-Encoding":"identity"})
            with urlopen(req,timeout=60) as f,p.open("wb") as dest:
                require(f.status==200 and urlsplit(f.geturl()).hostname=="ftp.ebi.ac.uk","source HTTP/redirect")
                for block in iter(lambda:f.read(1_048_576),b""):
                    n+=len(block);require(n<=size,"transfer exceeded exact size")
                    digest.update(block);sha.update(block);dest.write(block)
            require(n==size and digest.hexdigest()==expected,"complete source size/MD5 mismatch")
            out.update(full_source={"url":url,"bytes":n,"md5":digest.hexdigest(),"sha256":sha.hexdigest()})
            stage="header"
            targets=collections.defaultdict(list)
            for z in leads:
                v=z["variant"]["id"];ch,pos,ref,alt=v.split("_")
                targets[(ch,pos)].append((v,ref,alt))
            pool=collections.defaultdict(list);count=0;expanded=0;invalid=0
            with gzip.open(p,"rt",encoding="utf-8",newline="") as f:
                h=f.readline().rstrip("\r\n").lstrip("#").split("\t")
                require(len(set(h))==len(h),"duplicate GWAS header")
                out["header"]=h
                def choose(*keys):
                    for k in keys:
                        if k in h:return h.index(k)
                    raise ValueError("missing required field: "+"/".join(keys))
                ci=choose("hm_chrom","chromosome");pi=choose("hm_pos","base_pair_location")
                ea=choose("hm_effect_allele","effect_allele");oa=choose("hm_other_allele","other_allele")
                bi=choose("hm_beta","beta");si=choose("standard_error","hm_standard_error");pv=choose("p_value","hm_p_value")
                co=choose("hm_code","harmonisation_code")
                rs=next((h.index(k) for k in ("hm_rsid","rsid","rs_id") if k in h),None)
                stage="full_scan"
                maxcol=max(ci,pi)
                for line in f:
                    count+=1;expanded+=len(line)
                    require(expanded<=3_000_000_000 and count<=30_000_000,"decompressed scan cap")
                    a=line.split("\t",maxcol+1)
                    require(len(a)>maxcol,"short coordinate row")
                    hit=targets.get((a[ci].removeprefix("chr"),a[pi].strip()))
                    if hit is None:continue
                    row=line.rstrip("\r\n").split("\t");require(len(row)==len(h),"target row width")
                    for vid,ref,alt in hit:
                        x,y=row[ea].upper(),row[oa].upper()
                        if {x,y}!={ref,alt}:continue
                        z={"status":"HOLD_EFFECT","reason":"unresolved allele or statistic"}
                        try:
                            require(len(ref)==len(alt)==1 and {ref,alt} not in ({"A","T"},{"C","G"}),"palindromic or non-SNV lead")
                            require(int(row[co]) in (10,11,12,13),"unresolved harmonisation code")
                            b,se,pv0=float(row[bi]),float(row[si]),float(row[pv])
                            require(all(math.isfinite(v) for v in (b,se,pv0)) and se>0 and 0<=pv0<=1,"invalid beta/SE/p")
                            z={"status":"MATCHED_FULL_SOURCE_LEAD","beta_alt":b if x==alt else -b,"standard_error":se,"p_value":pv0,"hm_code":int(row[co])}
                            if rs is not None and re.fullmatch(r"rs\d+",row[rs]):z["rsid"]=row[rs]
                        except (ValueError,IndexError,TypeError) as e:z["reason"]=str(e)[:120]
                        pool[vid].append(z)
            require(count>=100_000,"incomplete GWAS row count")
            rows={}
            for z in leads:
                vid=z["variant"]["id"];v=pool[vid]
                rows[vid]=v[0] if len(v)==1 else {"status":"HOLD_MISSING_OR_NONUNIQUE_SOURCE_ROW","rows":len(v)}
            out.update(status="COMPLETE_MD5_VERIFIED_LEAD_SCAN",rows_scanned=count,lead_rows=rows,matched=sum(v["status"]=="MATCHED_FULL_SOURCE_LEAD" for v in rows.values()))
    except (OSError,EOFError,ValueError,TypeError,KeyError,IndexError,UnicodeError) as e:out.update(stage=stage,reason=str(e)[:200])
    out["limits"]=["Full file MD5/CRC and exact GRCh38 nonpalindromic lead matches required. Source-derived SE only. Same GWAS as Open Targets, not independent validation or LD-pruned MR. No estimate of missing SE from rounded P-values."]
    return out

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
    return {"status":"SKIN_DERIVED_CREDIBLE_SET_INVENTORY","accession":acc,"catalog":catalog,"study":s,"study_source":ms,"credible_set_source":ps,"count":p["count"],"credible_sets":rows,"full_source_readout":skin_source_contract(acc,catalog,rows),
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
    def placement(build,pos,ref,alt):
        return {"seq_id":"NC_000001.11" if build=="GRCh38" else "NC_000001.10",
                "placement_annot":{"seq_type":"refseq_chromosome","mol_type":"genomic","seq_id_traits_by_assembly":[{"assembly_name":build+".p13","assembly_accession":"GCF_000001405.25" if build=="GRCh37" else "GCF_000001405.39","is_top_level":True,"is_chromosome":True,"is_alt":False,"is_patch":False}]},
                "alleles":[{"allele":{"spdi":{"seq_id":"NC_000001.11" if build=="GRCh38" else "NC_000001.10","position":pos-1,"deleted_sequence":ref,"inserted_sequence":v}}} for v in (ref,alt)]}
    ob={"refsnp_id":"123","primary_snapshot_data":{"variant_type":"snv","placements_with_allele":[placement("GRCh37",10,"A","C"),placement("GRCh38",20,"C","A")]}}
    q=resolve_dbsnp(ob,r)
    assert q["pos"]==20 and q["grch37_position"]==10 and q["exposure_sign"]==1
    ob["primary_snapshot_data"]["placements_with_allele"][0]["alleles"][1]["allele"]["spdi"]["position"]=10
    try:resolve_dbsnp(ob,r)
    except ValueError:pass
    else:raise AssertionError("different SPDI positions accepted")
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
    except (OSError,ValueError,TypeError,KeyError,IndexError,subprocess.CalledProcessError) as e:x={"status":"HOLD_SOURCE_OR_CONTRACT","reason":str(e)[:220]}
    x["mode"]=a.mode;x["scientific_status"]="HOLD_INDEPENDENT_VALIDATION"
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
