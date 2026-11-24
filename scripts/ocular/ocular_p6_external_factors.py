#!/usr/bin/env python3
"""P6: TM external-factor risk audit, exploratory Hallmark GSEA, public microbiome ENA feasibility.

P6a and P6b read preexisting public gene/library aggregates, never donor-level data.
P6c uses ENA read_run *metadata only*. No FASTQ, organism abundance, patient identifiers.
Nothing here tests bacterial presence or steroid-induced ocular hypertension clinically.
"""
import argparse,csv,gzip,hashlib,io,json,math,random,urllib.request,urllib.parse
from collections import Counter,defaultdict
from pathlib import Path

TARGETS={"ERFE","CHAC1","PARM1","ADAMTS15","EGR3","LPAR1","ABHD3",
         "DDIT4L","TFPI","STK17B","RAB32","SDC1"}
HALLMARK={
 "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
 "HALLMARK_INFLAMMATORY_RESPONSE",
 "HALLMARK_TGF_BETA_SIGNALING",
 "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
 "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
 "HALLMARK_IL6_JAK_STAT3_SIGNALING",
 "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
 "HALLMARK_HYPOXIA",
 "HALLMARK_APOPTOSIS"
}
GMT_URL="https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt"
STUDIES=[
 {"accession":"PRJNA1145861","label":"ocular_surface_16s_glaucoma_dry_eye",
  "platform_claim":"16S amplicon, conjunctival swabs"},
 {"accession":"PRJEB74248","label":"ocular_surface_shotgun_tear_proteome",
  "platform_claim":"shotgun conjunctival metagenomics; tear proteome is a separate modality"},
 {"accession":"PRJNA1416399","label":"oral_16s_open_angle_glaucoma",
  "platform_claim":"16S amplicon, mouthwash"}
]
USER_AGENT="ocular-p6-public-meta-audit/1.0 (research reproducibility)"
def get(url,limit=2_000_000,timeout=70):
    request=urllib.request.Request(url,headers={"User-Agent":USER_AGENT})
    with urllib.request.urlopen(request,timeout=timeout) as h: data=h.read(limit+1)
    if len(data)>limit: raise ValueError("REMOTE_RESPONSE_EXCEEDS_LIMIT")
    return data
def save(out,name,body):
    (out/name).write_text(json.dumps(body,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf8")
def tsv(out,name,rows):
    if not rows:raise ValueError("NO_OUTPUT_ROWS")
    with (out/name).open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
def num(x):
    v=float(x)
    if not math.isfinite(v):raise ValueError("NONFINITE_NUMBER")
    return v

def qc_source(inp,out):
    j=json.loads(inp.read_text(encoding="utf8"))
    if j.get("dataset")!="GSE313736" or j.get("status")!="GO_57DONOR_FULL_RAW_PSEUDOBULK_QC":
        raise ValueError("P5Q_SOURCE_PROVENANCE_DRIFT")
    rows=j["per_library"]
    if len(rows)!=18 or len({x["gsm"] for x in rows})!=18:
        raise ValueError("EXPECTED_18_UNIQUE_LIBRARY_IDS")
    arms=Counter((x["pool"],x["condition"]) for x in rows)
    if arms!={(pool,arm):1 for pool in "ABCDEFGHI" for arm in ("vehicle","dex")}:
        raise ValueError("EXPECTED_9x2_SOURCE_POOLS")
    data=[]
    for pool in "ABCDEFGHI":
        r={arm:next(x for x in rows if x["pool"]==pool and x["condition"]==arm)
           for arm in ("vehicle","dex")}
        for arm in r:
            if r[arm]["status"]!="PASS_PSEUDOBULK" or r[arm]["matrix_info"]["gzip_crc"]!="PASS":
                raise ValueError("RAW_SOURCE_MATRIX_INTEGRITY_DRIFT")
        a,b=r["vehicle"],r["dex"]
        va=num(a["selected_barcodes"]);db=num(b["selected_barcodes"])
        vu=num(a["matrix_info"]["selected_umis"]);du=num(b["matrix_info"]["selected_umis"])
        if min(va,db,vu,du)<=0:raise ValueError("NONPOSITIVE_SOURCE_COUNTS")
        data.append({"pool":pool,"vehicle_gsm":a["gsm"],"dex_gsm":b["gsm"],
          "vehicle_selected_barcodes":int(va),"dex_selected_barcodes":int(db),
          "dex_vehicle_barcode_ratio":round(db/va,5),
          "vehicle_selected_umis":int(vu),"dex_selected_umis":int(du),
          "dex_vehicle_umi_ratio":round(du/vu,5),
          "treatment_aliased_with_specific_library_id":True,
          "technical_batch_randomization_documented_here":False})
    tsv(out,"P6a_9_library_pairs_public_qc.tsv",data)
    result={"status":"PASS_INPUT_INTEGRITY_HOLD_MYCO_AND_TECHNICAL_BATCH",
      "dataset":"GSE313736","observed_libraries":18,"observed_source_pools":9,
      "same_library_contains_both_treatments":False,
      "treatment_condition_mapped_to_distinct_library_per_pool":True,
      "biological_donor_pairing_confirmed_in_prior_stage":True,
      "mycoplasma_PCR_assay_status":"NOT_REPORTED_IN_THIS_PUBLIC_AGGREGATE",
      "culture_bacterial_test_status":"NOT_REPORTED_IN_THIS_PUBLIC_AGGREGATE",
      "endotoxin_LPS_assay_status":"NOT_REPORTED_IN_THIS_PUBLIC_AGGREGATE",
      "library_prep_cross_condition_randomization":"NOT_DOCUMENTED_IN_THIS_PUBLIC_AGGREGATE",
      "contamination_negative_claim_valid":False,
      "microbiome_inference_from_10x_3prime":False,
      "min_dex_vehicle_barcode_ratio":min(x["dex_vehicle_barcode_ratio"] for x in data),
      "max_dex_vehicle_barcode_ratio":max(x["dex_vehicle_barcode_ratio"] for x in data),
      "strict_interpretation":"No negative microbial test can be inferred from absence in QC metadata. Technical library identity and treatment are aliased at library level; effects may combine DEX response and preparation effects."}
    save(out,"P6a_contamination_batch_gate.json",result)
    print("P6A",result["status"],"libraries",len(rows),flush=True)

def bh(pvals):
    order=sorted(range(len(pvals)),key=lambda i:pvals[i])
    q=[1.0]*len(pvals);running=1.
    for k in range(len(order)-1,-1,-1):
        i=order[k]
        running=min(running,pvals[i]*len(pvals)/(k+1))
        q[i]=min(1.,running)
    return q

def es_numpy(weights,hit_indices):
    import numpy as np
    n=len(weights); nh=len(hit_indices)
    if nh<1 or nh==n:raise ValueError("DEGENERATE_SET")
    hits=np.zeros(n,dtype=bool);hits[hit_indices]=True
    sw=float(weights[hits].sum())
    if sw<=0:raise ValueError("ZERO_HIT_WEIGHT")
    running=np.cumsum(np.where(hits,weights/sw,-1./(n-nh)))
    ix=int(np.argmax(np.abs(running)))
    return float(running[ix])

def pathway(inp,out,permutations):
    import numpy as np
    with gzip.open(inp,"rt",encoding="utf8",newline="") as f:
        primary=list(csv.DictReader(f,delimiter="\t"))
    if len(primary)!=16050:raise ValueError("P5B_EXPECTED_16050_FILTERED_GENES")
    bysymbol={}
    duplicates=Counter()
    for row in primary:
        symbol=row["gene_symbol"].strip()
        p=num(row["PValue"]);fc=num(row["logFC"]);fdr=num(row["FDR"])
        if not (0<=p<=1 and 0<=fdr<=1):raise ValueError("INVALID_DE_P_OR_FDR")
        stat=(1. if fc>0 else -1. if fc<0 else 0.)*min(50.,-math.log10(max(p,1e-300)))
        candidate={"gene_symbol":symbol,"gene_id":row["gene_id"],"score":stat,
           "primary_log2fc":fc,"primary_FDR":fdr}
        if symbol in bysymbol:
            duplicates[symbol]+=1
            if (abs(candidate["score"]),candidate["gene_id"]) > (abs(bysymbol[symbol]["score"]),bysymbol[symbol]["gene_id"]):
                bysymbol[symbol]=candidate
        else: bysymbol[symbol]=candidate
    ordered=sorted(bysymbol.values(),key=lambda x:(-x["score"],x["gene_id"]))
    names=[x["gene_symbol"] for x in ordered]
    symbols=set(names)
    assert len(ordered)==len(symbols)
    response=get(GMT_URL,limit=150_000)
    digest=hashlib.sha256(response).hexdigest()
    sets={}
    for line in response.decode("utf8").splitlines():
        fields=line.split("\t")
        if fields[0] in HALLMARK and len(fields)>=4:
            sets[fields[0]]=set(fields[2:])
    if set(sets)!=HALLMARK:raise ValueError("MSIGDB_NINE_SELECTED_HALLMARK_SETS_MISSING")
    scores=np.array([x["score"] for x in ordered],dtype="float64")
    weights=np.maximum(np.abs(scores),1e-12)
    rng=np.random.default_rng(20261031)
    result=[]
    for setname in sorted(HALLMARK):
        hit=[i for i,s in enumerate(names) if s in sets[setname]]
        nh=len(hit)
        if nh<10 or nh>len(names)*.25:
            raise ValueError("GENE_SET_OVERLAP_TOO_SMALL_OR_LARGE "+setname+" "+str(nh))
        observed=es_numpy(weights,hit)
        # Competitive gene-label permutation, preserving set membership count. NOT a donor/sample permutation.
        nul=np.empty(permutations,dtype="float64")
        for j in range(permutations):
            ix=rng.choice(len(names),size=nh,replace=False)
            nul[j]=es_numpy(weights,ix)
        p=(1+int(np.sum(np.abs(nul)>=abs(observed))))/(permutations+1)
        gene_members=sorted(TARGETS&sets[setname])
        measured=sorted(TARGETS&sets[setname]&symbols)
        result.append({"hallmark":setname,"msigdb_set_size":len(sets[setname]),
           "tested_universe_overlap":nh,"ES":round(observed,7),
           "gene_label_permutation_n":permutations,
           "two_sided_empirical_p":p,
           "target_12_overlap_symbols":",".join(measured),
           "target_12_overlap_count":len(measured),
           "test_is_exploratory_and_noncausal":True})
    adj=bh([x["two_sided_empirical_p"] for x in result])
    for row,q in zip(result,adj):row["BH_FDR_across_9_hallmarks"]=round(q,8)
    tsv(out,"P6b_9_hallmark_ranked_gene_set_enrichment.tsv",result)
    subset=[{"gene_symbol":x,"primary_log2fc":bysymbol[x]["primary_log2fc"],
         "primary_BH_FDR":bysymbol[x]["primary_FDR"],
         "hallmark_membership_count":sum(x in v for v in sets.values()),
         "hallmark_membership":",".join(sorted(k for k,v in sets.items() if x in v))}
        for x in sorted(TARGETS) if x in bysymbol]
    if len(subset)!=12:raise ValueError("12_TARGETS_NOT_ALL_IN_PRIMARY_UNIVERSE")
    tsv(out,"P6b_12_candidate_hallmark_membership.tsv",subset)
    summary={"status":"PASS_9_PUBLISHED_HALLMARKS_EXPLORATORY_RANKED_ENRICHMENT",
      "primary_source":"P5b existing 16050-gene edgeR aggregate, not fresh sequencing",
      "msigdb_release":"2025.1.Hs","msigdb_url":GMT_URL,
      "msigdb_hallmark_sha256":digest,"tested_genes_before_symbol_collapse":len(primary),
      "tested_unique_gene_symbols":len(symbols),
      "duplicate_gene_symbols_collapsed":sum(duplicates.values()),
      "gene_sets":9,"gene_label_permutations":permutations,
      "sets_BH_FDR_below_0_05":sum(x["BH_FDR_across_9_hallmarks"]<.05 for x in result),
      "not_sample_level_phenotype_permutation":True,
      "not_dose_specific":True,"not_microbiome_association":True,
      "not_independent_replication":True,
      "confounding_risk":"DEX samples are correlated with separate library preps at source-pool level",
      "interpretation":"Pre-ranked gene label permutation tests aggregate ranked association, not donor-wise randomized treatment effects or LPS/bacterial exposure."}
    save(out,"P6b_ranked_hallmark_summary.json",summary)
    print("P6B",summary["status"],"FDR_lt_.05",summary["sets_BH_FDR_below_0_05"],flush=True)

def ena_audit(out):
    results=[]
    for study in STUDIES:
        acc=study["accession"]
        row={"accession":acc,"scope":study["label"],"study_design":study["platform_claim"],
             "is_TM_DEX_study":False,"patient_matched_TM_Dex_phenotype":False,
             "paired_tear_proteome_raw_files_verified":False,
             "has_low_biomass_blank_controls_verified":False,
             "clinical_dex_exposure_verified":False}
        fields="run_accession,sample_accession,library_strategy,library_layout,fastq_ftp,fastq_md5,read_count,base_count,study_accession"
        q=urllib.parse.urlencode({"accession":acc,"result":"read_run","fields":fields,"format":"tsv"})
        url="https://www.ebi.ac.uk/ena/portal/api/filereport?"+q
        try:
            raw=get(url,limit=7_000_000,timeout=95)
            reader=csv.DictReader(io.StringIO(raw.decode("utf8","replace")),delimiter="\t")
            if not {"run_accession","sample_accession","fastq_ftp"}.issubset(reader.fieldnames or []):
                raise ValueError("EXPECTED_ENA_REPORT_FIELDS_NOT_FOUND")
            data=list(reader)
            if not data:raise ValueError("NO_RUNS_RETURNED")
            runs=[r["run_accession"] for r in data]
            samples=[r["sample_accession"] for r in data]
            if len(set(runs))!=len(runs) or len(set(samples))<1:
                raise ValueError("DUPLICATED_RUN_IDS_OR_NO_SAMPLES")
            strategies=Counter((r.get("library_strategy") or "UNKNOWN").upper() for r in data)
            downloads=[bool(r.get("fastq_ftp")) for r in data]
            pairedfastqs=sum(r["fastq_ftp"].count(";")>=1 for r in data)
            row.update({"metadata_status":"PASS_READ_RUN_METADATA_ONLY",
                "run_count":len(runs),"unique_biosample_count":len(set(samples)),
                "library_strategy_counts":dict(strategies),
                "fastq_link_run_count":sum(downloads),
                "fastq_link_coverage":round(sum(downloads)/len(data),4),
                "runs_with_multiple_fastq_links":pairedfastqs,
                "raw_fastq_downloaded":False,
                "metadata_download_sha256":hashlib.sha256(raw).hexdigest(),
                "metadata_only_dataset_go":sum(downloads)>0,
                "microbiome_TM_mechanistic_inference":"HOLD_NOT_MATCHED_TO_DEX_TM"})
        except Exception as e:
            row.update({"metadata_status":"HOLD_ENA_FETCH_OR_SCHEMA",
                "reason":type(e).__name__+": "+str(e)[:180],
                "metadata_only_dataset_go":False,"raw_fastq_downloaded":False,
                "microbiome_TM_mechanistic_inference":"HOLD"})
        results.append(row)
    # Some studies can legitimately remain HOLD, which must never be mislabeled as verified.
    save(out,"P6c_3_study_ena_metadata_verification.json",{
      "status":"PASS_AUDIT_COMPLETED_NOT_EQUIVALENT_TO_ALL_STUDIES_GO",
      "studies":results,"n_confirmed_metadata":sum(x["metadata_status"]=="PASS_READ_RUN_METADATA_ONLY" for x in results),
      "n_with_fastq_links":sum(x["metadata_only_dataset_go"] for x in results),
      "mechanistic_inference":"HOLD",
      "clinical_steroid_induced_IOP_prediction":"STOP",
      "no_patient_or_sample_level_metadata_uploaded":True,
      "low_biomass_blank_controls_require_paper_or_sample_level_verification":True,
      "cannot_identify_microbes_from_human_10x_3prime_pseudobulk":True})
    # Keep row-specific statuses in JSON because mixed nested numbers don't map cleanly to TSV.
    print("P6C study registry statuses",[(r["accession"],r["metadata_status"],r.get("run_count")) for r in results],flush=True)
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--mode",required=True,choices=["p6a","p6b","p6c"])
    a.add_argument("--input",type=Path)
    a.add_argument("--out",required=True,type=Path)
    a.add_argument("--permutations",type=int,default=999)
    opt=a.parse_args()
    opt.out.mkdir(parents=True,exist_ok=True)
    if opt.mode=="p6a":
        if not opt.input:a.error("--input needed for p6a")
        qc_source(opt.input,opt.out)
    elif opt.mode=="p6b":
        if not opt.input:a.error("--input needed for p6b")
        if not 99<=opt.permutations<=9999:a.error("--permutations 99..9999")
        pathway(opt.input,opt.out,opt.permutations)
    else:
        ena_audit(opt.out)
if __name__=="__main__":
    main()
