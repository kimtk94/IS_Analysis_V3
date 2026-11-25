#!/usr/bin/env python3
"""Local QC of authorized aggregate external responses. No inference or network."""
import argparse,collections,csv,json,math,re
from pathlib import Path

MARKER_FIELDS="marker_id rsid requested_build chromosome position requested_effect_allele requested_other_allele historical_GRCh37_key".split()
FIELDS="marker_id rsid requested_build requested_chr requested_pos requested_effect_allele requested_other_allele phenotype returned_build returned_chr returned_pos returned_effect_allele returned_other_allele beta se p_two_sided effect_allele_frequency frequency_population imputation_metric imputation_value n_cases n_controls per_marker_n model_id status reason_if_unavailable".split()
PHENOTYPES=("AIS","LAS","CES","SVS")
STAT_FIELDS=("beta","se","p_two_sided","effect_allele_frequency")
UNAVAILABLE={"NOT_RETURNED","NOT_AVAILABLE","QC_FAIL","SUPPRESSED"}
def require(ok,msg):
    if not ok:raise ValueError(msg)
def read_tsv(path,fields,cap):
    require(path.stat().st_size<=2000000,"FILE_SIZE")
    with path.open(encoding="utf-8-sig",newline="") as f:
        t=csv.DictReader(f,delimiter="\t")
        require(t.fieldnames==fields,"HEADER_MISMATCH")
        rows=[]
        for d in t:
            require(len(rows)<cap and None not in d and all(isinstance(v,str) and len(v)<=1000 for v in d.values()),"ROW_WIDTH_OR_CAP")
            rows.append(d)
    return rows
def number(v,code):
    try:x=float(v)
    except (ValueError,TypeError):raise ValueError(code) from None
    require(math.isfinite(x),code);return x
def count(v,code):
    require(re.fullmatch(r"[1-9][0-9]{0,9}",v) is not None,code);return int(v)
def marker_contract(packet):
    rows=read_tsv(packet/"marker_request.tsv",MARKER_FIELDS,24)
    a=json.loads((packet/"anchor_readback.json").read_text())
    require(a["provenance"]["run_id"]==45534891070 and a["provenance"]["artifact_id"]==14986197928 and a["result"]["verified_count"]==24,"ANCHOR_PROVENANCE")
    raw=a["result"]["rows"]
    require(len(rows)==len(raw)==24,"MARKER_COUNT")
    out={};rsids=set();keys=set();chrom=set()
    for row in rows:
        k=row["marker_id"];rs=row["rsid"];c=row["chromosome"];p=row["position"];ea=row["requested_effect_allele"];oa=row["requested_other_allele"]
        require(k in raw and k not in out and rs not in rsids and re.fullmatch(r"rs[0-9]+",rs) is not None,"MARKER_IDENTITY")
        d=raw[k]
        require(row["requested_build"]=="GRCh38" and c==d["GRCh38_chr"] and p==str(d["GRCh38_pos"]) and rs==d["rsid"] and ea==d["historical_ALT"] and oa==d["historical_REF"] and row["historical_GRCh37_key"]==d["historical_parity"]["GRCh37_variant_key"] and d["historical_parity"]["joint_effect_p_eaf_reproduced"] is True,"ANCHOR_IDENTITY_DRIFT")
        require("".join(sorted(ea+oa)) in {"AC","AG","CT","GT"} and (c,p,ea,oa) not in keys,"ALLELE_OR_DUPLICATE")
        out[k]=[k,rs,"GRCh38",c,p,ea,oa];rsids.add(rs);keys.add((c,p,ea,oa));chrom.add(c)
    require(len(chrom)==16,"CHROMOSOME_COUNT")
    return out
def model_policy(path):
    if path is None:return {}
    require(path.stat().st_size<=100000,"POLICY_SIZE")
    d=json.loads(path.read_text(encoding="utf8"))
    require(isinstance(d,dict) and set(d)=={"models"} and isinstance(d["models"],dict) and len(d["models"])<=20,"POLICY_SCHEMA")
    for k,v in d["models"].items():
        require(isinstance(k,str) and re.fullmatch(r"[A-Za-z0-9_.-]{1,100}",k) and isinstance(v,dict) and set(v)=={"phenotype","imputation_metric","min_quality","min_cases","min_controls"},"MODEL_POLICY_SCHEMA")
        require(v["phenotype"] in PHENOTYPES and isinstance(v["imputation_metric"],str) and 0<len(v["imputation_metric"])<=100,"MODEL_POLICY_ID")
        q=v["min_quality"];require(type(q) in (int,float) and math.isfinite(q) and 0<=q<=1,"MODEL_POLICY_QUALITY")
        require(all(type(v[x]) is int and v[x]>0 for x in ("min_cases","min_controls")),"MODEL_POLICY_N")
    return d["models"]
def qc_row(d,models):
    status=d["status"]
    if status in UNAVAILABLE:
        require(all(not d[k] for k in STAT_FIELDS),"UNAVAILABLE_WITH_EFFECTS")
        if status=="NOT_RETURNED":
            require(all(not d[k] for k in FIELDS[8:24]),"UNRETURNED_WITH_DATA")
        else:require(bool(d["reason_if_unavailable"].strip()),"UNAVAILABLE_REASON_MISSING")
        return status,None
    require(status=="OK" and not d["reason_if_unavailable"].strip(),"STATUS_INVALID")
    require(d["returned_build"]==d["requested_build"] and d["returned_chr"]==d["requested_chr"] and d["returned_pos"]==d["requested_pos"],"RETURNED_COORDINATE_MISMATCH")
    pair=(d["returned_effect_allele"],d["returned_other_allele"]);target=(d["requested_effect_allele"],d["requested_other_allele"])
    require(pair in (target,target[::-1]),"RETURNED_ALLELE_MISMATCH")
    sign=1 if pair==target else -1
    b=number(d["beta"],"BETA_INVALID");se=number(d["se"],"SE_INVALID");p=number(d["p_two_sided"],"P_INVALID");eaf=number(d["effect_allele_frequency"],"EAF_INVALID")
    require(se>0,"SE_INVALID");require(0<p<=1,"P_INVALID_OR_UNDERFLOW");require(0<eaf<1,"EAF_INVALID")
    require(bool(d["frequency_population"].strip()),"FREQUENCY_POPULATION_MISSING")
    q=number(d["imputation_value"],"QUALITY_INVALID")
    require(bool(d["imputation_metric"].strip()) and 0<=q<=1,"QUALITY_INVALID")
    ncase=count(d["n_cases"],"CASE_COUNT_INVALID");nctrl=count(d["n_controls"],"CONTROL_COUNT_INVALID");n=count(d["per_marker_n"],"MARKER_N_INVALID")
    require(2<=n<=ncase+nctrl,"MARKER_N_EXCEEDS_TOTAL_OR_TOO_SMALL")
    require(re.fullmatch(r"[A-Za-z0-9_.-]{1,100}",d["model_id"]) is not None,"MODEL_ID_INVALID")
    policy=models.get(d["model_id"])
    label="STRUCTURAL_OK_POLICY_PENDING"
    if policy is not None:
        require(policy["phenotype"]==d["phenotype"] and policy["imputation_metric"]==d["imputation_metric"],"MODEL_POLICY_MISMATCH")
        require(q>=policy["min_quality"],"QUALITY_BELOW_PRELOCKED_THRESHOLD")
        require(ncase>=policy["min_cases"] and nctrl>=policy["min_controls"],"N_BELOW_PRELOCKED_THRESHOLD")
        label="STRUCTURAL_AND_CONFIGURED_QC_OK"
    normalized={"marker_id":d["marker_id"],"rsid":d["rsid"],"phenotype":d["phenotype"],"beta_requested_allele":sign*b,"se":se,"p_two_sided":p,"eaf_requested_allele":eaf if sign==1 else 1-eaf,"orientation":"DIRECT" if sign==1 else "REVERSED","frequency_population":d["frequency_population"],"model_id":d["model_id"],"per_marker_n":n}
    return label,normalized
def run(packet,response,policy_path=None):
    contract=marker_contract(packet);models=model_policy(policy_path)
    rows=read_tsv(response,FIELDS,96)
    expected={(k,p) for k in contract for p in PHENOTYPES};seen=set();tallies=collections.Counter();byphen={p:collections.Counter() for p in PHENOTYPES};issues=[];norm=[]
    for d in rows:
        key=(d["marker_id"],d["phenotype"])
        require(key in expected and key not in seen,"UNEXPECTED_OR_DUPLICATE_REQUEST_ROW")
        seen.add(key)
        require([d[k] for k in FIELDS[:7]]==contract[key[0]],"REQUEST_IDENTITY_DRIFT")
        try:
            state,n=qc_row(d,models)
            if n is not None:norm.append(n)
        except ValueError as e:
            state="ROW_QC_HOLD";issues.append({"marker_id":key[0],"phenotype":key[1],"reason_code":str(e)})
        tallies[state]+=1;byphen[key[1]][state]+=1
    require(seen==expected,"INCOMPLETE_FIXED_RESPONSE_GRID")
    report={"status":"AGGREGATE_STRUCTURE_QC_ONLY","expected_rows":96,"received_rows":len(rows),"counts":dict(tallies),"by_phenotype":{k:dict(v) for k,v in byphen.items()},"issues":issues,"configured_model_count":len(models),"reversed_alleles":sum(n["orientation"]=="REVERSED" for n in norm),"independent_validation_status":"HOLD","limits":["No phenotype, model, recruitment, overlap, authorization, effective power or full reference sequence certification by this parser.","No threshold inferred from external outcomes. Configured QC is conditional on an independently agreed model policy; this tool does not certify that policy.","No directional or significance test, no OR conversion, no meta-analysis. All 96 requested rows remain accounted for. The external response is not used by GitHub Actions."]}
    return report,norm
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--packet",type=Path,required=True);p.add_argument("--response",type=Path,required=True);p.add_argument("--policy",type=Path);p.add_argument("--out",type=Path,required=True)
    p.add_argument("--local-normalized-out",type=Path,help="Optional sensitive aggregate values. Only use approved local storage, never GitHub/Actions.")
    a=p.parse_args()
    try:r,n=run(a.packet,a.response,a.policy)
    except (OSError,ValueError,KeyError,TypeError,csv.Error) as e:
        # Do not echo untrusted record text or paths on failure.
        r={"status":"FILE_QC_HOLD","reason_code":str(e) if type(e) is ValueError and re.fullmatch(r"[A-Z0-9_]+",str(e)) else type(e).__name__,"independent_validation_status":"HOLD"};n=[]
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
    if a.local_normalized_out and r["status"]=="AGGREGATE_STRUCTURE_QC_ONLY":
        a.local_normalized_out.write_text(json.dumps(n,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":r["status"],"expected_rows":r.get("expected_rows"),"received_rows":r.get("received_rows"),"counts":r.get("counts"),"independent_validation_status":"HOLD"},indent=2))
    return int(r["status"]=="FILE_QC_HOLD")
if __name__=="__main__":raise SystemExit(main())
