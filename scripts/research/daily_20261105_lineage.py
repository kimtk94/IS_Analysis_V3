#!/usr/bin/env python3
"""Cross-day evidence lineage and scientific gate audit. GitHub Actions only."""
import argparse,datetime,hashlib,json,os,pathlib,subprocess
TRACKS=("is","glaucoma","hcc","radiomics","ocular","muscle","oral","skin")
NEXT={
"is":"Full GWAS checksum, allele harmonization, EUR/EAS LD and MR replication",
"glaucoma":"Donor-level raw counts, 9x leave-pool-out QL and independent strain",
"hcc":"Official 330 quant-column to aliquot mapping and RNA patient join",
"radiomics":"Original NIfTI spacing, clinical eligibility and MRI-RNA patient join",
"ocular":"Full FASTQ checksums and unique timepoint/strain mapping",
"muscle":"Resolve 2 missing loci and validate EAS LD and muscle eQTL",
"oral":"Actual expression matrix and 18 donor clinical identifiers",
"skin":"Complete GWAS SNP QC, genome build and MR instrument validation"}
def download(run,name,target):
    target.mkdir(parents=True,exist_ok=True)
    cmd=["gh","run","download",str(run),"--repo","kimtk94/IS_Analysis_V3","--name",name,"--dir",str(target)]
    p=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
    if p.returncode: return {"state":"HOLD_DOWNLOAD","error":p.stderr[-250:]}
    files=sorted(target.rglob("*.json"))
    if len(files)!=1:return {"state":"HOLD_JSON_CARDINALITY","count":len(files)}
    raw=files[0].read_bytes()
    try: obj=json.loads(raw)
    except ValueError:return {"state":"FAIL_MALFORMED_JSON"}
    if not isinstance(obj,dict):return {"state":"FAIL_NOT_OBJECT"}
    return {"state":"PASS_ARTIFACT_INTEGRITY","sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw),"content":obj}
def main():
    a=argparse.ArgumentParser();a.add_argument("--track",choices=TRACKS,required=True);v=a.parse_args().track
    if not os.environ.get("GITHUB_ACTIONS"):raise SystemExit("GitHub Actions only")
    out=pathlib.Path("results/research_operations/20261105");out.mkdir(parents=True,exist_ok=True)
    base=pathlib.Path("/tmp/evidence-lineage")/v
    old=download(41544001319,"research-20261103-"+v,base/"nov03")
    new=download(41742724187,"research-20261104-"+v,base/"nov04")
    issues=[];changes={}
    for day,d in (("20261103",old),("20261104",new)):
        if d["state"].startswith("FAIL"):issues.append(day+":"+d["state"])
        if d["state"].startswith("HOLD"):issues.append(day+":"+d["state"])
        if d["state"]=="PASS_ARTIFACT_INTEGRITY":
            obj=d.pop("content")
            # Preserve only non-identifying keys and declared statuses.
            changes[day]={"fields":sorted(obj.keys()),"scientific_gate":obj.get("scientific_gate"),
              "probe_status":obj.get("probe_status"),"status":obj.get("status"),
              "track":obj.get("track")}
            if obj.get("track")!=v:issues.append(day+":TRACK_ID_MISMATCH")
            if obj.get("scientific_gate")=="PASS":issues.append(day+":UNSUPPORTED_SCIENTIFIC_PASS")
    result={"date_kst":"2026-11-05","track":v,"task":"Independent cross-day artifact lineage and gate-drift check",
      "source_runs":[41544001319,41742724187],"artifact_20261103":old,
      "artifact_20261104":new,"comparison":changes,
      "evidence_status":"PASS" if not issues else "HOLD_REVIEW",
      "scientific_gate":"HOLD","issues":issues,"next_gate":NEXT[v],
      "note":"Artifact consistency is not independent biological validation."}
    path=out/(v+".json");path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"track":v,"evidence_status":result["evidence_status"],"issues":issues}))
    if any("FAIL" in i or "MISMATCH" in i for i in issues):raise SystemExit(1)
if __name__=="__main__":main()
