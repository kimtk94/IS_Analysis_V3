#!/usr/bin/env python3
"""GSE37474 fixed-panel paired rank score. No outcome-based exclusions."""
import argparse, collections, csv, json, math, re, statistics
from pathlib import Path
from ocular_exvivo_source_gate import anno, metadata, require, sha
from ocular_frozen_241_external import read_panel,PIN
from publication_ocular_ckd_sources import read_table,percentile
MATRIX="d348cb1add4a01ec1b17d82da75f4b2081de4d7c3ad2433a42261fdd9543c849"
ANNOT="d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394"

def title_pair(title):
    # Exact original title grammar independently read from the pinned GEO matrix
    # in Actions run 45535645739; no outcome data accessed in that readback.
    m=re.fullmatch(r"(Control|Dex treatment), donor ([1-5])",title.strip())
    require(m is not None,"Unrecognized donor/arm title")
    arm,n=m.groups()
    return int(n),"control" if arm=="Control" else "dex"

def pair_design(meta,ids):
    require(metadata(meta,"!Series_geo_accession")==[["GSE37474"]],"Wrong series")
    require(metadata(meta,"!Sample_geo_accession")==[ids] and metadata(meta,"!Sample_platform_id")==[["GPL570"]*10],"Wrong sample/platform order")
    titles=metadata(meta,"!Sample_title")
    require(len(titles)==1 and len(titles[0])==10,"Missing titles")
    chars={}
    for row in metadata(meta,"!Sample_characteristics_ch1"):
        require(len(row)==10,"Characteristics width")
        parts=[v.split(":",1) for v in row]
        require(all(len(x)==2 for x in parts),"Missing key/value characteristics")
        keys={x[0].strip().lower() for x in parts}
        require(len(keys)==1,"Mixed characteristic fields")
        k=next(iter(keys));require(k not in chars,"Duplicate characteristic")
        chars[k]=[x[1].strip().lower() for x in parts]
    require({"treatment","age","sex"}<=set(chars),"Missing characteristics")
    pairs=collections.defaultdict(dict)
    for i,title in enumerate(titles[0]):
        d,a=title_pair(title)
        require(a not in pairs[d],"Repeated donor arm")
        tr=chars["treatment"][i]
        expected="none" if a=="control" else "100 nm dexamethasone"
        require(tr==expected,"Treatment field does not corroborate title")
        pairs[d][a]=i
    require(set(pairs)==set(range(1,6)) and all(set(v)=={"control","dex"} for v in pairs.values()),"Expected five exact title-defined donor pairs")
    require(all(ids[v["control"]]==f"GSM{920734+2*(k-1)}" and ids[v["dex"]]==f"GSM{920735+2*(k-1)}" for k,v in pairs.items()),"Source GSM and donor/arm mapping drift")
    for v in pairs.values():
        for k in ("age","sex"):
            require(chars[k][v["control"]]==chars[k][v["dex"]] and chars[k][v["control"]] not in {"","na","n/a","unknown"},"Donor demographics inconsistent or missing")
    return [pairs[k] for k in sorted(pairs)]

def sign_p(values):
    pos=sum(v>0 for v in values);neg=sum(v<0 for v in values);n=pos+neg
    return min(1.,2*sum(math.comb(n,k) for k in range(min(pos,neg)+1))/2**n) if n else 1.
def holm(p):
    out={};running=0.
    for i,(k,v) in enumerate(sorted(p.items(),key=lambda x:x[1])):
        running=max(running,min(1.,(len(p)-i)*v));out[k]=running
    return out
def score(root):
    require(sha(root/"matrix.gz")==MATRIX and sha(root/"annot.gz")==ANNOT,"Previously verified inputs changed")
    panel=read_panel(root/"panel.zip")
    eligible,mapping,n=anno(root/"annot.gz",panel)
    meta,table=read_table(root/"matrix.gz","!series_matrix_table_begin","!series_matrix_table_end")
    h=next(table);ids=h[1:]
    require(h[0]=="ID_REF" and len(ids)==10 and set(ids)=={f"GSM{i}" for i in range(920734,920744)},"Exact GEO sample set")
    pairs=pair_design(meta,ids)
    values=[];seen=set();gene_probes=collections.defaultdict(list)
    for row in table:
        require(len(row)==len(h) and row[0] not in seen,"Duplicate/malformed row")
        probe=row[0];seen.add(probe)
        if probe not in eligible:continue
        vals=[float(v) for v in row[1:]]
        require(all(math.isfinite(v) for v in vals),"Unexpected nonfinite source value")
        if probe in mapping:gene_probes[mapping[probe]].append(len(values))
        values.append(vals)
    require(n==len(seen)==54675 and len(values)==len(eligible)==27920,"Blinded QC counts changed")
    ranks=[]
    for col in zip(*values):
        v=percentile(col)
        ranks.append({g:statistics.median(v[j] for j in ix) for g,ix in gene_probes.items()})
    out={"status":"PAIRED_TISSUE_DESCRIPTIVE","study":"GSE37474","pairing_verified_by_public_titles_treatment_age_sex":True,"pairs":5,"panel_sha256":PIN,"matrix_sha256":MATRIX,"annotation_sha256":ANNOT,"rank_universe":len(values),"panels":{},"limitations":["Mixed TM/corneoscleral tissue, not purified TM or IOP prediction. Prior literature has reanalyzed this dataset; independence from original panel selection and donor cohorts not fully certified.","Five donor pairs are the analysis units. Rank shift estimates are not log2FC; gene counts are descriptive, not independent tests. No pooled donor or gene inference.","Fixed 75% overall, 70% directional coverage gates; no retuning. Two-sided exact sign tests with Holm over two panels; no pooling with other studies to force significance."]}
    for name,genes,expected in [("all241",list(panel),199),("ECM16",[g for g,v in panel.items() if v["ecm_member"]],14)]:
        mapped=[g for g in genes if g in gene_probes]
        require(len(mapped)==expected and len(mapped)>=math.ceil(.75*len(genes)),"Coverage changed")
        for s in (-1,1):
            require(sum((panel[g]["acute_log2fc"]>0)==(s==1) for g in mapped)>=math.ceil(.7*sum((panel[g]["acute_log2fc"]>0)==(s==1) for g in genes)),"Sign coverage")
        effects=[];vectors=[]
        for v in pairs:
            d=[(ranks[v["dex"]][g]-ranks[v["control"]][g])*(1 if panel[g]["acute_log2fc"]>0 else -1) for g in mapped]
            effects.append(statistics.mean(d));vectors.append(d)
        loo=[statistics.mean(effects[:i]+effects[i+1:]) for i in range(5)]
        out["panels"][name]={"frozen":len(genes),"mapped":len(mapped),"positive_pairs":sum(v>0 for v in effects),"negative_pairs":sum(v<0 for v in effects),"ties":sum(v==0 for v in effects),"median_signed_percentile_shift":statistics.median(effects),"mean_signed_percentile_shift":statistics.mean(effects),"pair_score_distribution_sorted":sorted(effects),"leave_one_pair_out_mean_range":[min(loo),max(loo)],"positive_equal_pair_mean_genes":sum(statistics.mean(v)>0 for v in zip(*vectors)),"all_five_positive_genes":sum(all(x>0 for x in v) for v in zip(*vectors)),"two_sided_exact_sign_p":sign_p(effects)}
    adj=holm({k:v["two_sided_exact_sign_p"] for k,v in out["panels"].items()})
    for k,v in out["panels"].items():v["holm_p"]=adj[k]
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path);p.add_argument("--out",type=Path);p.add_argument("--selftest",action="store_true");a=p.parse_args()
    if a.selftest:
        assert title_pair("Control, donor 1")== (1,"control")
        assert title_pair("Dex treatment, donor 2")== (2,"dex")
        assert sign_p([1]*5)==.0625
        assert holm({"a":.0625,"b":.0625})=={"a":.125,"b":.125}
        assert percentile([4,4,1,9])==[.5,.5,.125,.875]
        print("Focused pairing/score/statistic checks passed; no biological result.");return
    try:x=score(a.root)
    except (OSError,ValueError,KeyError,TypeError,StopIteration,EOFError) as e:x={"status":"HOLD_SOURCE_OR_PAIRING","reason":str(e)[:160]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
