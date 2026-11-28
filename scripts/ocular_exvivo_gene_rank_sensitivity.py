#!/usr/bin/env python3
"""Post-readout sensitivity: collapse eligible Affymetrix probes per Entrez ID
before ranking across unique genes. Preserve primary result; no new independent
inference, score optimization, donor exclusions or cross-study pooling.
"""
import argparse,collections,csv,json,math,statistics
from pathlib import Path
from ocular_exvivo_source_gate import anno,require,sha
from ocular_exvivo_paired_score import pair_design,MATRIX,ANNOT,sign_p,holm
from ocular_frozen_241_external import read_panel,PIN
from publication_ocular_ckd_sources import read_table,percentile

def collapse_rank(by_gene):
    require(bool(by_gene),"Empty gene universe")
    keys=sorted(by_gene)
    arrays=[tuple(statistics.median(v) for v in zip(*by_gene[g])) for g in keys]
    return [dict(zip(keys,percentile(col))) for col in zip(*arrays)]

def run(root,primary):
    require(sha(root/"matrix.gz")==MATRIX and sha(root/"annot.gz")==ANNOT,"Input SHA changed")
    prior=json.loads(primary.read_text())
    require(prior["provenance"]["run_id"]==45535697612,"Primary run identity")
    p=prior["result"]
    require(p["status"]=="PAIRED_TISSUE_DESCRIPTIVE" and p["panel_sha256"]==PIN and p["matrix_sha256"]==MATRIX and p["annotation_sha256"]==ANNOT and p["pairs"]==5,"Primary evidence identity")
    panel=read_panel(root/"panel.zip");eligible,mapping,n=anno(root/"annot.gz",panel)
    meta,table=read_table(root/"matrix.gz","!series_matrix_table_begin","!series_matrix_table_end")
    h=next(table);ids=h[1:]
    require(h[0]=="ID_REF" and len(ids)==10 and set(ids)=={f"GSM{i}" for i in range(920734,920744)},"GEO sample set")
    pairs=pair_design(meta,ids)
    by_gene=collections.defaultdict(list);candidate_ids=collections.defaultdict(set);seen=set()
    for row in table:
        require(len(row)==len(h) and row[0] not in seen,"Matrix schema")
        probe=row[0];seen.add(probe)
        if probe not in eligible:continue
        vals=tuple(float(v) for v in row[1:]);require(all(math.isfinite(v) for v in vals),"Nonfinite probe")
        e=eligible[probe];by_gene[e].append(vals)
        if probe in mapping:candidate_ids[mapping[probe]].add(e)
    require(n==len(seen)==54675 and sum(len(v) for v in by_gene.values())==27920,"QC counts drift")
    require(all(len(v)==1 for v in candidate_ids.values()),"Candidate maps to multiple Entrez IDs")
    c={g:next(iter(v)) for g,v in candidate_ids.items()}
    require(len(set(c.values()))==len(c)==199,"Duplicate candidate Entrez or changed coverage")
    ranks=collapse_rank(by_gene);require(len(ranks)==10,"Rank columns")
    out={"status":"POST_READOUT_GENE_UNIVERSE_SENSITIVITY","study":"GSE37474","pairs":len(pairs),"panel_sha256":PIN,"matrix_sha256":MATRIX,"annotation_sha256":ANNOT,"primary_run":45535697612,"eligible_probes":27920,"unique_entrez_genes":len(by_gene),"probe_multiplicity":{"one_probe":sum(len(v)==1 for v in by_gene.values()),"multiple_probes":sum(len(v)>1 for v in by_gene.values()),"max_probes":max(map(len,by_gene.values()))},"panels":{},"limits":["Same 5 donor pairs and normalized array signals; not independent replication or an alternative confirmatory test.","Median within gene across strict eligible probe sets then average-rank across unique Entrez genes. Platform-specific probe affinity remains; no expression fold-change interpretation.","Preserve primary probe-rank score; this sensitivity was designed after primary result readout. Do not select method by effect magnitude or P. No additional inferential claim."]}
    for name,genes,expected in [("all241",list(panel),199),("ECM16",[g for g,v in panel.items() if v["ecm_member"]],14)]:
        mapped=[g for g in genes if g in c];require(len(mapped)==expected==p["panels"][name]["mapped"],"Fixed mapping drift")
        for s in (-1,1):
            require(sum((panel[g]["acute_log2fc"]>0)==(s==1) for g in mapped)>=math.ceil(.7*sum((panel[g]["acute_log2fc"]>0)==(s==1) for g in genes)),"Sign coverage")
        vec=[[(ranks[v["dex"]][c[g]]-ranks[v["control"]][c[g]])*(1 if panel[g]["acute_log2fc"]>0 else -1) for g in mapped] for v in pairs]
        eff=[statistics.mean(v) for v in vec];loo=[statistics.mean(eff[:i]+eff[i+1:]) for i in range(5)]
        out["panels"][name]={"mapped":len(mapped),"positive_pairs":sum(v>0 for v in eff),"negative_pairs":sum(v<0 for v in eff),"ties":sum(v==0 for v in eff),"median_signed_percentile_shift":statistics.median(eff),"mean_signed_percentile_shift":statistics.mean(eff),"pair_scores_sorted":sorted(eff),"leave_one_pair_out_mean_range":[min(loo),max(loo)],"positive_equal_pair_mean_genes":sum(statistics.mean(v)>0 for v in zip(*vec)),"all_five_positive_genes":sum(all(x>0 for x in v) for v in zip(*vec)),"primary_positive_pairs":p["panels"][name]["positive_pairs"],"primary_median_shift":p["panels"][name]["median_signed_percentile_shift"]}
    return out

def main():
    a=argparse.ArgumentParser();a.add_argument("--root",type=Path);a.add_argument("--primary",type=Path);a.add_argument("--out",type=Path);a.add_argument("--selftest",action="store_true");p=a.parse_args()
    if p.selftest:
        x=collapse_rank({"1":[(1.,3.),(5.,7.)],"2":[(7.,1.)]})
        assert x==[{"1":.25,"2":.75},{"1":.75,"2":.25}]
        assert collapse_rank({"1":[(1.,)],"2":[(1.,)]})==[{"1":.5,"2":.5}]
        print("Focused unique-gene rank tests passed; no biological result.");return
    try:x=run(p.root,p.primary)
    except (OSError,ValueError,KeyError,TypeError,StopIteration,EOFError) as e:x={"status":"HOLD_SENSITIVITY_SOURCE_OR_SCHEMA","reason":str(e)[:180]}
    p.out.parent.mkdir(parents=True,exist_ok=True);p.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
