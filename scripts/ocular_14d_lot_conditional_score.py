#!/usr/bin/env python3
"""GSE65240 conditional descriptive paired ranks under explicit depositor labels.
Primary source conflict remains HOLD; no donor independence or inferential P.
"""
import argparse,collections,json,math,statistics,re
from pathlib import Path
from ocular_14d_lot_source_gate import annotation,design,PLATFORMS,EXPECTED
from ocular_frozen_241_external import read_panel,PIN
from publication_ocular_ckd_sources import read_table,percentile
from ocular_exvivo_source_gate import metadata,require,sha
def score(root,evidence):
    old=json.loads((evidence/"ocular_14d_lot_source_gate_20261125.json").read_text())
    label=json.loads((evidence/"ocular_14d_lot_label_crosscheck_20261125.json").read_text())
    require(old["provenance"]["run_id"]==45538483133 and label["provenance"]["run_id"]==45561823908,"Pinned provenance")
    gate=old["result"];rec=label["filename_reconciliation"]
    require(rec["filename_matches"]==6 and rec["source_name_conflicts"]==4 and rec["primary_scoring"]=="HOLD","Conditional label contract")
    require(gate["panel_sha256"]==PIN and gate["source_name_conflict_count"]==4 and gate["same_platform_lot_pairs"]==3,"Prior blinded gate")
    samples={d["gsm"]:d for d in rec["samples"]}
    require(set(samples)==set(EXPECTED),"Exact sample set")
    panel=read_panel(root/"panel.zip")
    platforms={};mappings={};ranks={};pairs={}
    for plat in PLATFORMS:
        p=gate["platforms"][plat]
        mat=root/(plat+"-matrix.gz");ann=root/(plat+"-annot.gz")
        require(sha(mat)==p["matrix_sha256"] and sha(ann)==p["annotation_sha256"],"Original source changed "+plat)
        eligible,mp,n=annotation(ann,panel,plat)
        meta,table=read_table(mat,"!series_matrix_table_begin","!series_matrix_table_end")
        h=next(table);ids=h[1:]
        require(h[0]=="ID_REF" and len(ids)==len(set(ids)),"Matrix header")
        verified=design(meta,ids,plat);require(verified==p["samples"],"Blinded design drift")
        sup=metadata(meta,"!Sample_supplementary_file")
        require(len(sup)==1 and len(sup[0])==len(ids),"Supplementary filename width")
        for j,d in enumerate(verified):
            gsm=d["gsm"];z=samples[gsm];lot,arm=EXPECTED[gsm]
            name=sup[0][j].rsplit("/",1)[-1]
            m=re.fullmatch(re.escape(gsm)+r"_(2584|3423|4973)_?(control|DEX)\.txt\.gz",name,re.I)
            require(m and m[1]==lot and m[2].lower()==arm and z["filename"]==name and z["lot"]==lot and z["arm"]==arm and z["platform"]==plat,"Source label mismatch")
            require(z["source_conflict"]== (not d["source_name_lot_matches"]),"Source name conflict status drift")
            pair=pairs.setdefault(lot,{"platform":plat})
            require(pair["platform"]==plat and arm not in pair,"Cross-platform or duplicate arm")
            pair[arm]=(plat,j)
        values=[];seen=set();gp=collections.defaultdict(list);mapped=collections.defaultdict(set)
        for row in table:
            require(len(row)==len(h) and row[0] not in seen,"Matrix row schema")
            probe=row[0];seen.add(probe)
            if probe not in eligible:continue
            vals=[float(v) for v in row[1:]]
            require(all(math.isfinite(v) for v in vals),"Nonfinite source")
            if probe in mp:
                g=mp[probe];gp[g].append(len(values));mapped[g].add(eligible[probe])
            values.append(vals)
        require(len(seen)==p["matrix_rows"] and n==p["annotation_rows"] and len(values)==len(eligible)==p["complete_rank_universe"]==29831 and len(mapped)==p["mapped_candidates"]==220,"Blinded counts changed")
        require(all(len(v)==1 for v in mapped.values()),"Ambiguous candidate Entrez")
        mappings[plat]={g:next(iter(v)) for g,v in mapped.items()}
        require(len(set(mappings[plat].values()))==len(mapped),"Duplicate Entrez")
        ranks[plat]=[]
        for col in zip(*values):
            pc=percentile(col)
            ranks[plat].append({g:statistics.median(pc[j] for j in ix) for g,ix in gp.items()})
        platforms[plat]={"samples":ids,"rank_universe":len(values),"mapped_candidates":len(mapped),"matrix_sha256":sha(mat),"annotation_sha256":sha(ann)}
    require(set(pairs)=={"2584","3423","4973"} and all(set(v)=={"platform","control","dex"} for v in pairs.values()),"All fixed lot pairs required")
    a,b=[mappings[v] for v in PLATFORMS]
    common={g for g,v in a.items() if b.get(g)==v}
    out={"status":"CONDITIONAL_DEPOSITOR_LABEL_DESCRIPTIVE","study":"GSE65240","primary_source_gate":"HOLD_SOURCE_NAME_CONFLICT","source_conflicts":4,"commercial_lots":3,"donor_independence_certified":False,"panel_sha256":PIN,"platforms":platforms,"panels":{},"limits":["Title, explicit characteristics and original raw filenames corroborate all six intended sample labels; source_name conflicts in four. Conditional interpretation, not adjudicated chain of custody. No outcome-based label changes.","Use three fixed commercial lot pairs, each within a single platform, and conflict-free lot4973 separately. No cross-platform array pooling. Probe percentile ranks within each platform/sample; gene score median across mapped probes; common exact Entrez panel frozen before scoring.","No formal P values because independent biological units are not certified. Three commercial lots are not three independently verified donors. Conflict-free pair is one descriptive observation, not independent validation.","Same public study examined previously; no claim of blind validation, clinical prediction, drug equivalence, or pooled significance."]}
    for name,genes in [("all241",list(panel)),("ECM16",[g for g,v in panel.items() if v["ecm_member"]])]:
        gg=[g for g in genes if g in common]
        expected=gate["common_panels"][name]
        require(len(gg)==expected["mapped_common_exact_entrez"] and expected["coverage_gate"],"Frozen common mapping")
        counts={str(d):{"frozen":sum((panel[g]["acute_log2fc"]>0)==(d==1) for g in genes),"mapped":sum((panel[g]["acute_log2fc"]>0)==(d==1) for g in gg)} for d in (-1,1)}
        require(counts==expected["direction_coverage"],"Direction coverage drift")
        lot_scores={}
        for lot,pair in sorted(pairs.items()):
            pl=pair["platform"];c=ranks[pl][pair["control"][1]];t=ranks[pl][pair["dex"][1]]
            v=[(t[g]-c[g])*(1 if panel[g]["acute_log2fc"]>0 else -1) for g in gg]
            lot_scores[lot]={"platform":pl,"signed_percentile_shift":statistics.mean(v),"positive_direction_genes":sum(x>0 for x in v),"source_name_conflict":lot!="4973"}
        vals=[v["signed_percentile_shift"] for v in lot_scores.values()]
        out["panels"][name]={"frozen":len(genes),"mapped":len(gg),"direction_coverage":counts,"all_three_lots":lot_scores,"positive_lots":sum(v>0 for v in vals),"median_conditional_shift":statistics.median(vals),"mean_conditional_shift":statistics.mean(vals),"conflict_free_single_pair_shift":lot_scores["4973"]["signed_percentile_shift"],"inferential_p":None}
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--evidence",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:x=score(a.root,a.evidence)
    except (OSError,ValueError,KeyError,TypeError,StopIteration,EOFError,IndexError) as e:x={"status":"HOLD_SOURCE_OR_SCHEMA","reason":str(e)[:250]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n");print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
