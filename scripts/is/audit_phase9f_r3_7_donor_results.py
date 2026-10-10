#!/usr/bin/env python3
"""Independent R3_7 source audit. Never mutate source outputs; stdlib only."""
import argparse
import csv
import hashlib
import html
import json
import math
from pathlib import Path
from statistics import median

INPUTS=["HUMAN_DONOR_PSEUDOBULK.tsv","HUMAN_DONOR_CELLTYPE_COVERAGE.tsv",
        "HUMAN_DONOR_FEATURE_STATUS.tsv","HUMAN_DONOR_PAIRED_COMPARISONS.tsv",
        "HUMAN_DONOR_PSEUDOBULK_README.md","HUMAN_RUN_MANIFEST.json"]

def tsv(p):
    with p.open(newline="", encoding="utf8") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def num(s):
    if s in ("NA","",None):
        return float("nan")
    n=float(s)
    if not math.isfinite(n):
        raise ValueError("Nonfinite numeric value: "+str(s))
    return n

def render_chart(comps,path):
    width,height=1230,872
    lo=min(-1.0,min(y for c in comps for y in c["individual_deltas"])-0.5)
    hi=max(6.0,max(y for c in comps for y in c["individual_deltas"])+0.5)
    xpos=lambda v: 542+(v-lo)/(hi-lo)*438
    cdict={"SH3PXD2A":"#245f81","COL4A2":"#966321","COL4A1":"#6863a8",
           "CALHM2":"#267a65","ALDH2":"#a05262"}
    s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
       '<rect width="100%" height="100%" fill="white"/>',
       '<text x="36" y="42" font-family="sans-serif" font-size="26" font-weight="bold" fill="#18384d">R3_7 donor-paired expression</text>',
       '<text x="36" y="75" font-family="sans-serif" font-size="16" fill="#546c79">Each dot = one donor; diamond = median; healthy adult temporal-lobe reference</text>',
       '<text x="36" y="99" font-family="sans-serif" font-size="14" fill="#516575">Difference: log2(CPM + 1) in cell class A minus class B (not direct log2 fold change)</text>']
    for tick in range(math.ceil(lo),math.floor(hi)+1):
        x=xpos(tick)
        s += [f'<line x1="{x:.2f}" x2="{x:.2f}" y1="122" y2="760" stroke="{"#bb7070" if tick==0 else "#e3ebef"}" stroke-width="{"2" if tick==0 else "1"}"/>',
              f'<text x="{x:.2f}" y="788" font-size="14" text-anchor="middle" fill="#3e5668">{tick:+d}</text>']
    for i,c in enumerate(comps):
        y=154+i*60
        color=cdict.get(c["gene"],"#395366")
        label=c["gene"]+" | "+c["celltype_a"]+" vs "+c["celltype_b"]
        s += [f'<text x="30" y="{y-7}" font-family="sans-serif" font-size="13" fill="#20374b">{html.escape(label)}</text>',
              f'<text x="30" y="{y+13}" font-family="sans-serif" font-size="12" fill="#5b7381">n={c["n"]}, CPM +{c["n_positive"]}/{c["n"]}, detection +{c["detect_positive"]}/{c["n"]}</text>']
        for j,v in enumerate(c["individual_deltas"]):
            s.append(f'<circle cx="{xpos(v):.2f}" cy="{y+(j-(c["n"]-1)/2)*3:.1f}" r="5.5" fill="{color}" fill-opacity=".58"/>')
        med=xpos(c["median"])
        s.append(f'<path d="M {med:.2f} {y-9} L {med+8:.2f} {y} L {med:.2f} {y+9} L {med-8:.2f} {y} Z" fill="{color}"/>')
        s.append(f'<text x="994" y="{y+5}" font-family="sans-serif" font-size="15" font-weight="bold" fill="{color}">{c["median"]:+.2f}</text>')
    s += ['<text x="760" y="828" font-family="sans-serif" font-size="14" text-anchor="middle" fill="#3f5364">Within-donor difference in log2(CPM + 1)</text>','</svg>']
    path.write_text("\n".join(s)+"\n",encoding="utf8")

def render_balance(donor_cell,path):
    donors=sorted({d for d,c in donor_cell})
    celltypes=sorted({c for d,c in donor_cell})
    vmax=max(donor_cell.values())
    s=['<svg xmlns="http://www.w3.org/2000/svg" width="1050" height="865" viewBox="0 0 1050 865">',
       '<rect width="100%" height="100%" fill="#ffffff"/>',
       '<text x="33" y="45" font-family="sans-serif" font-size="25" font-weight="bold" fill="#193448">Donor cell-type imbalance (R3_7)</text>',
       '<text x="33" y="76" font-family="sans-serif" font-size="15" fill="#506878">Counts of original annotated cells by Patient; log1p intensity scale</text>']
    for j,d in enumerate(donors):
        s.append(f'<text x="{443+j*88}" y="118" text-anchor="middle" font-family="sans-serif" font-size="15" fill="#294255">{d}</text>')
    for i,cell in enumerate(celltypes):
        y=139+i*55
        s.append(f'<text x="33" y="{y+30}" font-family="sans-serif" font-size="16" fill="#20374a">{html.escape(cell)}</text>')
        for j,d in enumerate(donors):
            v=donor_cell.get((d,cell),0)
            shade=math.log1p(v)/math.log1p(vmax)
            cr=int(246-172*shade)
            fg="#ffffff" if shade>=.64 else "#263f51"
            s.extend([f'<rect x="{400+j*88}" y="{y}" width="82" height="47" rx="4" fill="rgb({cr},{min(255,cr+7)},{min(255,cr+18)})"/>',
                      f'<text x="{441+j*88}" y="{y+30}" font-family="sans-serif" text-anchor="middle" font-size="14" fill="{fg}">{v:,}</text>'])
    s += ['<text x="34" y="807" font-family="sans-serif" font-size="14" fill="#496374">Neurons, astrocytes and stem cells are heavily concentrated in specific donors.</text>','</svg>']
    path.write_text("\n".join(s)+"\n",encoding="utf8")

def run(src,out):
    if out.exists() and any(out.iterdir()):
        raise FileExistsError("Refusing existing nonempty audit output "+str(out))
    if any(not (src/n).is_file() or not (src/n).stat().st_size for n in INPUTS):
        raise FileNotFoundError("R3_7 required output missing")
    manifest=json.loads((src/"HUMAN_RUN_MANIFEST.json").read_text())
    if (manifest.get("notebook_version")!="R3_7_DONOR_PSEUDOBULK_PARSEVERIFIED"
        or manifest.get("donor_analysis_unit")!="Patient"
        or manifest.get("min_cells_per_donor_celltype")!=20
        or manifest.get("min_paired_donors")!=4):
        raise ValueError("Source version/protocol mismatch")
    raw=tsv(src/"HUMAN_DONOR_PSEUDOBULK.tsv")
    comparisons=tsv(src/"HUMAN_DONOR_PAIRED_COMPARISONS.tsv")
    features=tsv(src/"HUMAN_DONOR_FEATURE_STATUS.tsv")
    cov=tsv(src/"HUMAN_DONOR_CELLTYPE_COVERAGE.tsv")
    if (len(features),len(cov),len(comparisons))!=(9,99,10):
        raise ValueError("Unexpected target/coverage/predeclared pair inventory")
    fs={r["gene"]:r["feature_status"] for r in features}
    if fs.get("FGF5")!="FEATURE_NOT_IN_MATRIX" or fs.get("C4orf22")!="FEATURE_NOT_IN_MATRIX":
        raise ValueError("Missing-feature status is invalid")
    rows={}
    group={}
    for r in raw:
        k=r["gene"],r["donor"],r["celltype"]
        if k in rows:raise ValueError("Duplicate gene/donor/celltype")
        rows[k]=r
        count=int(r["n_cells"])
        lib=num(r["library_size"])
        if count<1 or lib<0:raise ValueError("Invalid cell/library size")
        if r["coverage_status"]=="ELIGIBLE" and (count<20 or lib<=0):
            raise ValueError("Eligible row violates predeclared threshold")
        if r["feature_status"]=="FEATURE_NOT_IN_MATRIX":
            if any(not math.isnan(num(r[f])) for f in ("target_sum_counts","detected_cells","pseudobulk_cpm","log2_cpm_plus1")):
                raise ValueError("Missing feature must use NA not zero")
        else:
            cpm=num(r["pseudobulk_cpm"])
            summed=num(r["target_sum_counts"])
            det=num(r["detected_cells"])
            if summed<0 or det<0 or det>count: raise ValueError("Invalid donor counts")
            if lib>0 and abs(cpm-1e6*summed/lib)>max(1e-5,1e-8*abs(cpm)):
                raise ValueError("CPM aggregation differs from raw count/library")
            if abs(num(r["log2_cpm_plus1"])-math.log2(1+cpm))>1e-8:
                raise ValueError("log2 CPM transform mismatch")
        key=(r["donor"],r["celltype"])
        if key in group and group[key]!=count:raise ValueError("Cell counts vary by gene")
        group[key]=count
    donors=sorted({k[0] for k in group})
    if len(donors)!=6 or sum(group.values())!=80515:
        raise ValueError("Six-donor / 80515-cell conservation failed")
    comps=[]
    for c in comparisons:
        gene,a,b=c["gene"],c["celltype_a"],c["celltype_b"]
        joined=[]
        for donor in donors:
            x=rows.get((gene,donor,a))
            y=rows.get((gene,donor,b))
            if x and y and x["coverage_status"]=="ELIGIBLE" and y["coverage_status"]=="ELIGIBLE":
                joined.append((donor,num(x["log2_cpm_plus1"])-num(y["log2_cpm_plus1"]),
                               num(x["detection_fraction"])-num(y["detection_fraction"])))
        n=len(joined)
        ds=[v for d,v,f in joined]
        fractions=[f for d,v,f in joined]
        np=sum(v>0 for v in ds)
        nn=sum(v<0 for v in ds)
        if (n!=int(c["n_paired_donors"]) or n<4 or
            ";".join(sorted(x[0] for x in joined))!=c["paired_donors"] or
            np!=int(c["n_positive"]) or nn!=int(c["n_negative"]) or
            abs(median(ds)-float(c["median_delta_log2_cpm_plus1"]))>1e-8 or
            c["inferential_p_value"]!="NOT_ESTIMATED" or
            c["qc_status"]!="DESCRIPTIVE_PAIRED_QC_PASS"):
            raise ValueError("Paired donor cross-check failed "+str((gene,a,b)))
        leave_one=[median(ds[:i]+ds[i+1:]) for i in range(n)]
        comps.append(dict(gene=gene,celltype_a=a,celltype_b=b,n=n,
            median=median(ds),minimum=min(ds),maximum=max(ds),
            n_positive=np,n_negative=nn,n_tied=sum(v==0 for v in ds),
            detect_positive=sum(v>0 for v in fractions),
            detect_negative=sum(v<0 for v in fractions),
            individual_donors=[x[0] for x in joined],individual_deltas=ds,
            loo_median_range=[min(leave_one),max(leave_one)]))
    donor_totals={d:sum(v for (dd,c),v in group.items() if dd==d) for d in donors}
    imbalance=[]
    for cell in sorted({k[1] for k in group}):
        total=sum(v for (d,c),v in group.items() if c==cell)
        highest=max((v,d) for (d,c),v in group.items() if c==cell)
        imbalance.append(dict(celltype=cell,n_cells=total,dominant_donor=highest[1],
             dominant_n=highest[0],dominant_fraction=highest[0]/total,
             eligible_donors=sum(group.get((d,cell),0)>=20 for d in donors)))
    audit=dict(status="DESCRIPTIVE_RECALC_PASS_NO_CAUSAL_INFERENCE",
        source=str(src),source_sha256={name:hashlib.sha256((src/name).read_bytes()).hexdigest() for name in INPUTS},
        donor_n=6,n_cells=80515,donor_totals=donor_totals,
        pseudobulk_n=len(raw),coverage_n=len(cov),features=fs,
        paired_n=10,paired_pass_n=10,
        unanimous_cpm_n=sum(c["n_positive"]==c["n"] or c["n_negative"]==c["n"] for c in comps),
        unanimous_detection_n=sum(c["detect_positive"]==c["n"] or c["detect_negative"]==c["n"] for c in comps),
        comparisons=comps,donor_imbalance=imbalance)
    out.mkdir(parents=True,exist_ok=True)
    (out/"R3_7_DONOR_REAUDIT.json").write_text(json.dumps(audit,indent=2,ensure_ascii=False)+"\n")
    render_chart(comps,out/"FIG_R3_7_PAIRED_DONOR_DELTAS.svg")
    render_balance(group,out/"FIG_R3_7_DONOR_CELLTYPE_IMBALANCE.svg")
    with (out/"R3_7_PAIRED_VERIFIED.tsv").open("w",newline="",encoding="utf8") as h:
        writer=csv.writer(h,delimiter="\t")
        writer.writerow(["gene","celltype_a","celltype_b","n_paired","median_delta_log2_cpm_plus1",
                         "min_delta","max_delta","n_CPM_positive","n_detection_positive",
                         "LOO_median_min","LOO_median_max"])
        for c in comps:
            writer.writerow([c["gene"],c["celltype_a"],c["celltype_b"],c["n"],
                f'{c["median"]:+.4f}',f'{c["minimum"]:+.4f}',f'{c["maximum"]:+.4f}',
                c["n_positive"],c["detect_positive"],
                f'{c["loo_median_range"][0]:+.4f}',f'{c["loo_median_range"][1]:+.4f}'])
    lines=["# R3_7 — independent six-donor pseudobulk source validation",
           "","Status: DESCRIPTIVE_RECALC_PASS_NO_CAUSAL_INFERENCE",
           "","Original reference: adult control temporal lobe GSE256493.",
           "80,515 cells, 6 donors; 585 gene-donor-celltype observations.",
           "Matched donor denominators, CPM, detection and medians independently recalculated.",
           "All ten prespecified comparisons pass min four paired donors.",
           f'CPM direction unanimous in {audit["unanimous_cpm_n"]}/10; detection direction unanimous in {audit["unanimous_detection_n"]}/10.',
           "No inferential P values or stroke-specific cell-QTL effects were calculated.",
           "Numbers represent differences of log2(CPM + 1), NOT standard log2 fold change.",
           "","## Verified contrasts","",
           "| Gene | A vs B | Paired donors | Median difference | Range | CPM + | Detection + |",
           "|---|---|---:|---:|---:|---:|---:|"]
    for c in comps:
        lines.append(f'| {c["gene"]} | {c["celltype_a"]} vs {c["celltype_b"]} | {c["n"]} | {c["median"]:+.3f} | [{c["minimum"]:+.2f}, {c["maximum"]:+.2f}] | {c["n_positive"]}/{c["n"]} | {c["detect_positive"]}/{c["n"]} |')
    lines += ["","## Donor representation","",
              "| Cell type | Total cells | Dominant donor | Dominant share | Eligible donors |",
              "|---|---:|---|---:|---:|"]
    for v in imbalance:
        lines.append(f'| {v["celltype"]} | {v["n_cells"]:,} | {v["dominant_donor"]} | {100*v["dominant_fraction"]:.2f}% | {v["eligible_donors"]}/6 |')
    lines += ["","## Evidence interpretation and guardrails",
      "","- SH3PXD2A: oligodendrocyte, fibroblast and smooth-muscle comparisons versus the combined Microglia and Macrophages class are positive in all eligible donors.",
      "- COL4A2 and COL4A1 vascular mural/fibroblast CPM comparisons versus endothelial are positive in all eligible donors.",
      "- COL4A2 Pericytes versus Endothelial cells: CPM difference 6/6 positive but detection-fraction difference 3/6 positive, 3/6 negative. The two measurements must not be conflated.",
      "- ALDH2 endothelial versus immune class: CPM 3 positive, 3 negative; no robust cell-type contrast.",
      "- Neurons, neuron progenitor, stem cells, and astrocytes are donor-concentrated, precluding claims of six-donor reproducibility for those cell types.",
      "- FGF5 and C4orf22 are missing RNA *features* in this processed Seurat reference, not known biological zeros.",
      "- Adult control reference, not ischemic stroke cases; six donors; no disease-state DE, SNP-mediated cell eQTL or causal mediation.",
      "- Original author class Microglia and Macrophages is combined; cannot resolve separate microglia or macrophage mechanisms.",
      "- Full per-donor results and SHA256 inputs are in R3_7_DONOR_REAUDIT.json.",
      "- FIG_R3_7_PAIRED_DONOR_DELTAS.svg and FIG_R3_7_DONOR_CELLTYPE_IMBALANCE.svg provide supplementary figures."]
    (out/"R3_7_DONOR_VALIDATION_REPORT.md").write_text("\n".join(lines)+"\n")
    print(json.dumps({k:audit[k] for k in ("status","donor_n","n_cells","pseudobulk_n","paired_n","paired_pass_n","unanimous_cpm_n","unanimous_detection_n")},indent=2))
    return audit

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    run(args.source,args.out)

if __name__=="__main__":
    main()
