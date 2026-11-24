#!/usr/bin/env python3
"""Read pinned published aggregate RData without installing the full upstream package."""
import argparse,hashlib,json,re,urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
import pyreadr

SOURCE="https://raw.githubusercontent.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/70461ea4a55716bb00f8467860521e8a8f94b9f9/data/"
PINS={"MUSCLE_TRNSCRPT_DA":(22711273,"0a4168c61102863487a4c4ecb95ec66bc250329e"),
      "HUMAN_FEATURE_TO_GENE":(20001757,"db4dd772a08ed10b1063201179386a2c113ecb55")}
SHORT="Resist.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)"
def require(ok,msg):
    if not ok:raise ValueError(msg)
def load(root,name):
    size,sha=PINS[name];p=root/(name+".rda");url=SOURCE+p.name
    req=urllib.request.Request(url,headers={"User-Agent":"Public-muscle-signature-audit/1.0"})
    with urllib.request.urlopen(req,timeout=45) as f,p.open("wb") as w:
        n=0
        while b:=f.read(1024*1024):
            n+=len(b);require(n<=size,"Publisher object size exceeded");w.write(b)
    b=p.read_bytes()
    require(len(b)==size and hashlib.sha1(("blob "+str(size)).encode()+bytes([0])+b).hexdigest()==sha,"Exact Git blob mismatch")
    h=hashlib.sha256(b).hexdigest();del b
    d=pyreadr.read_r(str(p))
    require(set(d)=={name},f"R object identity {list(d)}")
    x=d[name];require(isinstance(x,pd.DataFrame),"Expected public data.frame")
    print("SOURCE",name,len(x),list(x.columns),h,flush=True)
    return x,{"url":url,"bytes":size,"git_blob":sha,"sha256":h}
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.root.mkdir(parents=True,exist_ok=True);a.out.mkdir(parents=True,exist_ok=True)
    da,ds=load(a.root,"MUSCLE_TRNSCRPT_DA")
    cols=["tissue","assay","full_model","contrast","contrast_short","contrast_type","contrast_category","feature_id","logFC","p_value","adj_p_value"]
    require(set(cols)<=set(da),"DA schema")
    print("SOURCE_CONTRASTS",da["contrast_short"].value_counts().to_dict(),flush=True)
    mask=(da.tissue=="muscle")&(da.assay=="transcript-rna-seq")&(da.contrast_type=="exercise_with_controls")&(da.contrast_category=="RE-CON")&(da.contrast_short==SHORT)
    x=da.loc[mask,cols].copy();del da
    require(len(x)>1000 and x.full_model.nunique()==1 and x.contrast.nunique()==1 and not x.feature_id.duplicated().any(),"Fixed contrast/model identity")
    x["feature_id"]=x.feature_id.astype(str)
    total_rows=len(x)
    standard=x.feature_id.str.fullmatch(r"ENSG[0-9]+(?:\.[0-9]+)?")
    finite=np.isfinite(x[["logFC","p_value","adj_p_value"]].to_numpy(dtype=float)).all(axis=1)
    require((x.loc[finite,"p_value"].between(0,1)&x.loc[finite,"adj_p_value"].between(0,1)).all(),"P-value bounds")
    source_sig=finite & (x.adj_p_value<.05) & (x.logFC!=0)
    source_fdr_up=int((source_sig&(x.logFC>0)).sum());source_fdr_down=int((source_sig&(x.logFC<0)).sum())
    model=str(x.full_model.iloc[0]);contrast=str(x.contrast.iloc[0])
    print("FEATURE_ID_GATE",{"total":total_rows,"canonical":int(standard.sum()),"noncanonical":int((~standard).sum()),
          "noncanonical_examples":x.loc[~standard,"feature_id"].head(8).tolist(),"nonfinite":int((~finite).sum())},flush=True)
    # Exclude noncanonical/unestimable source rows from the mapping pool, without changing
    # source FDR family or threshold and before reading any validation expression.
    require(int(standard.sum())>=.9*total_rows and int(finite.sum())>=.9*total_rows,"Too few canonical/estimable features")
    id_gate={"canonical":int(standard.sum()),"noncanonical":int((~standard).sum()),"nonfinite":int((~finite).sum())}
    x=x.loc[standard&finite].copy()
    m,ms=load(a.root,"HUMAN_FEATURE_TO_GENE")
    mc=["assay","feature_id","entrez_gene","ensembl_gene","gene_symbol"]
    require(set(mc)<=set(m),"Feature map schema")
    y=m.loc[(m.assay=="transcript-rna-seq")&m.feature_id.astype(str).isin(x.feature_id),mc[1:]].copy();del m
    y=y.astype(str).drop_duplicates()
    y["ensembl_gene"]=y.ensembl_gene.str.replace(r"\.[0-9]+$","",regex=True)
    y=y.loc[y.entrez_gene.str.fullmatch(r"[0-9]+")&y.ensembl_gene.str.fullmatch(r"ENSG[0-9]+")&(y.ensembl_gene==y.feature_id.str.replace(r"\.[0-9]+$","",regex=True))]
    y=y.loc[~y.feature_id.duplicated(keep=False)]
    y=y.loc[~y.entrez_gene.duplicated(keep=False)&~y.ensembl_gene.duplicated(keep=False)]
    x=x.merge(y,how="left",on="feature_id",validate="one_to_one")
    sig=(x.adj_p_value<.05)&(x.logFC!=0)
    choices=x.loc[sig&x.entrez_gene.notna()].copy()
    choices["absFC"]=choices.logFC.abs()
    parts=[]
    for sign in (1,-1):
        v=choices.loc[choices.logFC*sign>0].sort_values(["p_value","absFC","ensembl_gene"],ascending=[True,False,True]).head(100).copy()
        require(len(v)>=20,"Fewer than 20 source significant unique genes in one direction")
        v["sign"]=sign;parts.append(v)
    panel=pd.concat(parts,ignore_index=True).rename(columns={"entrez_gene":"entrez","ensembl_gene":"ensembl","gene_symbol":"symbol"})
    require(not panel.entrez.duplicated().any() and not panel.ensembl.duplicated().any(),"Unique panel identity")
    output=a.out/"panel.tsv";panel[["entrez","ensembl","symbol","sign","logFC","p_value","adj_p_value"]].to_csv(output,sep="\t",index=False,na_rep="")
    summary={"status":"FROZEN_SOURCE_ONLY_NO_VALIDATION_OUTCOME","publisher_package":"2.0.8","source_collection":"c2.0","contrast_short":SHORT,
       "model":model,"contrast":contrast,"source_rows":total_rows,"feature_id_gate":id_gate,
       "source_fdr_up":source_fdr_up,"source_fdr_down":source_fdr_down,
       "unique_mapped":int(x.entrez_gene.notna().sum()),"eligible_fdr_up":int((choices.logFC>0).sum()),"eligible_fdr_down":int((choices.logFC<0).sum()),
       "panel_up":int((panel.sign==1).sum()),"panel_down":int((panel.sign==-1).sum()),
       "panel_sha256":hashlib.sha256(output.read_bytes()).hexdigest(),"sources":[ds,ms]}
    (a.out/"source.json").write_text(json.dumps(summary,indent=2)+"\n");print(json.dumps(summary,indent=2),flush=True)
if __name__=="__main__":main()
