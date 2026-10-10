#!/usr/bin/env python3
"""Explicit GRCh37 -> GRCh38 UCSC chain SNP overlap for existing EUR region panels.

Position chain only; true REF-base verification requires reference FASTA.
No study-cohort GTEx LD/causal colocalization claim.
"""
import argparse,bisect,csv,gzip,json,hashlib
from pathlib import Path
from collections import defaultdict
def parse_chain(path):
    out=defaultdict(list)
    with gzip.open(path,"rt") as f:
        lines=iter(f)
        for line in lines:
            if not line.startswith("chain "):continue
            h=line.split()
            if len(h)!=13:raise ValueError("Malformed chain")
            score=int(h[1]);tc=h[2];ts=h[4];qb=h[7];qs=h[9]
            t=int(h[5]);q=int(h[10]);blocks=[]
            for l in lines:
                if not l.strip():break
                x=list(map(int,l.split()))
                if len(x) not in (1,3):raise ValueError("Bad chain block")
                n=x[0]
                blocks.append((t,t+n,q,q+n))
                t+=n;q+=n
                if len(x)==3:t+=x[1];q+=x[2]
            if tc==qb and ts=="+" and qs=="+" and tc.startswith("chr") and tc[3:].isdigit():
                out[tc[3:]].append((score,blocks))
    indexed={}
    for chrom,choices in out.items():
        choices.sort(key=lambda x:-x[0])
        if len(choices)>1 and choices[0][0]==choices[1][0]:raise ValueError("Ambiguous chain top score")
        blocks=choices[0][1]
        indexed[chrom]=(blocks,[x[0] for x in blocks])
    return indexed
def lift(ref,chrom,pos):
    z=ref.get(chrom.replace("chr",""))
    if not z:return None
    b,starts=z
    idx=bisect.bisect_right(starts,int(pos)-1)-1
    if idx<0 or int(pos)-1>=b[idx][1]:return None
    return b[idx][2]+int(pos)-1-b[idx][0]+1
def read_tsv(path):
    with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def panel(path,mapper):
    pvar=path.with_suffix(".pvar")
    if not pvar.exists():raise ValueError("Missing pvar "+str(path))
    sites=set();cnt=0;nonmapped=0
    with pvar.open() as f:
        fields=None
        for line in f:
            if line.startswith("##"):continue
            if line.startswith("#"):
                fields=line.strip().lstrip("#").split("\t")
                continue
            if fields is None:raise ValueError("Missing pvar header")
            r=dict(zip(fields,line.strip().split("\t")))
            cnt+=1
            mapped=lift(mapper,r["CHROM"],r["POS"])
            if mapped is None:nonmapped+=1;continue
            sites.add((r["CHROM"].replace("chr",""),str(mapped),r["REF"].upper(),r["ALT"].upper()))
    return sites,cnt,nonmapped
def run(args):
    mapper=parse_chain(args.chain)
    allpanels=[]
    for p in sorted(args.eur.rglob("*.pgen")):
        if not p.with_suffix(".pvar").exists():continue
        allpanels.append((p,*panel(p,mapper)))
    if not allpanels:raise ValueError("No existing EUR reference panels")
    targets=[r for r in read_tsv(args.replay) if r["status"]=="PASS"]
    if len(targets)!=30:raise ValueError("Original PASS assay count drift")
    result=[]
    for x in targets:
        path=args.source/x["filename"]
        snps=read_tsv(path)
        q={(r["chromosome"].replace("chr",""),r["position"],r["ref"].upper(),r["alt"].upper()) for r in snps}
        ranked=sorted(((len(q&sites),str(p),n,unm) for p,sites,n,unm in allpanels),reverse=True)
        n,best,n_panel,nonmapped=ranked[0]
        result.append(dict(locus=x["locus"],dataset_key=x["dataset_key"],gene_base=x["gene_base"],
             source_SNPs=len(snps),best_EUR_reference_b37_to_b38_exact_allele_matches=n,
             percent_of_assay_SNPs=round(100*n/len(snps),3),best_panel=best,
             reference_panel_variants=n_panel,chain_unmapped_panel_SNPs=nonmapped,
             true_grch38_reference_fasta_verified=False,
             genotype_samples_exactly_GTex_donors=False,validated_multi_signal_coloc=False))
    args.out.mkdir(exist_ok=True,parents=True)
    with (args.out/"IS_30_EUR_REFERENCE_B37_TO_B38_LIFTED_SNP_OVERLAP.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter="\t");w.writeheader();w.writerows(result)
    summary=dict(audited_assays=len(result),eur_regional_panels=len(allpanels),
       assays_any_exact_allele_lifted_match=sum(r["best_EUR_reference_b37_to_b38_exact_allele_matches"]>0 for r in result),
       assays_over_50pct_overlap=sum(r["percent_of_assay_SNPs"]>=50 for r in result),
       maximum_matched_snps=max(r["best_EUR_reference_b37_to_b38_exact_allele_matches"] for r in result),
       chain_sha256=hashlib.sha256(args.chain.read_bytes()).hexdigest(),
       warning="Chain coordinate match plus REF/ALT string match is not GRCh38 FASTA or study-cohort GTEx LD verification",
       causal_gene_confirmed=False)
    (args.out/"IS_30_EUR_B37_LIFTOVER_AUDIT.json").write_text(json.dumps(summary,indent=2)+"\n")
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eur","chain","replay","source","out"):p.add_argument("--"+k,required=True,type=Path)
    print(json.dumps(run(p.parse_args()),indent=2))
