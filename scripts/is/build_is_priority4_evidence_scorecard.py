#!/usr/bin/env python3
"""Integrate IS P0.1/P0.2/P0.3 into a fail-closed evidence scorecard.

This is an evidence *status* ledger, not a causal-ranking or posterior calculator.
No unexecuted pipeline phase is promoted to complete.
"""
import argparse,csv,json,hashlib
from pathlib import Path

GENES=("CALHM2","NEURL","FGF5","C4orf22")
LOCUS={"CALHM2":"BBJ_IS_L002","NEURL":"BBJ_IS_L002",
       "FGF5":"BBJ_IS_L001","C4orf22":"BBJ_IS_L001"}
EXPECTED={"CALHM2":3,"NEURL":0,"FGF5":5,"C4orf22":0}
def read(path):
    with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()
def integrate(comparison,ledger,pairwise,out):
    out=Path(out)
    if out.exists():raise FileExistsError("Refusing to overwrite a prior audit")
    sources={"P03_gene_comparator":Path(comparison),
             "P01_CALHM2_3SNP_ledger":Path(ledger),
             "P02_CALHM2_pairwise_LD":Path(pairwise)}
    existence={name:path.is_file() for name,path in sources.items()}
    # P03 is needed to generate a gene-level file; missing P01/P02 are explicitly PENDING.
    if not existence["P03_gene_comparator"]:
        raise FileNotFoundError("P03 comparator must be executed first; no fake scorecard")
    rows=read(sources["P03_gene_comparator"])
    lookup={r["gene"]:r for r in rows}
    if set(lookup)!=set(GENES) or len(rows)!=4:
        raise ValueError("Gene scope mismatch or duplicated gene")
    for gene in GENES:
        r=lookup[gene]
        if r["locus"]!=LOCUS[gene]:raise ValueError("Gene-locus mismatch")
        if int(r["QTD_nominal_p_lt_1e4"])!=EXPECTED[gene]:
            raise ValueError("QTL source cardinality drift")
    nledger=None
    if existence["P01_CALHM2_3SNP_ledger"]:
        variants=read(sources["P01_CALHM2_3SNP_ledger"])
        if len(variants)!=3 or {r["gene"] for r in variants}!={"CALHM2"}:
            raise ValueError("P01 SNP ledger unexpected")
        if len({r["SNP_GRCh37"] for r in variants})!=3:
            raise ValueError("P01 duplicated genomic SNP")
        nledger=3
    npairs=None;cluster=None
    if existence["P02_CALHM2_pairwise_LD"]:
        pairs=read(sources["P02_CALHM2_pairwise_LD"])
        if len(pairs)!=30:raise ValueError("P02 3 populations x 10 pairs expected")
        population={"EAS504":0,"JPT104":0,"otherEAS400":0}
        for r in pairs:
            if r["population"] not in population:raise ValueError("Unknown population")
            population[r["population"]]+=1
        if any(n!=10 for n in population.values()):raise ValueError("Invalid pair matrix")
        npairs=len(pairs)
        if nledger is not None:
            snps={r["SNP_GRCh37"] for r in variants}
            cluster={}
            for pop in population:
                subset=[r for r in pairs if r["population"]==pop and
                    r["variant_1"] in snps and r["variant_2"] in snps]
                if len(subset)!=3:raise ValueError("Missing molecular SNP pair")
                cluster[pop]=all(float(r["r2"])>=.8 for r in subset)
    consolidated=[]
    for gene in GENES:
        r=lookup[gene]
        consolidated.append(dict(
            gene=gene,locus=LOCUS[gene],
            QTL_nominal_p_lt_1e4=int(r["QTD_nominal_p_lt_1e4"]),
            allele_exact_GWAS_SNPs=int(r["both_GWAS_allele_exact_nonpal_unique"]),
            same_SNP_BBJ_and_GIGA_GWS=int(r["BBJ_and_GIGA_GWS_same_SNP"]),
            triple_nominal_GWAS_overlap=int(r["triple_QTD_BBJ_GIGA_threshold"]),
            P03_status="VERIFIED_INPUT_PRESENT",
            P01_status=("VERIFIED_INPUT_PRESENT" if nledger==3 else "NOT_EXECUTED_OR_MISSING") if gene=="CALHM2" else "NOT_APPLICABLE",
            P02_status=("VERIFIED_INPUT_PRESENT" if npairs==30 else "NOT_EXECUTED_OR_MISSING") if gene=="CALHM2" else "NOT_APPLICABLE",
            EAS_three_QTL_SNPs_one_high_LD_cluster=(str(cluster["EAS504"]) if cluster is not None and gene=="CALHM2" else "NOT_EVALUATED"),
            JPT_three_QTL_SNPs_one_high_LD_cluster=(str(cluster["JPT104"]) if cluster is not None and gene=="CALHM2" else "NOT_EVALUATED"),
            evidence_level="DESCRIPTIVE_QTL_GWAS_PUBLIC_REFERENCE_ONLY",
            causal_gene_status="NOT_ESTABLISHED",
            formal_multisignal_coloc_status="BLOCKED_QTL_STUDY_LD_AND_N"))
    out.mkdir(parents=True)
    with (out/"IS_PRIORITY4_INTEGRATED_EVIDENCE_SCORECARD.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(consolidated[0]),delimiter="\t")
        w.writeheader();w.writerows(consolidated)
    manifest={"state":"INTEGRATED_EVIDENCE_NO_CAUSAL_RANKING",
        "sources":{name:{"path":str(path),"exists":existence[name],
          "sha256":sha(path) if existence[name] else None} for name,path in sources.items()},
        "p01_individual_calhm2_snp_count":nledger,
        "p02_valid_pair_rows":npairs,
        "p03_gene_count":4,
        "no_formal_coloc":True,"no_gwas_conditional_analysis":True,
        "original_2225_candidates_preserved":True}
    (out/"IS_PRIORITY4_INTEGRATED_EVIDENCE_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return consolidated

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--comparison",type=Path,required=True)
    p.add_argument("--ledger",type=Path,required=True)
    p.add_argument("--pairwise",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(integrate(a.comparison,a.ledger,a.pairwise,a.out),indent=2))
