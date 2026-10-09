#!/usr/bin/env python3
"""EBI eQTL Catalogue GTEx Artery Aorta nominal ALL-SNP QTL -> EAS AIS G0022.

The source is a *region-scoped* GRCh38 tabix extract, not full GTEx dataset.
Requires exact, unique hg38->hg19 allele conversion, no inferred strand flips.
Prepares coloc.abf INPUTS ONLY, not publication-ready causal genes.
"""
import argparse,csv,gzip,json,math,sys,hashlib
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_g0022_gtex_v8_dapg import variant_liftover,stable,VENDOR

GWAS=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz")
CHAIN=Path("/srv/is-analysis/data/reference/hg38ToHg19.over.chain")
QTL=Path("/srv/is-analysis/data/is/qtl/eqtl_catalogue/G0022_v1/QTD000131_artery_aorta_GRCh38_chr12_111650000_112370000.tsv.gz")
OUT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/eqtl_catalogue_aorta_v1")
ANNOTATION=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/IS_ALL_GENE_WINDOW_UNIVERSE.tsv")
EQTLCATALOGUE_COLUMNS=[
 "molecular_trait_id","chromosome","position","ref","alt","variant",
 "ma_samples","maf","pvalue","beta","se","type","ac","an","r2",
 "molecular_trait_object_id","gene_id","median_tpm","rsid"]
GENES=("ACAD10","ALDH2","NAA25","HECTD4","BRAP","PTPN11","IFT81","ATP2A2")
def gwas(path):
    vals={}
    with gzip.open(path,"rt") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            if r["group_id"]!="IS_XDATA_G0022" or r["dataset"]!="GCST90104545" or r["qc_status"] not in ("MATCH_ALT_EFFECT","MATCH_REF_EFFECT"):
                continue
            k=r["variant_id"]
            if k in vals:raise ValueError("GWAS duplicate variant: "+k)
            vals[k]=r
    if len(vals)<4500:raise ValueError("GWAS reference-aligned locus not available")
    return vals
def numeric(text):
    try:
        x=float(text)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None
def candidate_genes(path):
    with path.open() as f:
        rows=[r for r in csv.DictReader(f,delimiter="\t") if r["group_id"]=="IS_XDATA_G0022" and r["gene_symbol"] in GENES]
    if {r["gene_symbol"] for r in rows}!=set(GENES):raise ValueError("GENCODE candidates missing")
    return {stable(r["gene_id"]):r["gene_symbol"] for r in rows}
def parse_nominal(line):
    cols=line.rstrip("\n").split("\t")
    if len(cols)!=19:raise ValueError("Unexpected 19-column eQTL Catalogue nominal schema")
    return dict(zip(EQTLCATALOGUE_COLUMNS,cols))
def read_qtl(path,allele_converter,gene_index,source_gwas):
    qc=Counter();matched=defaultdict(dict)
    with gzip.open(path,"rt") as h:
        for line in h:
            if not line.strip():continue
            r=parse_nominal(line)
            qc["source_records"]+=1
            geneid=stable(r["gene_id"])
            if geneid not in gene_index:continue
            qc["target_genes_records"]+=1
            if r["chromosome"]!="12" or r["type"]!="SNP" or len(r["ref"])!=1 or len(r["alt"])!=1:
                qc["not_snp"]+=1;continue
            if str(r["variant"])!=f'chr12_{r["position"]}_{r["ref"]}_{r["alt"]}':
                raise ValueError("QTL variant tuple disagreement")
            key,status=variant_liftover(r["variant"]+"_b38",allele_converter)
            if key is None:
                qc[status]+=1;continue
            gwas_id,strand=key
            qtl_beta=numeric(r["beta"]);qtl_se=numeric(r["se"]);qtl_p=numeric(r["pvalue"])
            qtl_maf=numeric(r["maf"])
            if (qtl_beta is None or qtl_se is None or qtl_se<=0
                or qtl_p is None or not 0<=qtl_p<=1
                or qtl_maf is None or not 0<qtl_maf<=.5):
                qc["bad_qtl_effect"]+=1;continue
            g=source_gwas.get(gwas_id)
            if not g:
                qc["no_gwas_allele_pair"]+=1;continue
            if {g["ref"],g["alt"]} in ({"A","T"},{"C","G"}):
                qc["palindrome"]+=1;continue
            if r["ref"]==r["alt"]:
                qc["same_ref_alt"]+=1;continue
            ref_check=gwas_id.split(":")
            if ref_check[2]!=g["ref"] or ref_check[3]!=g["alt"]:
                raise ValueError("GWAS genomic reference mismatch")
            keygene=(geneid,gwas_id)
            record={
                "gene_symbol":gene_index[geneid],"gene_id":geneid,
                "variant_id":gwas_id,"qtl_grch38_variant":r["variant"],
                "liftover_strand":strand,"gwas_beta":g["alt_effect_beta"],
                "gwas_se":g["se"],"gwas_p":g["p"],"gwas_alt_eaf":g["alt_effect_eaf"],
                "qtl_beta":r["beta"],"qtl_se":r["se"],"qtl_p":r["pvalue"],"qtl_maf":r["maf"],
                "qtl_an":r["an"],"qtl_dataset_id":"QTD000131",
                "study_tissue":"GTEx_Artery_Aorta","qtl_assembly":"GRCh38",
                "gwas_assembly":"GRCh37","allele_match":"EXACT_REF_ALT_AFTER_LIFTOVER"}
            prev=matched[geneid].get(gwas_id)
            if prev is not None:
                for x in ("qtl_beta","qtl_se","qtl_p","qtl_maf"):
                    if prev[x]!=record[x]:
                        raise RuntimeError("Duplicate QTL SNP with different effects for "+str(keygene))
                qc["same_variant_multiple_rsid_duplicate"]+=1
                continue
            matched[geneid][gwas_id]=record
            qc["matched_unique_gene_snp"]+=1
    return matched,qc
def prepare(args):
    if not args.qtl.is_file():raise FileNotFoundError("Full-source GTEx artery aorta region required")
    from audit_g0022_gtex_v8_dapg import LiftOver
    if LiftOver is None:
        raise RuntimeError("pyliftover required for eQTL execution")
    converter=LiftOver(str(args.chain))
    src=gwas(args.gwas)
    gene=candidate_genes(args.annotation)
    matched,qc=read_qtl(args.qtl,converter,gene,src)
    args.out.mkdir(parents=True,exist_ok=True)
    audit=[]
    for gene_id,gene_symbol in sorted(gene.items(),key=lambda x:x[1]):
        group=matched.get(gene_id,{})
        rows=sorted(group.values(),key=lambda r:int(r["variant_id"].split(":")[1]))
        best_qtl=min((float(r["qtl_p"]) for r in rows),default=None)
        best_gwas=min((float(r["gwas_p"]) for r in rows),default=None)
        status=("READY_EXPLORATORY_ABF_COLOC" if
            len(rows)>=100 and best_qtl is not None and best_qtl<=1e-4 and
            best_gwas is not None and best_gwas<=5e-8
            else "BLOCKED_INSUFFICIENT_COMMON_GWAS_QTL_OR_SIGNAL")
        audit.append({"gene":gene_symbol,"gene_id":gene_id,
            "shared_allele_matched_unique_snps":len(rows),
            "minimum_qtl_p":best_qtl if best_qtl is not None else "",
            "minimum_ais_p":best_gwas if best_gwas is not None else "",
            "status":status,"source":"GTEx_eQTL_Catalogue_QTD000131_Aorta",
            "warning":"REGION_TRUNCATED_SINGLE_CAUSAL_ABF_ONLY_NOT_VALIDATED_COLOCALIZATION"})
        if not rows:continue
        with (args.out/f"G0022_{gene_symbol}_AIS_AORTA_ABF_INPUT.tsv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
            w.writeheader();w.writerows(rows)
    with (args.out/"G0022_EQTL_CATALOGUE_AORTA_INPUT_AUDIT.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(audit[0]),delimiter="\t")
        w.writeheader();w.writerows(audit)
    summary={"dataset_id":"QTD000131",
      "source_url":"https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000131/QTD000131.all.tsv.gz",
      "region_GRCh38":"12:111650000-112370000",
      "tissue":"GTEx_Artery_Aorta","study_n":387,
      "GRCh37_reconverted_and_exact_allele_matched":True,
      "gene_tests":len(audit),
      "gene_inputs_with_any_overlap":sum(x["shared_allele_matched_unique_snps"]>0 for x in audit),
      "abf_eligible_gene_inputs":sum(x["status"]=="READY_EXPLORATORY_ABF_COLOC" for x in audit),
      "qtl_qc":dict(qc),
      "qtl_region_extract_sha256":hashlib.sha256(args.qtl.read_bytes()).hexdigest(),
      "source_file_bytes":args.qtl.stat().st_size,
      "warning":"This is a region-only nominal-QTL extract; coloc single-signal assumption and study ancestry mismatch; no validated causal genes"}
    (args.out/"G0022_EQTL_CATALOGUE_AORTA_INPUT_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    for row in audit:print(row)
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--qtl",type=Path,default=QTL)
    p.add_argument("--gwas",type=Path,default=GWAS)
    p.add_argument("--chain",type=Path,default=CHAIN)
    p.add_argument("--annotation",type=Path,default=ANNOTATION)
    p.add_argument("--out",type=Path,default=OUT)
    args=p.parse_args()
    prepare(args)
