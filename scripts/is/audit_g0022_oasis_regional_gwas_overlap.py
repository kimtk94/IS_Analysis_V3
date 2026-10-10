#!/usr/bin/env python3
"""Read-only regional overlap: GRCh37 GIGASTROKE EAS AIS vs OASIS Bokeh GRCh38.

The UCSC highest-score chr12-to-chr12 positive-strand chain is used as a
single-valued *coordinate* reference. Genome REF/ALT is then matched exactly;
non-matches, positions in chain gaps and missing plot SNPs are kept as
distinct QC classes. QTL effect allele orientation is NOT established, and
nothing here constitutes coloc, signal fine-mapping or MR.
"""
import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

from audit_g0022_oasis_browser_qtl import parse
from audit_g0022_oasis_gwas_coverage_and_alleles import CS, GENES, CELLS, AIS_STUDY, AIS_TRAIT, split_variant


def best_primary_chr12_chain(chain):
    """Select a single high-score forward chr12 primary alignment, no alts."""
    candidates=[]
    with gzip.open(chain,"rt") as f:
        lines=iter(f)
        for line in lines:
            if not line.startswith("chain "):
                continue
            head=line.strip().split()
            if len(head)!=13:
                raise ValueError("Invalid UCSC chain header")
            score=int(head[1]);tname,tsize,tstrand,tbeg,tend=head[2:7]
            qname,qsize,qstrand,qbeg,qend=head[7:12]
            ident=head[12]
            blocks=[]
            tcur=int(tbeg);qcur=int(qbeg)
            for block in lines:
                if not block.strip():
                    break
                fields=[int(v) for v in block.split()]
                if len(fields) not in (1,3):
                    raise ValueError("Bad chain block")
                size=fields[0]
                if size<=0:
                    raise ValueError("Negative or empty chain block")
                blocks.append((tcur,tcur+size,qcur,qcur+size))
                tcur+=size;qcur+=size
                if len(fields)==3:
                    dt,dq=fields[1:]
                    if dt<0 or dq<0:
                        raise ValueError("Negative chain gap")
                    tcur+=dt;qcur+=dq
            if tcur!=int(tend) or qcur!=int(qend):
                raise ValueError("Chain end does not match block coordinates")
            if tname=="chr12" and qname=="chr12" and tstrand=="+" and qstrand=="+":
                candidates.append((score,ident,blocks))
    if not candidates:
        raise ValueError("No primary chr12 forward chain")
    ordered=sorted(candidates,key=lambda x:x[0],reverse=True)
    if len(ordered)>1 and ordered[0][0]==ordered[1][0]:
        raise ValueError("Highest-score chr12 chain not unique")
    score,ident,blocks=ordered[0]
    return {"score":score, "id":ident, "blocks":blocks,
            "chain_chromosome":"chr12", "chain_target_strand":"+",
            "strategy":"HIGHEST_SCORE_CHR12_PRIMARY_FORWARD_ONLY"}


def make_mapper(best):
    blocks=best["blocks"]
    if not blocks or blocks != sorted(blocks):
        raise ValueError("Unsorted primary chain blocks")
    starts=[block[0] for block in blocks]
    def convert(gr37_pos_1based):
        p=gr37_pos_1based-1
        if p<0: raise ValueError("Invalid one-based input position")
        idx=bisect.bisect_right(starts,p)-1
        if idx<0 or p>=blocks[idx][1]:
            return None
        return blocks[idx][2]+(p-blocks[idx][0])+1
    return convert


def read_ais_summary(source):
    rows={}
    total=0;skipped_not_QC=0
    with gzip.open(source,"rt",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        required={"dataset","group_id","variant_id","ref","alt","chr","pos",
                  "alt_effect_beta","alt_effect_eaf","se","p","qc_status","build","ancestry"}
        if not required.issubset(reader.fieldnames):
            raise ValueError("Incorrect original GWAS source schema")
        for item in reader:
            if item["group_id"]!="IS_XDATA_G0022" or item["dataset"]!=AIS_STUDY:
                continue
            total+=1
            c,p,ref,alt=split_variant(item["variant_id"],"GRCh37")
            if (c,p,ref,alt)!=(item["chr"],int(item["pos"]),item["ref"],item["alt"]):
                raise ValueError("AIS source variant-to-column mismatch")
            if item["build"]!="GRCh37" or item["ancestry"]!="EAS":
                raise ValueError("Invalid GWAS build or ancestry")
            if item["qc_status"] not in ("MATCH_ALT_EFFECT","MATCH_REF_EFFECT"):
                skipped_not_QC+=1
                continue
            if item["variant_id"] in rows:
                raise ValueError("Duplicate EAS AIS source variant")
            beta,se,pval,eaf=[float(item[v]) for v in
                              ("alt_effect_beta","se","p","alt_effect_eaf")]
            if not (math.isfinite(beta) and math.isfinite(se) and se>0
                    and math.isfinite(pval) and 0<=pval<=1
                    and math.isfinite(eaf) and 0<eaf<1):
                raise ValueError("Invalid AIS statistics")
            rows[item["variant_id"]]={
                "gwas_variant_grch37":item["variant_id"],"chromosome":c,
                "pos_grch37":p,"ref_grch37":ref,"alt_grch37":alt,
                "gwas_ALT_beta":beta,"gwas_se":se,"gwas_p":pval,
                "gwas_ALT_EAF":eaf, "gwas_source_QC":item["qc_status"],
            }
    if total<100 or len(rows)<100:
        raise ValueError("Suspiciously small region variant source")
    return rows,{"ais_source_rows":total,"ais_source_qc_variants":len(rows),
                 "ais_source_skipped_qc":skipped_not_QC}


def map_ais_variants(rows, lifter):
    coord_map={}
    status=Counter()
    for variant,item in rows.items():
        lifted=lifter(item["pos_grch37"])
        if lifted is None:
            status["NO_PRIMARY_CHAIN_MAPPING"]+=1
            continue
        if lifted in coord_map:
            raise ValueError("Two source variants mapped to same coordinate with separate records")
        coord_map[lifted]=dict(item,chromosome_grch38="chr12",pos_grch38=lifted)
        status["PRIMARY_CHAIN_MAPPED"]+=1
    return coord_map,status


def compare_one_page(gene,cell,html,coord_map):
    """Report all site-level positions; do not treat allele mismatch as a valid join."""
    try:
        qtl=parse(html.decode("utf-8","replace"),gene)
        page_state="BOKEH_PLOTTED_VARIANTS_AVAILABLE"
    except ValueError as e:
        if "Missing embedded Bokeh data" not in str(e):
            raise
        qtl=[]
        page_state="NOT_ASSESSED_NO_EMBEDDED_BOKEH_ASSOCIATIONS"
    seen={}
    accepted=[];qc=Counter()
    browser_af_greater_half=0
    for q in qtl:
        q_id=q["variant_id_hg38"]
        c,p,ref,alt=split_variant(q_id,"GRCh38")
        if c!="12":
            raise ValueError("QTL plot includes off-chromosome data")
        if (p,ref,alt) in seen:
            raise ValueError("Duplicated QTL variant allele in one source page")
        seen[(p,ref,alt)]=True
        af_in_bokeh=float(q["maf"])
        if not math.isfinite(af_in_bokeh) or not 0<=af_in_bokeh<=1:
            raise ValueError("Invalid original Bokeh allele frequency")
        browser_af_greater_half+=af_in_bokeh>.5
        original=coord_map.get(p)
        if original is None:
            qc["QTL_SNP_NOT_IN_THIS_GWAS_REGION_SOURCE"]+=1
            continue
        if (ref,alt)!=(original["ref_grch37"],original["alt_grch37"]):
            qc["QTL_AND_GWAS_ALLELE_MISMATCH_REJECTED"]+=1
            continue
        if len(ref)!=1 or len(alt)!=1:
            qc["INDEL_OR_MULTIBASE_VARIANT_EXCLUDED"]+=1
            continue
        beta=float(q["effect_size"]); se=float(q["effect_size_SE"])
        maf=float(q["maf"]);pval=float(q["pval_nominal"])
        if not (math.isfinite(beta) and math.isfinite(se) and se>0 and
                math.isfinite(maf) and 0<=maf<=1 and
                math.isfinite(pval) and 0<=pval<=1):
            raise ValueError("Invalid QTL per-SNP association source")
        accepted.append({
            "gene":gene,"cell":cell,
            "gwas_variant_grch37":original["gwas_variant_grch37"],
            "oasis_variant_grch38":q_id,"oasis_rsid":q["rsid"],
            "source_grch37_pos":original["pos_grch37"],
            "mapped_grch38_pos":p,
            "gwas_ALT_beta":original["gwas_ALT_beta"],
            "gwas_se":original["gwas_se"],
            "gwas_p":original["gwas_p"],
            "gwas_ALT_EAF":original["gwas_ALT_EAF"],
            "oasis_beta_UNHARMONIZED":beta,
            "oasis_se":se,"oasis_p":pval,
            "oasis_browser_maf_column_AS_REPORTED_NOT_TRUE_MAF":maf,
            "source_ref_alt_id_match":"EXACT_AFTER_PRIMARY_UCSC_CHAIN",
            "qtl_beta_alt_orientation_verified":False,
            "valid_for_full_locus_coloc":False,
        })
        qc["POSITION_REF_ALT_EXACT_SNVS"]+=1
    result={"gene":gene,"cell":cell,"source_state":page_state,
            "oasis_plot_snp_total":len(qtl),
            "bokeh_reported_maf_gt_0_5_count":browser_af_greater_half,
            "bokeh_maf_column_means_true_minor_frequency":False if browser_af_greater_half else None,
            "qtl_GWAS_position_ref_alt_snv_matches":len(accepted),
            "source_complete_full_cis_tested_variants_verified":False,
            "qtl_beta_orientation_verified":False,
            "plot_overlap_types":dict(qc)}
    return accepted,result


def run(args):
    best=best_primary_chr12_chain(args.chain)
    lifter=make_mapper(best)
    expected_pos={"rs11066015":111730205,"rs671":111803962,
                  "rs11066132":112030402,"rs77768175":112298314}
    checkpoints=[]
    for rsid,gr37,gr38,pip in CS:
        pos=int(gr37.split(":")[1])
        expected=expected_pos[rsid]
        mapped=lifter(pos)
        checkpoints.append({"rsid":rsid,"grch37":gr37,
                            "mapped_grch38_position":mapped,
                            "previous_source_grch38_position":expected,
                            "independent_chain_coordinate_agrees":mapped==expected})
        if mapped!=expected:
            raise ValueError("UCSC primary chain disagrees with prior source CS coordinates")
    gwas,origin=read_ais_summary(args.gwas)
    mapped,chain_counts=map_ais_variants(gwas,lifter)
    allrows=[];pages=[]
    for gene in GENES:
        for cell in CELLS:
            path=args.oasis_browser/f"OASIS_{gene}_{cell}_original.html"
            source=path.read_bytes()
            rows,summ=compare_one_page(gene,cell,source,mapped)
            summ["original_html_sha256"]=hashlib.sha256(source).hexdigest()
            summ["original_html_path"]=str(path)
            summ["overlap_gwas_chromosome_window"]="CHR12_G0022_ONLY"
            allrows.extend(rows)
            pages.append(summ)
    overview={
        "audit":"IS_G0022_OASIS_5876_AIS_REGIONAL_REAL_VARIANT_OVERLAP_V1",
        "science_gate":"REGIONAL_ALLELE_IDENTITY_QC_NOT_SIGNED_QTL_HARMONIZATION_OR_COLOC",
        "source_ais_study":AIS_STUDY,"trait":AIS_TRAIT,
        "source_GWAS_sha256":hashlib.sha256(args.gwas.read_bytes()).hexdigest(),
        "liftover_chain_sha256":hashlib.sha256(args.chain.read_bytes()).hexdigest(),
        "liftover_chain_primary_id":best["id"],"liftover_chain_primary_score":best["score"],
        "liftover_strategy":best["strategy"],
        "all_four_CS_source_build_checkpoints":checkpoints,
        **origin,"chain_mapping_results":dict(chain_counts),
        "gwas_unique_lifted_positions":len(mapped),
        "browser_page_count":len(pages),"pages":pages,
        "joined_gene_cell_snp_association_rows":len(allrows),
        "distinct_GWAS_SNVS_matched_any_browser_gene_cell":len({
            item["gwas_variant_grch37"] for item in allrows}),
        "unique_genes_in_joint_data":sorted({item["gene"] for item in allrows}),
        "unique_celltypes_in_joint_data":sorted({item["cell"] for item in allrows}),
        "effective_qtl_cell_sample_n":"NOT_VERIFIED",
        "qtl_effect_allele_definition":"NOT_VERIFIED",
        "oasis_bokeh_maf_label_warning":"Browser maf named column has values >0.5: cannot be minor allele frequency; original E-GEAD-1054 cis file names column af. Keep frequency unoriented until official author confirmation.",
        "total_bokeh_maf_gt_half":sum(x["bokeh_reported_maf_gt_0_5_count"] for x in pages),
        "source_complete_all_tested_cis_universe":"NOT_VERIFIED",
        "matched_reference_LD":"NOT_OBTAINED",
        "valid_colocalizations":0,"causal_genes_established":0,
        "interpretation":"All matched SNPs share exact aligned GRCh38 coordinate and REF/ALT strings, NOT proven QTL-beta effect allele or matched LD; cannot claim causal gene."
    }
    args.outdir.mkdir(parents=True,exist_ok=True)
    if allrows:
        with (args.outdir/"G0022_OASIS_EAS_AIS_REGIONAL_ALLELE_MATCHED_SNVS.tsv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter="\t")
            w.writeheader();w.writerows(allrows)
    with (args.outdir/"G0022_OASIS_EAS_AIS_SOURCE_COVERAGE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(pages[0]),delimiter="\t",extrasaction="ignore")
        w.writeheader()
        for r in pages:
            d=dict(r)
            d["plot_overlap_types"]=json.dumps(d["plot_overlap_types"],sort_keys=True)
            w.writerow(d)
    (args.outdir/"G0022_OASIS_EAS_AIS_REGIONAL_OVERLAP_SUMMARY.json").write_text(
        json.dumps(overview,indent=2,ensure_ascii=False)+"\n")
    return overview


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chain",type=Path,required=True)
    ap.add_argument("--gwas",type=Path,required=True)
    ap.add_argument("--oasis-browser",type=Path,required=True)
    ap.add_argument("--outdir",type=Path,required=True)
    result=run(ap.parse_args())
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
