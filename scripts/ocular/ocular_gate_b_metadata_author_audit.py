#!/usr/bin/env python3
"""Gate B metadata-only cell-selection and author-DEG inventory audit.

Reads public GEO cell metadata + SOFT. Reports aggregate threshold changes,
not donor identities or individual phenotype/genotype. Never infers an
original 80-gene list from the GEO abstract alone.
"""
import argparse,csv,gzip,hashlib,json,re,tempfile
from collections import Counter,defaultdict
from pathlib import Path
from ocular_p4_fullcohort import collect_metadata
from ocular_p1_integrity import BASE,download
from ocular_p0_audit import fetch,root

N=169398
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ocular-gate-b-") as td:
        fp=Path(td)/"cell_metadata.tsv.gz"
        d=download(BASE+"series/GSE313nnn/GSE313736/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",
                   fp,12*1024*1024)
        selections={}
        pools={}
        pool_dist={}
        for score in (0.50,0.75):
            meta,lib_count,donor_pools=collect_metadata(fp,score)
            by=defaultdict(Counter)
            for (lib,arm),barcodes in meta.items():
                for donor in barcodes.values():by[donor][arm]+=1
            for donor in by:
                if donor not in donor_pools or len(donor_pools[donor])!=1:
                    raise ValueError("MULTIPOOL_DONOR")
            n30={donor for donor,c in by.items() if min(c["dex"],c["no-dex"])>=30}
            n50={donor for donor,c in by.items() if min(c["dex"],c["no-dex"])>=50}
            n10={donor for donor,c in by.items() if min(c["dex"],c["no-dex"])>=10}
            selections[score]={"eligible_10":n10,"eligible_30":n30,"eligible_50":n50}
            pool_dist[score]=dict(Counter(next(iter(donor_pools[d])) for d in n30))
            pools[score]={"score":score,"metadata_paired_donors_10cells":len(n10),
              "metadata_paired_donors_30cells":len(n30),
              "metadata_paired_donors_50cells":len(n50),
              "selected_barcodes_30_eligible":{
                 arm:sum(c[arm] for d,c in by.items() if d in n30)
                 for arm in ("no-dex","dex")},
              "eligibility_source":"metadata only; actual 10x barcode intersection separately validated at cutoff 0.75",
              "source_library_count":len(lib_count)}
        low=selections[.50]["eligible_30"];high=selections[.75]["eligible_30"]
        if len(high)!=57 or not high.issubset(low):
            raise ValueError("BASELINE_57_OR_NESTED_CONFIDENCE_ELIGIBILITY_DRIFT")
        raw,tx=fetch(root("GSE313736")+"soft/GSE313736_family.soft.gz",
                     limit=10*1024*1024)
        plain=gzip.decompress(raw).decode("utf8","replace")
        series_lines=[line for line in plain.splitlines() if
                      line.startswith(("!Series_supplementary_file","!Series_relation","!Series_pubmed_id",
                                      "!Series_summary","!Series_overall_design"))]
        files=[v.strip() for line in series_lines if line.startswith("!Series_supplementary_file")
               for v in [line.split("=",1)[1]]]
        deg_like=[f for f in files if re.search(r"(?:differential|DEGs?|edgeR|fold[_-]?change|logFC|DE[_-]?results)",f,re.I)]
        # A series abstract's 63+17 is a count, not a gene list. No hand-coded reconstruction.
        original_reported={"up":63,"down":17,"total":80,
             "ambiguous_reported_downregulation_logFC_threshold":True,
             "exact_list_accessible_in_series_supplement_inventory":bool(deg_like),
             "potential_DEG_list_supplement_urls":deg_like,
             "original_80_gene_identity_match":"HOLD_SOURCE_AUTHOR_GENE_LIST_UNAVAILABLE"}
        result={"status":"PASS_GATE_B_PUBLIC_METADATA_QC_AUTHOR_LIST_HOLD",
                "data":"GSE313736","total_GEO_cell_metadata_rows":N,
                "source_SOFT_sha256":hashlib.sha256(raw).hexdigest(),
                "public_series_supplement_names":[f.rsplit("/",1)[-1] for f in files],
                "public_series_supplement_count":len(files),
                "score0_50":pools[.50],"score0_75":pools[.75],
                "score_0_75_donors_subset_of_0_50":high.issubset(low),
                "donors_added_at_0_50_vs_0_75":len(low-high),
                "donors_lost_at_0_50_vs_0_75":len(high-low),
                "source_pools_with_selected_donors_score0_50":pool_dist[.50],
                "source_pools_with_selected_donors_score0_75":pool_dist[.75],
                "original_publication_abstract_claim":original_reported,
                "official_GEO_sample_matrix_18_complete_in_prior_P4":True,
                "no_GEX_expression_difference_estimated_in_this_job":True,
                "no_genotypes_or_donor_identifiers_uploaded":True,
                "next_gate":"Raw 10x pseudo-count reaggregation at score .50 vs fixed .75, then pair-adjusted edgeR; never use metadata-only cell fractions as gene-effect replication"}
        (args.out/"GATE_B_BEAMA_050_075_ELIGIBILITY_AND_AUTHOR_SOURCE.json").write_text(
            json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
        print("GATE_B",result["status"],"0.50 n30",len(low),"0.75 n30",len(high),
              "original_de_gene_list_in_GEO_supplement",bool(deg_like),flush=True)
if __name__=="__main__":main()
