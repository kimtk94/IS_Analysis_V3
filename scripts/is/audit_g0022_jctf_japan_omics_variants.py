#!/usr/bin/env python3
"""JCTF Japan Omics Browser actual rs671/4-GWAS-CS variant QTL records.

Retains raw public HTML SHA-256, exact query identity, ALDH2 eQTL / pQTL
effects and SuSiE/FINEMAP variant PIP; not full locus molecular colocalization.
Browser QTL may be affected by protein assay epitope binding and by COVID
case sampling. Never assert causal ALDH2--stroke mediation from marginal QTL.
"""
import argparse,csv,hashlib,json,math,re
from html.parser import HTMLParser
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
SRC=Path("/srv/is-analysis/data/is/qtl/japan_omics/G0022_v1")
OUT=ROOT/"jctf_japan_omics_variants_v1"
PAGES={
 "12:112168009:G:A":("rs11066015","rs11066015.html"),
 "12:112241766:G:A":("rs671","rs671.html"),
 "12:112468206:C:T":("chr12:112030402:C:T","chr12_112030402_C_T.html"),
 "12:112736118:A:G":("chr12:112298314:A:G","chr12_112298314_A_G.html")
}
HEADERS=["Gene ID","Gene Symbol","TSS Distance","Category","P-value",
    "Effect Size","Effect Size (se)","PIP (SuSiE)","PIP (FINEMAP)","MAF"]
class HtmlTable(HTMLParser):
    def __init__(self,table_id="target"):
        super().__init__(convert_charrefs=True)
        self.table_id=table_id
        self.active=False;self.cell=None;self.row=None;self.rows=[]
    def handle_starttag(self,tag,attrs):
        if tag=="table" and dict(attrs).get("id")==self.table_id:
            if self.active:raise ValueError("Nested target table")
            self.active=True
        elif self.active and tag=="tr":
            self.row=[]
        elif self.active and tag in ("td","th") and self.row is not None:
            self.cell=[]
    def handle_data(self,data):
        if self.cell is not None:self.cell.append(data)
    def handle_endtag(self,tag):
        if self.active and tag in ("td","th") and self.cell is not None:
            self.row.append(" ".join(" ".join(self.cell).split()))
            self.cell=None
        elif self.active and tag=="tr" and self.row is not None:
            if self.row:self.rows.append(self.row)
            self.row=None
        elif self.active and tag=="table":
            self.active=False
def parse_html(text,query):
    if f"Genes within 1Mb ({query})" not in text:
        raise ValueError("Japan Omics variant lookup did not return requested identifier")
    parser=HtmlTable();parser.feed(text)
    if len(parser.rows)<2 or parser.rows[0]!=HEADERS:
        raise ValueError("Japan Omics Browser target table schema changed")
    out=[]
    for row in parser.rows[1:]:
        if len(row)!=len(HEADERS):
            raise ValueError("Malformed source table row")
        item=dict(zip(HEADERS,row))
        if not re.fullmatch(r"ENSG\d{11}",item["Gene ID"]):
            raise ValueError("Invalid QTL gene ID")
        if item["Category"] not in ("eQTL","pQTL"):
            raise ValueError("Unexpected QTL assay type")
        for fld in ("P-value","Effect Size","Effect Size (se)","PIP (SuSiE)",
                    "PIP (FINEMAP)","MAF"):
            try:v=float(item[fld])
            except (ValueError,TypeError):
                raise ValueError("Invalid source numeric field "+fld)
            if not math.isfinite(v):raise ValueError("Nonfinite molecular QTL result")
        if not (0<=float(item["P-value"])<=1 and
                0<=float(item["PIP (SuSiE)"])<=1 and
                0<=float(item["PIP (FINEMAP)"])<=1 and
                0<=float(item["MAF"])<=.5 and float(item["Effect Size (se)"])>0):
            raise ValueError("Invalid QTL posterior/frequency/effect SE")
        out.append(item)
    if len({(r["Gene ID"],r["Category"]) for r in out})!=len(out):
        raise ValueError("Gene-category duplicate ambiguous JCTF result")
    return out
def build(source,out,root):
    if not source.is_dir():raise FileNotFoundError(source)
    gw=json.loads((root/"G0022_FULL_AIS_SUSIE_SUMMARY.json").read_text())
    if gw.get("status")!="FULL_LOCUS_EXPLORATORY_ONLY":
        raise ValueError("GWAS full-locus evidence status changed")
    cs={v for arr in gw["credible_set_variant_ids"].values() for v in arr}
    if cs!=set(PAGES):raise ValueError("JCTF four source variants disagree with GWAS CS")
    with (root/"G0022_FULL_AIS_PIP.tsv").open() as f:
        model_pip={r["variant_id"]:float(r["pip"]) for r in csv.DictReader(
          f,delimiter="\t") if r["variant_id"] in cs}
    data=[];provenance=[]
    for variant,(query,name) in PAGES.items():
        path=source/name
        if not path.is_file():raise FileNotFoundError(path)
        content=path.read_text()
        rows=parse_html(content,query)
        provenance.append({
            "variant_grch37":variant,"jctf_query":query,
            "query_url":"https://japan-omics.jp/variant?input_value="+query,
            "raw_filename":name,"source_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "source_rows":len(rows),
            "source_status":"PUBLIC_VARIANT_TABLE_READ_ONLY"})
        for row in rows:
            data.append({
              "GWAS_variant_GRCh37":variant,"JOB_query_variant":query,
              "gwas_ais_full_locus_pip":model_pip[variant],
              "gene_id":row["Gene ID"],"gene_symbol":row["Gene Symbol"],
              "qtl_category":row["Category"],"gene_tss_distance_source_bp":row["TSS Distance"],
              "jctf_allele_association_p":float(row["P-value"]),
              "jctf_beta_effect_allele_not_revalidated":float(row["Effect Size"]),
              "jctf_beta_se":float(row["Effect Size (se)"]),
              "jctf_susie_variant_pip":float(row["PIP (SuSiE)"]),
              "jctf_finemap_variant_pip":float(row["PIP (FINEMAP)"]),
              "jctf_reported_maf":float(row["MAF"]),
              "qtl_study_source":"JCTF_Japan_Omics_Browser",
              "qtl_effect_allele_orientation_verified":0,
              "causal_stroke_coloc_status":"NOT_TESTED_FULL_QTL_SUMSTATS_REQUIRED"})
    if not data:raise ValueError("No JCTF browser source rows parsed")
    out.mkdir(parents=True,exist_ok=True)
    for p,items in [
      (out/"G0022_JCTF_4SNP_JAPANESE_EQTL_PQTL_EVIDENCE.tsv",data),
      (out/"G0022_JCTF_4SNP_RAW_SOURCE_MANIFEST.tsv",provenance)]:
        with p.open("w",newline="") as f:
            w=csv.DictWriter(f,delimiter="\t",fieldnames=list(items[0]))
            w.writeheader();w.writerows(items)
    aldh2=[r for r in data if r["gene_symbol"]=="ALDH2"]
    paired={(v,c):r for r in aldh2 for v,c in
            [(r["GWAS_variant_GRCh37"],r["qtl_category"])]}
    if not all((var,cat) in paired for var in cs for cat in ("eQTL","pQTL")):
        raise RuntimeError("Incomplete Japanese ALDH2 QTL evidence for all CS")
    result={
      "source":"Japan_Omics_Browser_JCTF_Wang_NatureGenetics_2024",
      "is_East_Asian_QTL":True,
      "clinical_sampling":"Japanese_COVID19_TaskForce_cases",
      "source_variant_pages":4,"source_table_records":len(data),
      "ALDH2_eQTL_variant_records":sum(x["qtl_category"]=="eQTL" for x in aldh2),
      "ALDH2_pQTL_variant_records":sum(x["qtl_category"]=="pQTL" for x in aldh2),
      "rs671_ALDH2_eQTL_p":paired[("12:112241766:G:A","eQTL")]["jctf_allele_association_p"],
      "rs671_ALDH2_eQTL_susie_pip":paired[("12:112241766:G:A","eQTL")]["jctf_susie_variant_pip"],
      "rs671_ALDH2_pQTL_p":paired[("12:112241766:G:A","pQTL")]["jctf_allele_association_p"],
      "rs671_ALDH2_pQTL_susie_pip":paired[("12:112241766:G:A","pQTL")]["jctf_susie_variant_pip"],
      "colocalization_numerically_calculated":False,
      "validated_causal_stroke_genes":0,
      "source_data_limitation":"Variant-wise public browser values; full summary/LD and assay effect allele contract required for multi-signal coloc; COVID sampling and pQTL epitope effects.",
      "status":"EAS_JCTF_VARIANT_QTL_EVIDENCE_EXTRACTED_NO_COLOC_CLAIM"}
    (out/"G0022_JCTF_4SNP_QTL_EVIDENCE_SUMMARY.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    for x in aldh2:print(x["GWAS_variant_GRCh37"],x["qtl_category"],x["jctf_allele_association_p"],x["jctf_susie_variant_pip"])
    return result
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,default=SRC)
    p.add_argument("--out",type=Path,default=OUT)
    p.add_argument("--root",type=Path,default=ROOT)
    a=p.parse_args()
    build(a.source,a.out,a.root)
