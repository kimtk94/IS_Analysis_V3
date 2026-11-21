#!/usr/bin/env python3
"""Evidence gate for binary effective N for TPMI Phecode433.21 vs GIGASTROKE EAS AIS.

Do not substitute phenotype portal counts for PLINK2 variant-level OBS_CT.
Do not estimate case/control among OBS_CT. Shared samples and case-definition
equivalence must be demonstrated before LDSC, MR or colocalisation.
"""
from __future__ import annotations
import argparse,html,json,re,urllib.request
from html.parser import HTMLParser
from pathlib import Path

TPMI_URL="https://pheweb.ibms.sinica.edu.tw/pheno/433.21"
AIS_URL="https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90104545"
ORIGINAL_TPMI_OBS_CT=204099
PAGE_EXPECTED_CASES=9249
PAGE_EXPECTED_CONTROLS=304660
AIS_DISCOVERY_CASES=19032
AIS_DISCOVERY_CONTROLS=237242
TPMI_PHENO="Cerebral artery occlusion, with cerebral infarction"

class Visible(HTMLParser):
    def __init__(self):
        super().__init__()
        self.text=[]
        self.skip=[]
    def handle_starttag(self,tag,attrs):
        if tag in ("script","style"):self.skip.append(tag)
    def handle_endtag(self,tag):
        if self.skip and self.skip[-1]==tag:self.skip.pop()
    def handle_data(self,data):
        if not self.skip:self.text.append(data)

def get(url,accept):
    req=urllib.request.Request(
        url,headers={"User-Agent":"TPMI-43321-EAS-AIS-Neff-evidence-gate/1.0",
                     "Accept":accept})
    with urllib.request.urlopen(req,timeout=45) as resp:
        if resp.status!=200 or resp.geturl()!=url:
            raise ValueError("source HTTP/redirect changed")
        b=resp.read(300001)
        if len(b)>300000: raise ValueError("source response too large")
        return b.decode("utf-8-sig")

def audit():
    page=get(TPMI_URL,"text/html")
    parser=Visible();parser.feed(page)
    text=" ".join(" ".join(parser.text).split())
    match=re.search(r"(\d[\d,]*)\s+cases,\s+(\d[\d,]*)\s+controls",text,re.I)
    if not match or TPMI_PHENO not in text:
        raise ValueError("TPMI source phenotype or sample counts missing")
    cases=int(match.group(1).replace(",",""))
    controls=int(match.group(2).replace(",",""))
    if (cases,controls)!=(PAGE_EXPECTED_CASES,PAGE_EXPECTED_CONTROLS):
        raise ValueError("TPMI official phenotype case/control counts changed")
    dat=json.loads(get(AIS_URL,"application/json"))
    if dat.get("accession_id")!="GCST90104545":
        raise ValueError("GIGASTROKE AIS catalogue identity mismatch")
    traits=dat.get("efo_traits",[])
    if "HP_0002140" not in {t.get("efo_id") for t in traits if isinstance(t,dict)}:
        raise ValueError("GIGASTROKE AIS ontology mismatch")
    initial=dat.get("initial_sample_size","")
    expected=(f"{AIS_DISCOVERY_CASES:,} East Asian ancestry cases, "
              f"{AIS_DISCOVERY_CONTROLS:,} East Asian ancestry controls")
    if initial!=expected:raise ValueError("GIGASTROKE EAS reported counts changed")
    page_n=cases+controls
    obs=ORIGINAL_TPMI_OBS_CT
    return {
      "schema":"TPMI-Phecode43321-EAS-AIS-effective-N-source-audit-v1",
      "source_status":"PASS_OFFICIAL_PHENOTYPE_AND_CATALOG_METADATA",
      "TPMI_phenotype_code":"433.21",
      "TPMI_phenotype_label":TPMI_PHENO,
      "TPMI_PheWeb_page_cases":cases,
      "TPMI_PheWeb_page_controls":controls,
      "TPMI_PheWeb_page_total":page_n,
      "TPMI_original_PLINK2_additive_OBS_CT":obs,
      "TPMI_phenotype_total_minus_tested_OBS_CT":page_n-obs,
      "TPMI_case_control_breakdown_within_OBS_CT":None,
      "TPMI_source_per_SNP_Neff_verified":False,
      "EAS_AIS_discovery_cases_reported":AIS_DISCOVERY_CASES,
      "EAS_AIS_discovery_controls_reported":AIS_DISCOVERY_CONTROLS,
      "EAS_AIS_per_variant_case_control_counts":None,
      "EAS_AIS_effective_N_verified":False,
      "TPMI_and_GIGASTROKE_cohort_independence_verified":False,
      "phenotype_operational_case_definition_equivalent":False,
      "LDSC_h2_neff_gate":"HOLD_UNVERIFIED_BINARY_CASE_CONTROL_IN_REGRESSION",
      "MR_coloc_neff_gate":"HOLD_REQUIRED_SNP_LEVEL_INPUTS",
      "LDSC_MR_coloc_executed":False,
      "effective_N_imputed":False,
      "scientific_status":"HOLD",
      "evidence_sources":[TPMI_URL,AIS_URL],
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    try:result=audit()
    except Exception as e:
        result={"schema":"TPMI-Phecode43321-EAS-AIS-effective-N-source-audit-v1",
                "source_status":"HOLD_SOURCE_METADATA",
                "error_type":type(e).__name__,"error":str(e)[:280],
                "effective_N_imputed":False,"scientific_status":"HOLD"}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"status":result["source_status"],
                      "TPMI_PheWeb_N":result.get("TPMI_PheWeb_page_total"),
                      "TPMI_PLINK2_OBS_CT":result.get("TPMI_original_PLINK2_additive_OBS_CT"),
                      "N_difference":result.get("TPMI_phenotype_total_minus_tested_OBS_CT"),
                      "approved_Neff":result.get("TPMI_source_per_SNP_Neff_verified",False),
                      "error":result.get("error"),
                      "scientific_status":"HOLD"},ensure_ascii=False))
if __name__=="__main__":main()
