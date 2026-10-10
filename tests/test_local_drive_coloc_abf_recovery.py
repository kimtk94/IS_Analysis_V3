"""Offline backup manifest, allele/SNP schema and fail-closed recovery tests."""
from __future__ import annotations
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

SOURCE=Path(__file__).resolve().parents[1]/"scripts/is/audit_local_drive_coloc_abf_recovery.py"
spec=importlib.util.spec_from_file_location("coloc_recovery",SOURCE)
mod=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=mod
spec.loader.exec_module(mod)

HEADER=["locus","dataset_key","gene_base","eqtl_beta","eqtl_se",
        "gwas_beta","gwas_se","match_key","qtl_n_scalar","harmonization"]
def record(locus="BBJ_IS_L001",key="GTEx_V8__Artery_Aorta",gid="ENSG00000138675"):
    return dict(zip(HEADER,[locus,key,gid,"0.05","0.02","0.08","0.03",
                            "4:100:T:C","175","EXACT_REF_ALT_GRCH38"]))

class TestColocSourceRecovery(unittest.TestCase):
    def test_qc_pass_exact_variant_and_stats(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"input.tsv"
            r=record()
            p.write_text("\t".join(HEADER)+"\n"+"\t".join(r[k] for k in HEADER)+"\n")
            ok,n,status=mod.inspect(p,{"locus":r["locus"],
                "dataset_key":r["dataset_key"],"gene_base":r["gene_base"],
                "expected_nsnps":"1"})
            self.assertTrue(ok)
            self.assertEqual(n,1)
            self.assertEqual(status,"SCHEMA_SNP_COUNT_STATS_QC_PASS")

    def test_fail_on_duplicate_snps(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"input.tsv"
            r=record()
            line="\t".join(r[k] for k in HEADER)
            p.write_text("\t".join(HEADER)+"\n"+line+"\n"+line+"\n")
            ok,n,status=mod.inspect(p,{"locus":r["locus"],
                "dataset_key":r["dataset_key"],"gene_base":r["gene_base"],
                "expected_nsnps":"2"})
            self.assertFalse(ok)
            self.assertEqual(status,"DUPLICATE_OR_MISSING_VARIANT")

    def test_646_filename_universe_and_missing_copies(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            local=root/"local";local.mkdir()
            recovery=root/"recovery";recovery.mkdir()
            master=[];replay=[];drive=set()
            for i in range(646):
                locus="BBJ_IS_L001"
                ds="GTEx_V8__Artery_Aorta"
                gene=f"ENSG{i:011d}"
                file=f"{locus}__{ds}__{gene}.tsv"
                master.append({"locus":locus,"dataset_key":ds,"gene_base":gene})
                replay.append({"locus":locus,"dataset_key":ds,"gene_base":gene,
                               "filename":file,"expected_nsnps":"1","status":"MISSING_INPUT"})
                drive.add(file)
            # One existing locally verified input.
            first=replay[0]
            r=record(gid=first["gene_base"])
            (local/first["filename"]).write_text(
                "\t".join(HEADER)+"\n"+"\t".join(r[k] for k in HEADER)+"\n")
            result=mod.audit(master,replay,sorted(drive),local,recovery)
            self.assertEqual(len(result),646)
            self.assertEqual(sum(r["source_QC_status"]=="LOCAL_SOURCE_QC_PASS" for r in result),1)
            self.assertEqual(sum(r["source_QC_status"]=="DRIVE_PRESENT_NOT_YET_COPIED" for r in result),645)
            with self.assertRaises(ValueError):
                mod.audit(master,replay,sorted(drive)[:-1],local,recovery)

if __name__=="__main__":
    unittest.main()
