"""Fixture-only fail-closed tests for ImmuNexUT/JCTF competing gene audit."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/is"))
from audit_g0022_immunexut_competing_genes import (
    VARIANTS, parse_variant_page, load_jctf, build_summary, acquire_source
)

def csv_file(path, records):
    with path.open("w") as f:
        w=csv.DictWriter(f, fieldnames=list(records[0]), delimiter="\t")
        w.writeheader();w.writerows(records)

class TestCompetingGenes(unittest.TestCase):
    def html(self, rsid="rs671", gene="BRAP", p=1e-20):
        return ("<snps-eqtl-table-component eqtl-data ='"+
                json.dumps([dict(snp_id=rsid, gene_symbol=gene,cell_type="Neu",
                     chromosome="chr12",position=VARIANTS[rsid][1],
                     eqtl_pval1=p,eqtl_effect_beta1=-0.3)],separators=(",",":"))+
                "'></snps-eqtl-table-component>")

    def test_variant_identity_and_gene(self):
        rows=parse_variant_page(self.html(),"rs671")
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]["gene"],"BRAP")
        self.assertFalse(rows[0]["source_effect_allele_known"])
        self.assertEqual(rows[0]["qtl_status"],"SUPPORTED_MARGINAL_QTL_NOT_COLOC")

    def test_wrong_snp_position_or_duplicate_blocked(self):
        with self.assertRaises(ValueError):
            parse_variant_page(self.html("rs671"),"rs11066015")
        broken=self.html().replace("111803962","111803961")
        with self.assertRaises(ValueError):
            parse_variant_page(broken,"rs671")
        repeated=self.html().replace("}]", "},"+json.dumps(dict(snp_id="rs671",gene_symbol="BRAP",cell_type="Neu",
            chromosome="chr12",position=111803962,eqtl_pval1=1e-12,
            eqtl_effect_beta1=-0.1),separators=(",",":"))[1:])
        with self.assertRaises((ValueError, TypeError)):
            parse_variant_page(repeated,"rs671")

    def test_jctf_source_filters_and_ensemble(self):
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/"source.tsv"
            cs=VARIANTS["rs671"][0]
            csv_file(f,[
                dict(gene_symbol="BRAP",GWAS_variant_GRCh37=cs,qtl_category="pQTL",
                     jctf_allele_association_p="2.61e-5",jctf_susie_variant_pip=".0767"),
                dict(gene_symbol="RPH3A",GWAS_variant_GRCh37=cs,qtl_category="eQTL",
                     jctf_allele_association_p="1.6e-21",jctf_susie_variant_pip="0")])
            j=load_jctf(f)
            rows, status=build_summary(parse_variant_page(self.html(),"rs671"),j)
            by={r["gene"]:r for r in rows}
            self.assertEqual(by["BRAP"]["jctf_pqtl_marker_count"],1)
            self.assertEqual(by["RPH3A"]["jctf_eqtl_marker_count"],1)
            self.assertFalse(by["BRAP"]["causal_gene_established"])
            self.assertEqual(status["valid_colocs"],0)
            self.assertEqual(by["ALDH2"]["immunexut_variants_associated"],0)

    def test_existing_source_requires_explicit_overwrite(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/"rs671.html"
            p.write_text(self.html())
            with self.assertRaises(FileExistsError):
                acquire_source(p,"rs671",fetch=True,overwrite=False)
            self.assertIn(b"BRAP",acquire_source(p,"rs671",fetch=False,overwrite=False))

if __name__=="__main__":
    unittest.main()
