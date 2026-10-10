"""JCTF source-stamped rs671 and four GWAS-CS molecular-QTL parser tests."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_jctf_japan_omics_variants import parse_html,build,PAGES,HEADERS

def source_html(query,rows):
    def htmlrow(fields,tag):
        return "<tr>"+"".join("<"+tag+">"+str(x)+"</"+tag+">" for x in fields)+"</tr>"
    return ('<html><body><h2>Genes within 1Mb ('+query+')</h2>'
       '<table id="target"><thead>'+htmlrow(HEADERS,"th")+'</thead><tbody>'
       +"".join(htmlrow(row,"td") for row in rows)
       +'</tbody></table></body></html>')
class JCTF(unittest.TestCase):
    def fixture_row(self,category):
        return ["ENSG00000111275","ALDH2","30000",category,"1e-9",
                "-0.2","0.03","0.35","0.3","0.24"]
    def test_parser_accepts_eQTL_pQTL_and_rejects_unexpected_lookup(self):
        example=source_html("rs671",[self.fixture_row("eQTL"),
                                      self.fixture_row("pQTL")])
        a=parse_html(example,"rs671")
        self.assertEqual(len(a),2)
        self.assertEqual(a[1]["Category"],"pQTL")
        with self.assertRaisesRegex(ValueError,"did not return requested"):
            parse_html(example,"rs11066015")
    def test_rejects_bad_pip_and_gene_duplicate(self):
        a=self.fixture_row("eQTL")
        invalid=a[:];invalid[7]="1.5"
        with self.assertRaisesRegex(ValueError,"posterior"):
            parse_html(source_html("rs671",[invalid]),"rs671")
        with self.assertRaisesRegex(ValueError,"duplicate"):
            parse_html(source_html("rs671",[a,a]),"rs671")
    def test_source_file_provenance_and_no_false_coloc(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td);root=base/"full";root.mkdir()
            src=base/"source";src.mkdir()
            target=base/"out"
            variants=list(PAGES)
            (root/"G0022_FULL_AIS_SUSIE_SUMMARY.json").write_text(json.dumps({
                "status":"FULL_LOCUS_EXPLORATORY_ONLY",
                "credible_set_variant_ids":{"L1":variants}}))
            with (root/"G0022_FULL_AIS_PIP.tsv").open("w") as f:
                w=csv.DictWriter(f,fieldnames=["variant_id","pip"],delimiter="\t")
                w.writeheader()
                w.writerows({"variant_id":v,"pip":.25} for v in variants)
            for var,(query,filename) in PAGES.items():
                (src/filename).write_text(source_html(query,[self.fixture_row("eQTL"),
                                                           self.fixture_row("pQTL")]))
            res=build(src,target,root)
            self.assertEqual(res["ALDH2_eQTL_variant_records"],4)
            self.assertEqual(res["ALDH2_pQTL_variant_records"],4)
            self.assertFalse(res["colocalization_numerically_calculated"])
            self.assertEqual(res["validated_causal_stroke_genes"],0)
if __name__=="__main__":unittest.main()
