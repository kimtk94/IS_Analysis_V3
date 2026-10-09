"""GENCODE window mapping fixture, protein and noncoding genes retained."""
import csv,gzip,json,subprocess,sys,tempfile,unittest
from pathlib import Path
SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/map_all_genes_gencode_v19.py"
class TestGeneMap(unittest.TestCase):
 def test_coding_noncoding_overlap_and_boundary(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);f=root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv"
   cols=["group_id","chr","window_start","window_end","lead_variant","region_start","region_end"]
   with f.open("w") as h:
    w=csv.DictWriter(h,fieldnames=cols,delimiter="\t");w.writeheader()
    w.writerow(dict(group_id="A",chr="1",window_start=100,window_end=300,lead_variant="1:150:A:G",region_start=140,region_end=160))
   gtf=root/"tiny.gtf.gz"
   with gzip.open(gtf,"wt") as h:
    h.write('1\tsrc\tgene\t101\t180\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE_A"; gene_type "protein_coding";\n')
    h.write('1\tsrc\tgene\t300\t400\t.\t-\t.\tgene_id "GENE2"; gene_name "LNC_A"; gene_type "lincRNA";\n')
    h.write('1\tsrc\tgene\t301\t500\t.\t+\t.\tgene_id "GENE3"; gene_name "OUT"; gene_type "protein_coding";\n')
   r=subprocess.run([sys.executable,str(SCRIPT),"--root",str(root),"--gtf",str(gtf)],capture_output=True,text=True)
   self.assertEqual(r.returncode,0,r.stderr)
   report=json.loads((root/"IS_ALL_GENE_UNIVERSE_SUMMARY.json").read_text())
   self.assertEqual(report["region_gene_associations"],2)
   self.assertEqual(report["coding_gene_ids"],1)
   self.assertEqual(report["noncoding_or_other_gene_ids"],1)
if __name__=="__main__":unittest.main()
