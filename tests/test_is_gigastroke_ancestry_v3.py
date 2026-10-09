"""Regression: GIGASTROKE per-file YAML, never inferred ancestry by accession."""
import csv,json,tempfile,unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from verify_gigastroke_ancestry_v3 import build

def add(path,kind,n):
    path.parent.mkdir(parents=True,exist_ok=True)
    parts=[f'  - sample_ancestry_category:\n      - {k}\n    sample_size: {v}\n' for k,v in zip(kind,n)]
    path.write_text('genome_assembly: GRCh37\nsamples:\n'+''.join(parts)+'data_file_md5sum: 0123456789abcdef0123456789abcdef\n')

class StudyYamlTest(unittest.TestCase):
    def test_cross_ancestry_vs_european_and_preserve_mismatch(self):
        with tempfile.TemporaryDirectory() as t:
            base=Path(t);meta=base/"meta";meta.mkdir()
            rows=[dict(accession="GCST00001",phenotype_class="AIS",sample_size="100",ancestry_class="AFR_AMR"),
                  dict(accession="GCST00002",phenotype_class="CES",sample_size="200",ancestry_class="UNKNOWN")]
            with (meta/"GIGASTROKE_STUDY_MAP_V2.tsv").open("w") as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
            add(meta/"yaml/GCST00001_buildGRCh37.tsv.gz-meta.yaml",
                ["East Asian","European"],[100,200])
            add(meta/"yaml/GCST00002_buildGRCh37.tsv.gz-meta.yaml",
                ["European"],[200])
            s=build(meta,base/"data",base/"out")
            self.assertEqual(s["studies"],2)
            self.assertEqual(s["n_sample_size_mismatches"],1)
            self.assertEqual(s["ancestry_breakdown"],{"EUR":1,"MULTI_ANCESTRY":1})
            with (base/"out/GIGASTROKE_SOURCE_VERIFIED_V3.tsv").open() as f:
                records=list(csv.DictReader(f,delimiter="\t"))
            self.assertEqual(records[0]["subpopulation_n"],"EAS:100;EUR:200")
            self.assertEqual(records[1]["ancestry_verified"],"EUR")
            self.assertEqual(records[0]["independent_replication"],"NO_OVERLAP_UNKNOWN")
if __name__=="__main__":
    unittest.main()
