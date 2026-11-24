#!/usr/bin/env python3
"""Synthetic 22-autosome fixtures; never misrepresent unit tests as real GWAS."""
import gzip
import io
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_gwas_1000g_bim_all22 import (
    AUTOSOMES,EXPECTED_GWAS_ROWS,PREVIOUS_CHR12,CATEGORIES,
    stage_bim_by_chr,stage_gwas_by_chr,single_chr_overlap,
    finish_counts,check_previous_chr12,
)


def prepare(root, invalid=None):
    bim=root/"test_ld.tgz"
    records=[]
    with tarfile.open(bim,"w:gz") as tar:
        for chrom in sorted(AUTOSOMES):
            data="".join(f"{chrom} rs{chrom}_{idx} 0 {1000+idx} A G\n" for idx in range(120))
            content=data.encode("utf-8")
            info=tarfile.TarInfo(f"1000G/ref_chr{chrom}.bim")
            info.size=len(content)
            tar.addfile(info,io.BytesIO(content))
    compressed=root/"gwas.tsv.gz"
    rows=0
    with gzip.open(compressed,"wt",encoding="utf-8") as f:
        f.write("chromosome\tbase_pair_location\teffect_allele\tother_allele\n")
        for chrom in sorted(AUTOSOMES):
            for idx in range(120):
                f.write(f"{chrom}\t{1000+idx}\tG\tA\n")
                rows+=1
            f.write(f"{chrom}\t99999\tG\tA\n")
            rows+=1
        if invalid is not None:
            f.write(invalid)
            rows+=1
    return bim,compressed,rows


class All22SourceTests(unittest.TestCase):
    def test_explicit_22_autosomes(self):
        self.assertEqual(len(AUTOSOMES),22)
        self.assertNotIn(23,AUTOSOMES)

    def test_full_synthetic_22_sequential_mem_bounded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            ref,gwas,rows=prepare(root)
            bim,bt,paths=stage_bim_by_chr(ref,root)
            self.assertEqual(len(bim),22)
            self.assertEqual(len(paths),22)
            self.assertEqual(bt["all_autosomal_bim_rows"],22*120)
            with patch.dict(EXPECTED_GWAS_ROWS,{"EUR":rows}):
                snps,gt=stage_gwas_by_chr(gwas,root,"EUR")
            self.assertEqual(gt["all_original_gwas_rows"],22*121)
            allrows={}
            for c in AUTOSOMES:
                allrows[str(c)]=single_chr_overlap(c,bim[c],snps[c])
                self.assertEqual(allrows[str(c)]["UNORDERED_ALLELE_SET_EXACT"],120)
                self.assertEqual(allrows[str(c)]["COORDINATE_NOT_IN_LD_BIM"],1)
            total=finish_counts(allrows,gt)
            self.assertEqual(total["all_22_autosome_gwas_rows"],22*121)
            self.assertEqual(total["all_22_autosome_BIM_rows"],22*120)
            self.assertEqual(total["UNORDERED_ALLELE_SET_EXACT"],22*120)
            self.assertEqual(total["22_autosome_coordinate_overlap_pct"],round(100*120/121,4))

    def test_gwas_publisher_rows_cannot_silently_drift(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            _,gwas,_=prepare(root)
            with self.assertRaisesRegex(ValueError,"row count drift"):
                stage_gwas_by_chr(gwas,root,"EUR")

    def test_previous_chr12_real_source_expected_counts_guard(self):
        fake={str(c):{"gwas_rows":100} for c in (1,2)}
        with self.assertRaisesRegex(ValueError,"previous verified CI mismatch"):
            check_previous_chr12("EUR",fake)

    def test_reconcile_previous_real_chr12_counts(self):
        e=PREVIOUS_CHR12["EUR"]
        one={"gwas_rows":e["selected_chr_GWAS_variant_rows"]}
        for k in CATEGORIES:
            one[k]=e.get(k,0)
        two={"gwas_rows":0}
        for k in CATEGORIES:
            two[k]=0
        self.assertEqual(check_previous_chr12("EUR",{"1":one,"2":two})["status"],"PASS")

    def test_non_numeric_chr_coordinate_is_fatal(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            _,gwas,rows=prepare(root,invalid="1\tNaN\tA\tG\n")
            with patch.dict(EXPECTED_GWAS_ROWS,{"EUR":rows}):
                with self.assertRaisesRegex(ValueError,"noninteger"):
                    stage_gwas_by_chr(gwas,root,"EUR")

    def test_no_false_success_for_gwas_missing_autosome(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            path=root/"small.tsv.gz"
            with gzip.open(path,"wt",encoding="utf-8") as f:
                f.write("chromosome\tbase_pair_location\teffect_allele\tother_allele\n")
                for i in range(120):
                    f.write(f"1\t{1000+i}\tA\tG\n")
            with patch.dict(EXPECTED_GWAS_ROWS,{"EUR":120}):
                with self.assertRaisesRegex(ValueError,"few GWAS rows"):
                    stage_gwas_by_chr(path,root,"EUR")

if __name__=="__main__":
    unittest.main()
