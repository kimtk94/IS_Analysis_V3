#!/usr/bin/env python3
"""No live sources. Synthetic BIM/GWAS fixtures, including allele safety gates."""
import gzip
import io
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_gwas_1000g_bim_overlap import (
    canonical_chrom,pair,complement_pair,classify_variant,
    reference_bim_coordinates,classify_gwas,SOURCES,
)


class GwasBimOverlapTests(unittest.TestCase):
    def test_expected_original_source_ids(self):
        self.assertEqual(SOURCES["EUR"]["gwas_accession"],"GCST90104540")
        self.assertEqual(SOURCES["EAS"]["gwas_accession"],"GCST90104545")
        self.assertEqual(SOURCES["EAS"]["ld_archive"],"1000G_Phase3_EAS_plinkfiles.tgz")

    def test_chromosome_and_variant_normalization(self):
        self.assertEqual(canonical_chrom("chr1"),1)
        self.assertIsNone(canonical_chrom("X"))
        self.assertIsNone(canonical_chrom("25"))
        self.assertEqual(pair("A","G"),frozenset(("A","G")))
        self.assertIsNone(pair("AT","G"))
        self.assertEqual(complement_pair(frozenset(("A","G"))),frozenset(("T","C")))

    def test_allele_orientation_not_automatically_harmonized(self):
        ref={
            (1,101):{frozenset(("A","G"))},
            (1,102):{frozenset(("A","T"))},
            (1,103):{frozenset(("C","T"))},
            (1,104):{frozenset(("A","C")),frozenset(("G","T"))},
        }
        self.assertEqual(classify_variant(1,101,"G","A",ref),"UNORDERED_ALLELE_SET_EXACT")
        self.assertEqual(classify_variant(1,101,"T","C",ref),"STRAND_COMPLEMENT_UNRESOLVED")
        self.assertEqual(classify_variant(1,102,"A","T",ref),"PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION")
        self.assertEqual(classify_variant(1,103,"G","T",ref),"COORDINATE_PRESENT_ALLELES_DISCORDANT")
        self.assertEqual(classify_variant(1,104,"A","C",ref),"MULTI_ALLELIC_REFERENCE_POSITION_REVIEW")
        self.assertEqual(classify_variant(1,105,"A","G",ref),"COORDINATE_NOT_IN_LD_BIM")
        self.assertEqual(classify_variant(1,101,"AT","G",ref),"UNSUPPORTED_NON_SNP_OR_ALLELE")

    def test_synthetic_tar_bim_and_gwas_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            archive=root/"ref.tgz"
            bim=(b"1 rs1 0 101 A G\n1 rs2 0 102 A T\n"
                 b"2 rs3 0 103 G T\n3 rs4 0 105 A G\n")
            with tarfile.open(archive,"w:gz") as tgz:
                info=tarfile.TarInfo("1000G/test.bim")
                info.size=len(bim)
                tgz.addfile(info,io.BytesIO(bim))
            ref,counts,names=reference_bim_coordinates(archive,{1,2})
            self.assertEqual(counts["selected_BIM_rows"],3)
            self.assertEqual(counts["unique_selected_chr_pos"],3)
            self.assertEqual(len(names),1)
            gwas=root/"gwas.tsv.gz"
            with gzip.open(gwas,"wt",encoding="utf-8") as f:
                f.write("chromosome\tbase_pair_location\teffect_allele\tother_allele\n")
                f.write("1\t101\tG\tA\n")
                f.write("1\t102\tT\tA\n")
                f.write("2\t103\tC\tA\n")
                f.write("2\t200\tA\tG\n")
                f.write("3\t105\tA\tG\n")
            matched,by_chr=classify_gwas(gwas,ref,{1,2})
            self.assertEqual(matched["all_original_GWAS_variant_rows"],5)
            self.assertEqual(matched["selected_chr_GWAS_variant_rows"],4)
            self.assertEqual(matched["UNORDERED_ALLELE_SET_EXACT"],1)
            self.assertEqual(matched["PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION"],1)
            self.assertEqual(matched["STRAND_COMPLEMENT_UNRESOLVED"],1)
            self.assertEqual(matched["COORDINATE_NOT_IN_LD_BIM"],1)
            self.assertEqual(matched["coordinate_level_overlap"],3)
            self.assertEqual(by_chr["1"]["selected_GWAS_rows"],2)
            self.assertEqual(matched["selected_chr_coordinate_overlap_pct"],75.0)

    def test_truncated_gzip_does_not_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"broken.tsv.gz"
            path.write_bytes(b"not-gzip")
            with self.assertRaises((OSError,EOFError)):
                classify_gwas(path,{}, {1})


if __name__=="__main__":
    unittest.main()
