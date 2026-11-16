"""Offline bounded coordinate concordance and hold tests."""
import gzip
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.ckm_stroke.audit_bmi_ldscore_coordinates import (
    original_candidate_coordinates, match_ldscore_stream, read_candidate_ids, audit,
)


class CoordinatesTests(unittest.TestCase):
    def test_original_source_candidate_rsids_matched(self):
        pos={"rs1":None,"rs2":None}
        rows=iter([
            (1,["1","123","rs1","A","C","0.2","0.1","0.03","0.1","100"]),
            (2,["2","345","rs2","G","T","0.2","0.1","0.03","0.1","100"]),
            (0,["__SOURCE_COMPLETE_SHA256_PASS__","456"])
        ])
        counts=original_candidate_coordinates(pos,rows)
        self.assertEqual(pos,{"rs1":(1,123),"rs2":(2,345)})
        self.assertEqual(counts["original_sha256_verified"],1)

    def test_duplicate_original_coordinate_discordant_aborts(self):
        pos={"rs1":None}
        rows=iter([
            (1,["1","123","rs1"]+["A"]*7),
            (2,["1","124","rs1"]+["A"]*7),
            (0,["SENTINEL","200"])
        ])
        with self.assertRaisesRegex(ValueError,"multiple coordinates"):
            original_candidate_coordinates(pos,rows)

    def test_real_chrpos_concordance_and_mismatch(self):
        pos={"rs1":(1,123),"rs2":(1,346),"rs3":(2,200)}
        raw=("CHR SNP BP L2\n1 rs1 123 1.1\n1 rs2 345 1.5\n").encode()
        from collections import Counter
        seen=set()
        c=Counter()
        examples=[]
        import io
        match_ldscore_stream(io.BytesIO(gzip.compress(raw)),1,pos,seen,c,examples)
        self.assertEqual(c["candidate_coord_concordant"],1)
        self.assertEqual(c["candidate_coord_discordant"],1)
        self.assertEqual(examples[0]["LD_BP"],345)
        self.assertEqual(seen,{"rs1","rs2"})

    def test_repeated_ldscore_marker_counted(self):
        pos={"rs1":(1,123)}
        raw=b"CHR SNP BP L2\n1 rs1 123 1.1\n1 rs1 123 1.2\n"
        import io
        from collections import Counter
        c=Counter()
        match_ldscore_stream(io.BytesIO(gzip.compress(raw)),1,pos,set(),c,[])
        self.assertEqual(c["candidate_repeated_ldscore_rsid_rows"],1)

    def test_candidate_file_requires_unique_rsids(self):
        with TemporaryDirectory() as d:
            f=Path(d)/"dup.gz"
            f.write_bytes(gzip.compress(b"SNP\tN\tZ\tA1\tA2\nrs1\t1\t1\tA\tC\nrs1\t2\t1\tA\tG\n"))
            with self.assertRaisesRegex(ValueError,"duplicate"):
                read_candidate_ids(f)

    def test_missing_candidate_fails_closed(self):
        with TemporaryDirectory() as d:
            base=Path(d)
            output=audit(base/"not_found.gz",base/"missing.json")
            self.assertEqual(output["technical_status"],"HOLD_SOURCE_OR_QC_ERROR")
            self.assertFalse(output["ready_for_ldsc"])
            self.assertFalse(output["ldsc_h2_rg_performed"])


if __name__=="__main__":
    unittest.main()
