"""Native BBJ allele/chain G1 strict input audit synthetic fixtures."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

SRC=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_g1_bbj_native_allele_chain.py"
spec=importlib.util.spec_from_file_location("bbj_native_allele_chain",SRC)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class LO:
    def __init__(self,coords=None):
        self.coords=coords if coords is not None else [('chr4',79759932,'+',1)]
    def convert_coordinate(self,c,p):
        return self.coords

def synthetic():
    r=dict(locus='BBJ_IS_L001',dataset_key='GTEx_V8__Artery_Aorta',
        gene_base='ENSG00000152784',variant_id='4:80681087:G:C',
        match_key='4:79759933:G:C',gwas_beta='0.0089170152351476',
        gwas_eaf='0.770149969929239')
    n={'v':'4:80681087:G:C','CHR':'4','POS':'80681087',
       'Allele1':'G','Allele2':'C','BETA':r['gwas_beta'],
       'AF_Allele2':r['gwas_eaf']}
    return r,n

class TestNativeBBJ(unittest.TestCase):
    def test_exact_native_and_chain_pass(self):
        r,n=synthetic()
        z=mod.score(r,n,LO())
        self.assertEqual(z['status'],'PASS_NATIVE_EFFECT_ALLELE_PLUS_CHAIN')
        self.assertEqual(z['gwas_effect_vs_id'],'ALT')

    def test_beta_mismatch_fail(self):
        r,n=synthetic();n['BETA']='-0.11'
        z=mod.score(r,n,LO())
        self.assertEqual(z['status'],'FAIL')
        self.assertIn('BETA_DIVERGES',z['reason'])

    def test_eaf_mismatch_fail(self):
        r,n=synthetic();n['AF_Allele2']='0.12'
        z=mod.score(r,n,LO())
        self.assertEqual(z['status'],'FAIL')
        self.assertIn('EAF_DIVERGES',z['reason'])

    def test_missing_native_never_pass(self):
        r,n=synthetic()
        self.assertEqual(mod.score(r,None,LO())['status'],'NOT_FOUND_NATIVE')

    def test_multiple_liftover_positions_fail(self):
        r,n=synthetic()
        z=mod.score(r,n,LO([('chr4',79759932,'+',1),('chr4',79759933,'+',2)]))
        self.assertEqual(z['status'],'FAIL')

    def test_genome_build_mismatch_fail(self):
        r,n=synthetic()
        z=mod.score(r,n,LO([('chr4',8000000,'+',1)]))
        self.assertEqual(z['status'],'FAIL')
        self.assertIn('CHAIN_GRCH38',z['reason'])

    def test_allele_set_swap_sources_ref_does_not_fail(self):
        r,n=synthetic()
        # Where Allele2 is REF, the input GWAS beta is defined on REF.
        n['Allele1'],n['Allele2']=n['Allele2'],n['Allele1']
        z=mod.score(r,n,LO())
        self.assertEqual(z['status'],'PASS_NATIVE_EFFECT_ALLELE_PLUS_CHAIN')
        self.assertEqual(z['gwas_effect_vs_id'],'REF')

    def test_existing_output_fails_closed(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'out';p.mkdir()
            sentinel=p/'preserved.txt';sentinel.write_text('ORIGINAL')
            with self.assertRaises(FileExistsError):
                mod.run(Path(t)/'input',Path(t)/'archive',
                        Path(t)/'chain',Path(t)/'lib',p)
            self.assertEqual(sentinel.read_text(),'ORIGINAL')

if __name__=='__main__':
    unittest.main()
