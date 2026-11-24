import csv,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"run_master_instrument_qc.py"
def write(p,h,rows):
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(h);w.writerows(rows)
def test_instrument_qc(tmp_path):
    src=tmp_path/"norm.tsv";out=tmp_path/"qc.tsv"
    h=["protein_id","gene_symbol","rsid","variant_id","chrom","pos","effect_allele","other_allele","beta","se","eaf","p","n","ancestry","genome_build","platform","cis_trans"]
    write(src,h,[
      ["P1","G1","rs1","","1","100","A","G",0.2,0.02,0.3,1e-8,1000,"EUR","GRCh37","Olink","cis"],
      ["P1","G1","rs2","","1","200","C","T",0.01,0.02,0.3,1e-8,1000,"EUR","GRCh37","Olink","cis"],
      ["P2","G2","rs3","","6","30000000","A","T",0.3,0.02,0.4,1e-9,1000,"EUR","GRCh37","Olink","cis"],
    ])
    p=subprocess.run([sys.executable,str(SCRIPT),"--input",str(src),"--output",str(out),"--exclude-mhc","--min-fstat","10"],cwd=ROOT,capture_output=True,text=True)
    assert p.returncode==0,p.stdout+"\n"+p.stderr
    with out.open(encoding="utf-8") as f:rows=list(csv.DictReader(f,delimiter="\t"))
    assert len(rows)==1
    assert rows[0]["rsid"]=="rs1"
    assert float(rows[0]["f_stat"])>10
