#!/usr/bin/env python3
"""Network-backed, bounded P0 ocular public-data audit. No controlled/genotype downloads."""
import argparse, gzip, hashlib, io, json, re, time, urllib.parse, urllib.request, zipfile
from collections import defaultdict
from html.parser import HTMLParser
from pathlib import Path
from datetime import datetime, timezone

LIMIT=16*1024*1024
class Links(HTMLParser):
    def __init__(self): super().__init__(); self.urls=[]
    def handle_starttag(self,t,a):
        if t=='a':
            h=dict(a).get('href','')
            if h and not h.startswith('?') and not h.endswith('/'): self.urls.append(h)
def fetch(url,limit=LIMIT):
    err=None
    for attempt in range(3):
        try:
            q=urllib.request.Request(url,headers={'User-Agent':'ocular-p0-public-audit/1.0'})
            with urllib.request.urlopen(q,timeout=35) as r:
                b=r.read(limit+1)
                if len(b)>limit: raise ValueError('EXCEEDS_DOWNLOAD_CAP')
                if b.lstrip().startswith(b'<!DOCTYPE html') and url.endswith(('.gz','.docx')): raise ValueError('HTML_INSTEAD_OF_DATA')
                return b,{'http':r.status,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
        except Exception as e:
            err=f'{type(e).__name__}: {str(e)[:170]}'
            if attempt<2:time.sleep(attempt+1)
    raise RuntimeError(err)
def root(acc): return f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/'
def soft_samples(raw):
    out=[]; row=None
    for line in gzip.decompress(raw).decode('utf8','replace').splitlines():
        if line.startswith('^SAMPLE ='):
            if row:out.append(row)
            row={'gsm':line.partition('=')[2].strip(),'title':'','desc':[],'links':[]}
        elif row and line.startswith('!Sample_') and ' = ' in line:
            k,v=line.split(' = ',1)
            if k=='!Sample_title':row['title']=v
            if k=='!Sample_description':row['desc'].append(v)
            if k=='!Sample_supplementary_file':row['links'].append(v)
    if row:out.append(row)
    return out
def conditions(samples,acc):
    donor=defaultdict(set); count=defaultdict(int); per=[]
    for s in samples:
        text=s['title']
        low=text.lower()
        cond='DEX' if re.search(r'\b(dexamethasone|dex[-_a-z]*)\b',low) else ('CTRL' if re.search(r'\b(vehicle|control|untreated)\b',low) else 'OTHER')
        if acc=='GSE313736': ids=set(re.findall(r'\bHTMC[- ]?(\d+)\b',' '.join(s['desc']),re.I))
        else: ids=set(re.findall(r'\bTM\s*(\d+)\b',s['title'],re.I))
        donor[cond].update(ids);count[cond]+=1
        per.append({'gsm':s['gsm'],'condition':cond,'donor_count_in_sample_description':len(ids),'title':s['title']})
    return {'gsm_total':len(samples),'samples_by_condition':dict(count),'distinct_donors_by_condition':{k:len(v) for k,v in donor.items()},'donor_overlap_ctrl_dex':len(donor['CTRL']&donor['DEX']),'donor_overlap_is_not_cell_level_demultiplexing':True,'sample_preview':per[:5]}
def metadata_summary(blob,name):
    """Only aggregate donor statistics. Never publish donor identifiers or genotypes."""
    import csv
    with gzip.open(io.BytesIO(blob),'rt',encoding='utf8',errors='replace') as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        cols=reader.fieldnames or []
        id_field=next((c for c in ['donor_id_demuxlet','donor_id','sampleid'] if c in cols),None)
        cond_field=next((c for c in ['DEX Treatment','dex_treatment','treatment','Treatment'] if c in cols),None)
        count=0; ids=set(); by_cond=defaultdict(set); sample_rows=defaultdict(int); missing=0
        for r in reader:
            count+=1
            key=(r.get(id_field) or '').strip() if id_field else ''
            cond=(r.get(cond_field) or '').strip() if cond_field else ''
            sample_rows[cond]+=1
            if key and key.lower() not in ('na','nan','none','.','unassigned','doublet'):
                ids.add(key);by_cond[cond].add(key)
            else:missing+=1
            if count>=500000:break
    treatment_values=dict(sorted(sample_rows.items(),key=lambda x:-x[1])[:12])
    out={'rows_checked':count,'columns':cols,'donor_id_col':id_field,
         'donor_unique_count':len(ids),'missing_or_invalid_donor_rows':missing,
         'treatment_col':cond_field,'treatment_counts':treatment_values,
         'unique_donors_per_treatment':{k:len(v) for k,v in by_cond.items()}}
    if cond_field:
        def is_dex(x):return bool(re.search(r'dex|^1$|^yes$|^true$',x,re.I)) and not bool(re.search(r'no[-_ ]?dex|non[-_ ]?dex|without',x,re.I))
        def is_control(x):return bool(re.search(r'ctrl|vehicle|control|false|^0$|no|untreated',x,re.I))
        d=set().union(*(v for k,v in by_cond.items() if is_dex(k))) if by_cond else set()
        c=set().union(*(v for k,v in by_cond.items() if is_control(k))) if by_cond else set()
        out['donor_overlap_dex_control_validated_metadata']=len(d & c)
        out['donor_count_dex']=len(d);out['donor_count_control']=len(c)
        if not d or not c:out['warning']='TREATMENT_VALUES_REQUIRE_MANUAL_REVIEW'
    return out

def probe_gzip(blob,name):
    if blob[:2]!=b'\x1f\x8b':raise ValueError('BAD_GZIP')
    with gzip.open(io.BytesIO(blob),'rt',encoding='utf8',errors='replace') as f:
        first=f.readline().rstrip('\n');second=f.readline().rstrip('\n')
    n=len(first.split('\t'))
    if n<2:raise ValueError('INVALID_TABLE')
    if 'narrowpeak' in name.lower():
        fields=first.split('\t');assert len(fields)>=6
        int(fields[1]);int(fields[2])
    return {'gzip':'PASS','fields_on_first_line':n,'first_fields':first.split('\t')[:6] if 'narrowpeak' not in name.lower() else ['chr','start','end'],'second_row_present':bool(second)}
def geo(acc,patterns):
    r={'accession':acc,'status':'HOLD','soft':{},'supplementary':{},'files':[]}
    samples=[]
    try:
        raw,tx=fetch(root(acc)+'soft/'+acc+'_family.soft.gz',limit=9*1024*1024)
        samples=soft_samples(raw)
        r['soft']={'transfer':tx,'design':conditions(samples,acc)}
    except Exception as e:r['soft']['error']=str(e)
    try:
        raw,tx=fetch(root(acc)+'suppl/',limit=2*1024*1024)
        h=Links();h.feed(raw.decode('utf8','replace'))
        urls=[urllib.parse.urljoin(root(acc)+'suppl/',x) for x in h.urls]
        r['supplementary']={'transfer':tx,'filenames':[urllib.parse.unquote(x.split('/')[-1]) for x in urls]}
    except Exception as e:
        r['supplementary']['error']=str(e);urls=[]
    # Also inspect GSM-attached links from SOFT when series-level listing has no small files.
    for s in samples:
        for l in s['links']:
            if l.startswith('ftp://ftp.ncbi.nlm.nih.gov/'):l='https://ftp.ncbi.nlm.nih.gov/'+l.split('ftp.ncbi.nlm.nih.gov/',1)[1]
            if l.startswith('https://ftp.ncbi.nlm.nih.gov/'):urls.append(l)
    urls=list(dict.fromkeys(urls))
    for pattern in patterns:
        matches=[u for u in urls if re.search(pattern,urllib.parse.unquote(u.split('/')[-1]),re.I)]
        if not matches:r['files'].append({'pattern':pattern,'status':'NOT_FOUND_IN_SOFT_OR_LISTING'})
        for url in matches[:4]:
            fn=urllib.parse.unquote(url.split('/')[-1]);f={'name':fn,'url':url}
            if re.search(r'\.vcf|\.mtx|genotype|\.bam|\.fastq|\.sra',fn,re.I):f['status']='SKIPPED_LARGE_OR_GENOTYPE'
            else:
                try:
                    b,tx=fetch(url)
                    f.update({'status':'TRANSFER_PASS','transfer':tx})
                    if fn.endswith('.gz'):
                        f['content_probe']=probe_gzip(b,fn)
                        if 'cell_metadata' in fn.lower() or 'donor_metadata' in fn.lower():
                            f['metadata_summary']=metadata_summary(b,fn)
                except Exception as e:f.update({'status':'FAIL','reason':str(e)})
            r['files'].append(f)
    if samples and any(x.get('status')=='TRANSFER_PASS' for x in r['files']):r['status']='PARTIAL_G0_PASS'
    return r
def mendeley():
    base='https://api.data.mendeley.com/datasets/nszmnfcwhs'
    r={'accession':'nszmnfcwhs_v1','status':'HOLD','files':[],'api':{}}
    for typ,url in [('dataset',base+'?version=1'),('file_list',base+'/files?version=1')]:
        try:
            b,tx=fetch(url,limit=2*1024*1024)
            j=json.loads(b);r['api'][typ]={'transfer':tx}
            if typ=='file_list':
                items=j if isinstance(j,list) else j.get('results',j.get('files',[]))
                for d in items:
                    r['files'].append({'name':d.get('filename',d.get('name','')),'id':d.get('id'),'size':d.get('size'),'status':'LISTED'})
        except Exception as e:r['api'][typ]={'error':str(e)}
    for f in r['files']:
        if not f['name'].lower().endswith('.docx') or not f.get('id'):continue
        try:
            b,tx=fetch(base+'/files/'+f['id']+'/file_downloaded?version=1',limit=8*1024*1024)
            z=zipfile.ZipFile(io.BytesIO(b));xml=z.read('word/document.xml').decode('utf8','replace')
            f.update({'status':'TRANSFER_PASS','transfer':tx,'docx_valid':True,'tables':xml.count('<w:tbl>')})
        except Exception as e:f.update({'status':'FAIL','reason':str(e)})
    if r['files'] and all(f['status']=='TRANSFER_PASS' for f in r['files'] if f['name'].endswith('.docx')):r['status']='PUBLIC_TABLES_TRANSFERRED_NO_SUBJECT_MATRIX_VERIFIED'
    return r
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);a=ap.parse_args()
    datasets=[
      geo('GSE313736',[r'cell_metadata.*\.gz$',r'donor_metadata.*\.gz$',r'\.dmx\.best\.gz$',r'_barcodes\.tsv\.gz$']),
      geo('GSE301517',[r'Control_vs_Dex_Gene_CPM_.*\.gz$']),
      geo('GSE301516',[r'ATAC_(Control|Dex)\.narrowPeak\.gz$']),
      mendeley(),
      {'accession':'phs004037.v1.p1','status':'HOLD_CONTROLLED','reason':'dbGaP authorized access required; individual genotypes and phenotypes were not fetched'}]
    audit={'run_utc':datetime.now(timezone.utc).isoformat(),'policy':'No private genotype/phenotype, raw FASTQ, or full scRNA MTX transferred. Public audit metadata only.','datasets':datasets}
    p=Path(a.out);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n')
    for d in datasets:print(d['accession'],d['status'],'transferred',sum(f.get('status')=='TRANSFER_PASS' for f in d.get('files',[])))
    print('OUTPUT',p)
if __name__=='__main__':main()
