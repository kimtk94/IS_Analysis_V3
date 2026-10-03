"""Read-only input/config/environment readiness report; never scientific PASS."""
import importlib.metadata
import json
from pathlib import Path
import re
import shutil
import subprocess
import pandas as pd

def doctor(recipe_path,registry_path):
    recipe=json.loads(Path(recipe_path).read_text());registry=json.loads(Path(registry_path).read_text())
    issues=[];checked=set()
    if not recipe.get('units'):issues.append({'unit':'project','reason':'At least one analysis unit required'})
    def issue(unit,reason):issues.append({'unit':unit,'reason':reason})
    def exists(unit,path):
        if not path or not Path(path).is_file():issue(unit,'Missing file: '+str(path));return False
        return True
    for unit in recipe.get('units',[]):
        uid=unit.get('id','UNKNOWN')
        if not (isinstance(unit.get('start'),int) and isinstance(unit.get('end'),int) and 0<unit['start']<=unit['end']):issue(uid,'Positive integer cis interval required')
        if str(unit.get('chr')) not in [str(i) for i in range(1,23)]+['X','Y']:issue(uid,'Validated chromosome required')
        for role in ['exposure','outcome']:
            did=unit.get(role);spec=registry.get('datasets',{}).get(did)
            if spec is None:issue(uid,'Unregistered dataset: '+str(did));continue
            if spec.get('build')!=unit.get('build'):issue(uid,'Dataset/unit genome build mismatch: '+str(did))
            ld=unit.get(role+'_ld',{})
            if ld.get('ancestry')!=spec.get('ancestry') or ld.get('build')!=spec.get('build'):issue(uid,'LD ancestry/build metadata mismatch: '+role)
            for field in ['matrix','variants']:exists(uid,ld.get(field))
            if did in checked:continue
            checked.add(did)
            if not re.fullmatch('[0-9a-fA-F]{64}',spec.get('sha256','')):issue(did,'Pin the actual source SHA256')
            if spec.get('type')=='quant' and not (isinstance(spec.get('sdY'),(float,int)) and spec['sdY']>0):issue(did,'Explicit trait SD required')
            if spec.get('type')=='cc' and not (isinstance(spec.get('case_fraction'),(float,int)) and 0<spec['case_fraction']<1):issue(did,'Case fraction in (0,1) required')
            if spec.get('effect_scale')!='beta':issue(did,'Explicit beta/logOR effect scale required')
            columns=spec.get('columns',{})
            required={'chr','pos','effect_allele','other_allele','beta','se','pval','eaf','n'}
            if not required.issubset(columns):issue(did,'Missing canonical adapter mappings')
            if exists(did,spec.get('path')):
                try:
                    # Read only the header; do not hash/download/inspect full raw data here.
                    source=pd.read_csv(spec['path'],sep=spec.get('separator','\t'),nrows=0)
                    if not set(columns.values()).issubset(source.columns):issue(did,'Mapped source columns absent')
                except Exception as e:issue(did,'Header read failed: '+str(e))
    if recipe.get('annotation'):exists('annotation',recipe['annotation'])
    if recipe.get('cohort'):
        c=recipe['cohort']
        for field in ['panel','dosages','weights','phenotypes','incident_panel']:
            if c.get(field):exists('cohort',c[field])
        if not c.get('pcs'):issue('cohort','Explicit ancestry PCs required')
    packages=['numpy','pandas','scipy'];runtime={}
    for package in packages:
        try:runtime[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:issue('environment','Missing Python package: '+package)
    required_r=[]
    if recipe.get('advanced_mr',True):required_r+=['TwoSampleMR']
    if recipe.get('regional_analysis',True):required_r+=['coloc','susieR','jsonlite']
    if recipe.get('cohort'):required_r+=['lme4','survival']
    if required_r:
        if not shutil.which('Rscript'):issue('environment','Rscript unavailable')
        else:
            expression='for(p in c('+','.join(json.dumps(p) for p in sorted(set(required_r)))+')) {if(!requireNamespace(p,quietly=TRUE)) stop(p); cat(p,as.character(packageVersion(p)),"\\n")}'
            try:runtime['R']=subprocess.check_output(['Rscript','-e',expression],text=True,stderr=subprocess.STDOUT,timeout=60)
            except (OSError,subprocess.SubprocessError) as e:issue('environment','R dependency check failed: '+str(e))
    return {'project':recipe.get('project'),'status':'BLOCKED_INPUTS' if issues else 'READY_FOR_SCIENTIFIC_GATES',
            'scientific_validation':'NOT_EXECUTED','datasets_checked':len(checked),'runtime':runtime,'issues':issues}
