import argparse
import json
from pathlib import Path
import subprocess
import sys


def main():
    if len(sys.argv)>1 and sys.argv[1]=='blueprint':
        from .architecture import main as blueprint_main
        return blueprint_main(sys.argv[2:])
    if len(sys.argv)>1 and sys.argv[1]=='resources':
        from .resources import main as resources_main
        return resources_main(sys.argv[2:])
    from . import science
    from .engine import execute,atomic_json,sha
    from .compile import compile_project
    parser=argparse.ArgumentParser(description='MasterOmics central workflow')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('run');p.add_argument('projects',nargs='+');p.add_argument('--registry',required=True);p.add_argument('--jobs',type=int,default=1);p.add_argument('--plan',action='store_true')
    p=sub.add_parser('dag');p.add_argument('config');p.add_argument('--jobs',type=int,default=1);p.add_argument('--plan',action='store_true')
    p=sub.add_parser('doctor');p.add_argument('project');p.add_argument('--registry',required=True);p.add_argument('--output')
    p=sub.add_parser('inventory');p.add_argument('root');p.add_argument('output')
    for command in ['acquire','normalize','instruments','harmonize','mr','ld','regional_gate','annotate','report','evidence','score','cohort_panel']:
        p=sub.add_parser(command);p.add_argument('args',nargs='+')
    a=parser.parse_args();root=Path(__file__).resolve().parents[1]
    try:
        if a.command=='doctor':
            from .doctor import doctor
            result=doctor(a.project,a.registry)
            if a.output:atomic_json(a.output,result)
            print(json.dumps(result,indent=2));return 1 if result['issues'] else 0
        if a.command=='run':
            if a.jobs<1:raise ValueError('jobs must be positive')
            code=0
            for recipe in a.projects:
                cfg=Path(recipe).resolve().parent/'compiled'/f'{Path(recipe).stem}.dag.json'
                compile_project(recipe,a.registry,cfg)
                code=max(code,execute(cfg,root,a.jobs,a.plan))
            return code
        if a.command=='dag':return execute(a.config,root,a.jobs,a.plan)
        if a.command=='inventory':
            rootpath=Path(a.root)
            if not rootpath.is_dir():raise ValueError('Inventory root missing')
            records=[]
            for p in sorted(rootpath.rglob('*')):
                if p.is_file() and p.suffix in ['.tsv','.gz','.bed','.bim','.fam','.chain','.csv']:
                    records.append({'path':str(p),'size_bytes':p.stat().st_size,'sha256':sha(p) if p.stat().st_size<50*1024*1024 else None})
            atomic_json(a.output,{'root':str(rootpath),'files':records,'large_files':'metadata_only'});return 0
        args=a.args
        if a.command in ['acquire','normalize']:
            getattr(science,a.command)(json.loads(Path(args[0]).read_text()),args[1])
        elif a.command=='instruments':
            u=json.loads(Path(args[0]).read_text());science.instruments(args[1],u,u['exposure_ld'],args[2],args[3])
        elif a.command=='harmonize':
            u=json.loads(Path(args[0]).read_text());science.harmonize(args[1],args[2],u,args[3],args[4])
        elif a.command=='mr':science.mr(*args)
        elif a.command=='score':science.score(*args)
        elif a.command=='cohort_panel':science.cohort_panel(*args)
        elif a.command=='ld':
            u=json.loads(Path(args[0]).read_text());science.regional_ld(args[1],u['exposure_ld'],u['outcome_ld'],args[2],args[3])
        elif a.command=='regional_gate':
            u=json.loads(Path(args[0]).read_text());qc=json.loads(Path(args[2]).read_text())
            if qc['retained']<u['min_regional_snps'] or qc['retained']/min(qc['exposure_region'],qc['outcome_region'])<u.get('min_overlap_fraction',0.5):
                raise ValueError('Insufficient full-region overlap; cannot run coloc')
            subprocess.run(['Rscript',str(root/'masteromics/analysis.R'),'coloc',args[1],args[0],args[3]],check=True)
        elif a.command=='annotate':science.annotate(*args)
        elif a.command=='report':science.report(args[1:],args[0])
        elif a.command=='evidence':science.evidence(args[1:],args[0])
        return 0
    except Exception as e:
        print(f'MasterOmics error: {e}',file=sys.stderr);return 1
if __name__=='__main__':sys.exit(main())
