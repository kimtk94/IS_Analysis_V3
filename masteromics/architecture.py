"""Data-independent shared research blueprint and strict adapter binding compiler.

A blueprint is not an analysis result. Unbound stages block execution; no fake
scientific tables or placeholder success checkpoints are ever produced.
"""
import argparse
import copy
import json
from pathlib import Path
import sys
from .engine import atomic_json,check_graph,execute

# CKD remains pQTL/MR-first. IS is the current BBJ locus-first thesis design.
# Both catalogs compile to the same engine and output-contract machinery.
CKD_CATALOG=[
 ('acquisition',[], 'raw_sources','Local/archive/download adapters; explicit source identity'),
 ('source_qc',['acquisition'],'source_integrity','Checksums, consent/access, ancestry and build'),
 ('normalize',['source_qc'],'canonical_statistics','Map source columns once; numeric/allele contracts'),
 ('cis_regions',['normalize'],'full_cis_statistics','Validated coordinates/build and full regional statistics'),
 ('ld_qc',['cis_regions'],'signed_ld','Population-matched LD, variant order/sign, PSD QC'),
 ('instruments',['cis_regions','ld_qc'],'selected_instruments','Significance/F/LD selection, distinct from full regions'),
 ('harmonize',['instruments','cis_regions'],'aligned_statistics','Chromosome/position/alleles/build; palindromic policy'),
 ('discovery_mr',['harmonize'],'discovery_results','Wald/IVW with explicit effect scale and testing family'),
 ('sensitivity',['discovery_mr'],'sensitivity_results','Eligible median/Egger/intercept/Q; no unavailable method success'),
 ('regional_coloc',['harmonize','cis_regions'],'coloc_results','Full overlap, prior sensitivity, trait SD/case fraction'),
 ('finemap',['regional_coloc','ld_qc'],'finemap_results','SuSiE convergence/credible sets; unresolved is not shared'),
 ('replication',['discovery_mr','sensitivity','finemap'],'replication_evidence','Define independence, EUR/EAS and dual-LD boundaries'),
 ('annotation',['replication'],'functional_evidence','Kidney tissue/single-cell/eQTL and source provenance'),
 ('score',['replication'],'individual_scores','Approved genotype QC, allele-aligned weights and ancestry PCs'),
 ('longitudinal',['score'],'longitudinal_results','Repeated phenotype mixed models with declared missingness policy'),
 ('incident',['score'],'incident_results','Exclude prevalent disease; endpoint definition and Cox diagnostics'),
 ('evidence',['annotation','longitudinal','incident'],'final_evidence','All requested stage coverage, provenance and unresolved findings'),
 ('report',['evidence'],'research_report','Evidence-backed tables/report; completion is not biological validity')]

IS_CATALOG=[
 ('acquisition',[], 'raw_sources','Pinned BBJ/GIGASTROKE and functional evidence source identity'),
 ('source_qc',['acquisition'],'source_integrity','Checksums, ancestry/build, sample and schema audit'),
 ('normalize',['source_qc'],'canonical_gwas','Canonical BBJ/GIGASTROKE GWAS statistics without forced pQTL semantics'),
 ('gwas_loci',['normalize'],'discovery_loci','Primary BBJ locus definition and exclusion rules'),
 ('finemap',['gwas_loci'],'finemap_results','BBJ locus fine-mapping and credible-set evidence'),
 ('cross_ancestry',['finemap'],'replication_evidence','GIGASTROKE/EAS direction and credible-set replication'),
 ('molecular_coloc',['finemap'],'molecular_evidence','eQTL/pQTL molecular colocalization for locus mechanism branches'),
 ('mechanism',['cross_ancestry','molecular_coloc'],'mechanism_evidence','FGF5/ALDH2/SH3PXD2A/COL4A2 branch integration'),
 ('celltype',['mechanism'],'celltype_evidence','Mouse/human cell-type localization with explicit inference limits'),
 ('human_annotation',['celltype'],'human_celltype_evidence','Author-validated GSE256493 vascular annotation'),
 ('evidence',['cross_ancestry','molecular_coloc','celltype','human_annotation'],'final_evidence','Frozen locus-first evidence matrix with unresolved findings'),
 ('report',['evidence'],'research_report','Evidence-backed IS thesis report and manuscript figures')]

CATALOGS={'ckd':CKD_CATALOG,'ischemic_stroke':IS_CATALOG}
OPTIONAL_BY_PROJECT={'ckd':{'score','longitudinal','incident'},'ischemic_stroke':set()}
CANONICAL=['dataset_id','trait_id','gene_symbol','protein_id','assay_id','ancestry','build','chr','pos','variant_id','rsid','effect_allele','other_allele','eaf','beta','se','pval','n','source_file','source_sha256']


def catalog_for(project):
    try:return CATALOGS[project]
    except KeyError:raise ValueError('Project must be CKD or ischemic_stroke')


def optional_for(project):
    catalog_for(project)
    return OPTIONAL_BY_PROJECT[project]


def template(project,workspace):
    catalog=catalog_for(project);optional=optional_for(project)
    policy={'replication_definition':None,'genome_build':None,'palindromic_policy':None,'testing_family':None}
    if project=='ckd':
        policy.update(discovery_ancestry='EUR',outcome_validation_ancestry='EAS',cis_window_bp=None)
    else:
        policy.update(discovery_ancestry='Japanese',outcome_validation_ancestry='EAS',locus_definition=None)
    return {'schema_version':1,'project':project,'workspace':str(Path(workspace).resolve()),
        'data_review_status':'DEFERRED','scope':'structure_only_until_adapters_and_data_are_bound',
        'dataset_registry':{},'analysis_policy':policy,
        'stages':[{'id':sid,'depends':deps,'enabled':sid not in optional,'required':sid not in optional,
                   'artifact_type':artifact,'description':description,'binding':None}
                  for sid,deps,artifact,description in catalog]}


def validate(cfg):
    if cfg.get('schema_version')!=1:raise ValueError('Unsupported blueprint version')
    project=cfg.get('project');catalog=catalog_for(project);optional=optional_for(project)
    workspace=Path(cfg['workspace'])
    if not workspace.is_absolute():raise ValueError('Workspace must be absolute')
    stages=cfg['stages'];expected={x[0] for x in catalog}
    ids=[x['id'] for x in stages]
    if len(ids)!=len(set(ids)) or set(ids)!=expected:raise ValueError('Exactly one of each catalog stage is required')
    seen=set()
    for stage in stages:
        sid=stage['id']
        if not set(stage['depends']).issubset(seen):raise ValueError('Invalid stage dependency order: '+sid)
        prescribed=next(x[1] for x in catalog if x[0]==sid)
        if set(stage['depends'])!=set(prescribed):raise ValueError('Scientific dependency contract changed: '+sid)
        if sid not in optional and not stage.get('enabled'):raise ValueError('Required stage disabled: '+sid)
        seen.add(sid)
    enabled={x['id'] for x in stages if x['enabled']}
    for stage in stages:
        if stage['enabled'] and stage['id'] in optional and not set(stage['depends']).issubset(enabled):
            raise ValueError('Optional cohort stage requires enabled upstream: '+stage['id'])
    return stages


def inspect(cfg):
    stages=validate(cfg);records=[]
    for s in stages:
        status='NOT_REQUESTED' if not s['enabled'] else 'UNBOUND_ADAPTER' if not s.get('binding') else 'BOUND_NOT_DATA_VALIDATED'
        records.append({'id':s['id'],'depends':s['depends'],'status':status,'artifact_type':s['artifact_type']})
    missing=[s['id'] for s in records if s['status']=='UNBOUND_ADAPTER']
    return {'project':cfg['project'],'structure_status':'VALID','execution_status':'BLOCKED_UNBOUND_ADAPTERS' if missing else 'BINDINGS_PRESENT_DATA_NOT_VALIDATED','scientific_status':'NOT_EXECUTED','data_review_status':cfg.get('data_review_status','UNREVIEWED'),'unbound_stages':missing,'stages':records}


def compile_bindings(cfg):
    stages=validate(cfg);assessment=inspect(cfg)
    if assessment['unbound_stages']:raise ValueError('Execution blocked by unbound stages: '+', '.join(assessment['unbound_stages']))
    policy=cfg.get('analysis_policy',{})
    if policy.get('genome_build') not in ['GRCh37','GRCh38']:
        raise ValueError('Explicit supported genome build required')
    if not all(policy.get(k) for k in ['replication_definition','palindromic_policy','testing_family']):
        raise ValueError('Explicit replication, allele and testing-family policies required')
    if cfg['project']=='ckd':
        if not isinstance(policy.get('cis_window_bp'),int) or policy['cis_window_bp']<=0:
            raise ValueError('CKD requires a positive cis window')
    elif not policy.get('locus_definition'):
        raise ValueError('Ischemic stroke requires an explicit locus definition')
    active={s['id'] for s in stages if s['enabled']};dag=[];used_outputs=set();workspace=Path(cfg['workspace']).resolve()
    for stage in stages:
        if not stage['enabled']:continue
        b=stage['binding']
        if not isinstance(b.get('argv'),list) or not b['argv'] or not all(isinstance(x,str) for x in b['argv']):raise ValueError('Binding needs explicit argv: '+stage['id'])
        if not isinstance(b.get('inputs'),list) or not b['inputs'] or not all(isinstance(x,str) for x in b['inputs']):raise ValueError('Binding needs explicit input provenance: '+stage['id'])
        outputs=b.get('outputs',[])
        if not outputs:raise ValueError('Binding needs output contracts: '+stage['id'])
        for o in outputs:
            path=Path(o['path'])
            if not path.is_absolute() or workspace not in path.resolve().parents:raise ValueError('Outputs must be within project workspace')
            resolved=str(path.resolve())
            if resolved in used_outputs:raise ValueError('Duplicate output owner: '+resolved)
            used_outputs.add(resolved)
            if o.get('kind') not in ['tsv','json']:raise ValueError('Only validated TSV/JSON contracts supported')
            if o['kind']=='tsv' and (not o.get('columns') or o.get('min_rows',1)<1):raise ValueError('Nonempty explicit TSV schema required')
            if o['kind']=='json' and not o.get('keys'):raise ValueError('Explicit JSON contract keys required')
        if not all('@out'+str(i) in b['argv'] for i in range(len(outputs))):raise ValueError('Adapter must write engine temporary @out placeholders')
        dag.append({**copy.deepcopy(b),'id':stage['id'],'depends':[d for d in stage['depends'] if d in active]})
    check_graph(dag)
    return {'project':cfg['project'],'output_root':str(workspace/'results'),'stages':dag}


def initialize(root):
    root=Path(root).resolve()
    paths=[root/project/'project.json' for project in ['ckd','ischemic_stroke']]
    if any(p.exists() for p in paths):raise ValueError('Project config exists; initialization never overwrites it')
    for project,path in zip(['ckd','ischemic_stroke'],paths):
        for name in ['config','raw','reference','interim','standardized','instruments','regional','results','state','logs','reports']:(path.parent/name).mkdir(parents=True,exist_ok=True)
        cfg=template(project,path.parent);atomic_json(path,cfg)
        atomic_json(path.parent/'config/dataset_registry.json',{'schema_version':1,'datasets':{},'status':'DATA_REVIEW_DEFERRED'})
        atomic_json(path.parent/'config/schema_contract.json',{'schema_version':1,'canonical_summary_statistics':CANONICAL,'full_regions_separate_from_instruments':True})
        atomic_json(path.parent/'state/structure_status.json',inspect(cfg))
    return [str(p) for p in paths]


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action',required=True)
    p=sub.add_parser('init');p.add_argument('--root',required=True,type=Path)
    p=sub.add_parser('inspect');p.add_argument('config',type=Path)
    p=sub.add_parser('run');p.add_argument('config',type=Path);p.add_argument('--jobs',type=int,default=1);p.add_argument('--plan',action='store_true')
    args=parser.parse_args(argv)
    try:
        if args.action=='init':print(json.dumps({'created':initialize(args.root),'scientific_status':'NOT_EXECUTED'},indent=2));return 0
        cfg=json.loads(args.config.read_text())
        if args.action=='inspect' or args.plan:
            print(json.dumps(inspect(cfg),indent=2));return 0
        dag=compile_bindings(cfg)
        if args.jobs<1:raise ValueError('jobs must be positive')
        compiled=Path(cfg['workspace'])/'config/compiled.dag.json';atomic_json(compiled,dag)
        return execute(compiled,Path(__file__).resolve().parents[1],args.jobs)
    except (ValueError,KeyError,OSError) as e:print('Blueprint error: '+str(e),file=sys.stderr);return 1
if __name__=='__main__':sys.exit(main())
