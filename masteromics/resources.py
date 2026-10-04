"""Offline resource catalog validation and explicit acquisition planning.

Catalog entries describe access, not validated analytical inputs. This module never
fetches files, signs in, or interprets catalog availability as scientific validity.
"""
import argparse
from datetime import date
import json
from pathlib import Path
import re
import sys
from urllib.parse import urlsplit

DEFAULT = Path(__file__).resolve().parents[1] / 'data/metadata/masteromics_resource_registry.json'
CATEGORIES = {'gwas', 'ld_reference', 'human_stroke_scrna', 'human_reference_scrna',
              'animal_stroke_scrna', 'eqtl', 'pqtl', 'sqtl', 'atac', 'spatial', 'proteomics', 'annotation'}
ACCESS = {'open', 'registration', 'controlled', 'publication_only', 'unknown'}
FIELDS = {'dataset_id', 'resource_name', 'category', 'species', 'condition', 'tissue',
          'modality', 'accession', 'source_url', 'publication_doi', 'publication_pmid',
          'genome_build', 'ancestry', 'case_n', 'control_n', 'donor_n', 'sample_n',
          'cell_n', 'raw_available', 'processed_available', 'access_type',
          'download_status', 'checksum_status', 'license_or_dua', 'primary_role',
          'limitations', 'last_verified', 'projects', 'priority', 'verification_scope',
          'artifact_url', 'file_name', 'release'}
GENETIC = {'gwas', 'ld_reference', 'eqtl', 'pqtl', 'sqtl'}


def web_url(value):
    if not isinstance(value, str):
        return False
    try:
        parsed = urlsplit(value)
        return parsed.scheme == 'https' and bool(parsed.hostname) and not parsed.username and not parsed.password
    except ValueError:
        return False


def validate(catalog):
    if catalog.get('schema_version') != 1 or not isinstance(catalog.get('resources'), list) or not catalog['resources']:
        raise ValueError('Nonempty schema_version=1 resource catalog required')
    seen = set()
    for r in catalog['resources']:
        if not FIELDS <= r.keys():
            raise ValueError(f"Missing catalog fields: {sorted(FIELDS-r.keys())}")
        key = r['dataset_id']
        if not isinstance(key, str) or not re.fullmatch(r'[a-z0-9][a-z0-9_-]*', key) or key in seen:
            raise ValueError('Resource IDs must be unique safe identifiers')
        seen.add(key)
        if r['category'] not in CATEGORIES or r['access_type'] not in ACCESS:
            raise ValueError(f'{key}: invalid category/access')
        if not web_url(r['source_url']) or (r['artifact_url'] is not None and not web_url(r['artifact_url'])):
            raise ValueError(f'{key}: invalid HTTPS URL')
        if bool(r['artifact_url']) != bool(r['file_name']):
            raise ValueError(f'{key}: artifact URL and file name must be paired')
        if r['file_name'] and (Path(r['file_name']).name != r['file_name'] or r['file_name'] in {'.', '..'}):
            raise ValueError(f'{key}: unsafe file name')
        if r['download_status'] not in {'not_checked', 'available', 'partial', 'restricted', 'missing', 'failed'}:
            raise ValueError(f'{key}: invalid download status')
        if r['checksum_status'] != 'not_checked':
            raise ValueError(f'{key}: catalog does not certify downloaded checksums')
        if r['genome_build'] not in {None, 'GRCh37', 'GRCh38'}:
            raise ValueError(f'{key}: invalid build')
        if r['verification_scope'] not in {'official_portal_review', 'existing_manifest_not_probed'}:
            raise ValueError(f'{key}: invalid verification scope')
        date.fromisoformat(r['last_verified'])
        if not r['projects'] or not set(r['projects']) <= {'ckd', 'ischemic_stroke'} or r['priority'] not in {1, 2, 3}:
            raise ValueError(f'{key}: invalid projects/priority')
        for field in ['case_n', 'control_n', 'donor_n', 'sample_n', 'cell_n']:
            if r[field] is not None and (type(r[field]) is not int or r[field] < 1):
                raise ValueError(f'{key}: invalid {field}')
        for field in ['raw_available', 'processed_available']:
            if r[field] is not None and type(r[field]) is not bool:
                raise ValueError(f'{key}: invalid availability flag')
        if r['category'] in {'human_stroke_scrna', 'human_reference_scrna'} and r['species'] != 'Homo sapiens':
            raise ValueError(f'{key}: human category requires human species')
        if r['category'] == 'animal_stroke_scrna' and r['species'] == 'Homo sapiens':
            raise ValueError(f'{key}: animal category cannot represent human disease')
    return catalog


def plan(catalog, ids, bindings=None):
    validate(catalog)
    if not ids or len(ids) != len(set(ids)):
        raise ValueError('Select explicit unique resource IDs')
    bindings = bindings or {}
    indexed = {r['dataset_id']: r for r in catalog['resources']}
    if set(ids) - indexed.keys():
        raise ValueError(f'Unknown resource IDs: {sorted(set(ids)-indexed.keys())}')
    if set(bindings) - set(ids):
        raise ValueError('Bindings must belong to the explicitly selected resources')
    items = []
    for key in ids:
        r = indexed[key]
        b = bindings.get(key, {})
        if not isinstance(b, dict):
            raise ValueError(f'{key}: binding must be an object')
        artifact = {field: b.get(field, r.get(field)) for field in ['artifact_url', 'file_name', 'release', 'genome_build', 'ancestry']}
        blocks = []
        # A binding must never override resource access class.
        if r['access_type'] != 'open':
            blocks.append({'registration': 'AUTHENTICATION_REQUIRED', 'controlled': 'DATA_ACCESS_APPROVAL_REQUIRED',
                           'publication_only': 'ARTIFACT_UNVERIFIED', 'unknown': 'ACCESS_UNVERIFIED'}[r['access_type']])
        if not artifact['artifact_url'] or not artifact['file_name']:
            blocks.append('SELECT_EXACT_ARTIFACT')
        elif not web_url(artifact['artifact_url']) or Path(artifact['file_name']).name != artifact['file_name'] or artifact['file_name'] in {'.', '..'}:
            raise ValueError(f'{key}: invalid artifact URL/file name')
        if not isinstance(artifact['release'], str) or not artifact['release'].strip():
            blocks.append('PIN_RELEASE')
        checksum = b.get('sha256')
        if not isinstance(checksum, str) or not re.fullmatch(r'[0-9a-f]{64}', checksum):
            blocks.append('PIN_SHA256')
        if b.get('license_reviewed') is not True:
            blocks.append('REVIEW_LICENSE_OR_DUA')
        if r['category'] in GENETIC:
            if artifact['genome_build'] not in {'GRCh37', 'GRCh38'}:
                blocks.append('PIN_GENOME_BUILD')
            if not isinstance(artifact['ancestry'], str) or not artifact['ancestry'].strip():
                blocks.append('PIN_ANCESTRY')
        items.append({'dataset_id': key, 'access_type': r['access_type'], 'status': 'BLOCKED' if blocks else 'READY_FOR_ACQUISITION',
                      'blocks': blocks, 'artifact': artifact, 'sha256': checksum,
                      'scientific_status': 'NOT_REVIEWED', 'independence_status': 'NOT_REVIEWED'})
    return {'mode': 'OFFLINE_PLAN_ONLY', 'downloads_executed': 0, 'items': items}


def main(argv=None):
    parser = argparse.ArgumentParser(description='MasterOmics resource catalog; metadata only, never downloads')
    parser.add_argument('--registry', default=str(DEFAULT))
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('validate')
    p = sub.add_parser('list'); p.add_argument('--project', choices=['ckd', 'ischemic_stroke'])
    p = sub.add_parser('plan'); p.add_argument('--ids', nargs='+', required=True); p.add_argument('--bindings')
    a = parser.parse_args(argv)
    try:
        catalog = validate(json.loads(Path(a.registry).read_text()))
        if a.command == 'validate':
            result = {'status': 'VALID', 'resources': len(catalog['resources']), 'network_requests': 0, 'scientific_status': 'NOT_REVIEWED'}
        elif a.command == 'list':
            result = [r for r in catalog['resources'] if not a.project or a.project in r['projects']]
        else:
            result = plan(catalog, a.ids, json.loads(Path(a.bindings).read_text()) if a.bindings else None)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        # Blocked acquisition plans are actionable failures, not successful downloads.
        return 2 if a.command == 'plan' and any(i['blocks'] for i in result['items']) else 0
    except (ValueError, TypeError, KeyError, OSError) as error:
        print(f'Resource catalog error: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
