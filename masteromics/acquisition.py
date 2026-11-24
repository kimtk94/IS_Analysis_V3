"""Checksum-pinned transfers, resumable staging and per-artifact audit history.

No credentials or signup flows. Scientific input validation is a later stage.
"""
from contextlib import contextmanager
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import time
import urllib.error
import urllib.request

from .resources import web_url

CHUNK = 1024 * 1024


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(CHUNK), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    fd, name = tempfile.mkstemp(prefix=path.name+'.', suffix='.tmp', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(value, stream, ensure_ascii=False, indent=2)
            stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def now():
    return datetime.now(timezone.utc).isoformat()


@contextmanager
def lock(path):
    # Never guess whether a previous process is dead, especially across hosts.
    fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write(json.dumps({'pid': os.getpid(), 'started_at': now()}))
        yield
    finally:
        Path(path).unlink(missing_ok=True)


def verify(path, expected, expected_bytes=None, compressed=False):
    size = Path(path).stat().st_size
    if size <= 0:
        raise ValueError('Empty artifact')
    if expected_bytes is not None and size != expected_bytes:
        raise ValueError('Artifact byte count mismatch')
    if digest(path) != expected:
        raise ValueError('Artifact SHA256 mismatch')
    if compressed:
        count = 0
        with gzip.open(path, 'rb') as stream:
            for block in iter(lambda: stream.read(CHUNK), b''):
                count += len(block)
        if count == 0:
            raise ValueError('Empty gzip payload')
    return size


def transfer(url, partial, metadata, expected_bytes, timeout, opener):
    offset = partial.stat().st_size if partial.exists() else 0
    headers = {'Accept-Encoding': 'identity', 'User-Agent': 'MasterOmics/1 acquisition'}
    if offset:
        headers['Range'] = f'bytes={offset}-'
        if metadata.get('etag'):
            headers['If-Range'] = metadata['etag']
    request = urllib.request.Request(url, headers=headers)
    with opener(request, timeout=timeout) as response:
        if not web_url(response.geturl()):
            raise ValueError('HTTPS redirect required')
        status = response.getcode()
        if response.headers.get('Content-Encoding', 'identity') != 'identity':
            raise ValueError('Encoded HTTP body cannot represent pinned file bytes')
        length = response.headers.get('Content-Length')
        length = int(length) if length is not None else None
        if length is not None and length < 0:
            raise ValueError('Invalid Content-Length')
        if status == 206:
            match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
            if not offset or not match:
                raise ValueError('Invalid partial response')
            start, end, total = map(int, match.groups())
            if start != offset or end < start or end >= total or (length is not None and length != end-start+1):
                raise ValueError('Content-Range does not match requested offset')
            if expected_bytes is not None and total != expected_bytes:
                raise ValueError('Remote byte count differs from pinned size')
            if metadata.get('etag') and response.headers.get('ETag') not in {None, metadata['etag']}:
                raise ValueError('Entity changed during resume')
            length = end-start+1
            mode = 'ab'
        elif status == 200:
            # Servers may ignore Range or change entity: restart, never append.
            mode = 'wb'; offset = 0
            if expected_bytes is not None and length is not None and length != expected_bytes:
                raise ValueError('Remote byte count differs from pinned size')
        else:
            raise ValueError(f'Unexpected HTTP status {status}')
        etag = response.headers.get('ETag')
        if etag and '\n' not in etag and '\r' not in etag:
            metadata['etag'] = etag
        atomic_json(partial.with_name(partial.name+'.json'), metadata)
        read = 0
        with partial.open(mode) as stream:
            while True:
                block = response.read(CHUNK)
                if not block:
                    break
                stream.write(block); read += len(block)
                if length is not None and read > length:
                    raise ValueError('HTTP body exceeds declared range/length')
            stream.flush(); os.fsync(stream.fileno())
        if length is not None and read != length:
            raise OSError('Interrupted HTTP body; retain partial for resume')
        if status == 206 and end != total-1:
            raise OSError('Partial range ended before complete artifact')


def stage_file(target, sha256, *, url=None, local_path=None, identity=None,
               expected_bytes=None, retries=2, timeout=60, opener=None):
    """Stage one explicitly reviewed source. Existing targets are never overwritten."""
    if not isinstance(sha256, str) or not re.fullmatch(r'[0-9a-f]{64}', sha256):
        raise ValueError('Pinned SHA256 required')
    if url is not None and not web_url(url):
        raise ValueError('Public download requires credential-free HTTPS URL')
    if expected_bytes is not None and (type(expected_bytes) is not int or expected_bytes <= 0):
        raise ValueError('Expected byte count must be positive')
    if type(retries) is not int or not 0 <= retries <= 10 or timeout <= 0:
        raise ValueError('Invalid retry/timeout settings')
    target = Path(target).absolute()
    local = Path(local_path).resolve() if local_path else None
    if local is not None and url is not None:
        raise ValueError('Choose either local or remote source')
    target.parent.mkdir(parents=True, exist_ok=True)
    compressed = target.name.endswith(('.gz', '.bgz'))
    spec = {'target': str(target), 'source': str(local) if local else url, 'sha256': sha256,
            'expected_bytes': expected_bytes, 'identity': identity or {}}
    config_hash = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()
    record = target.with_name(target.name+'.acquisition.json')
    history = target.with_name(target.name+'.history.jsonl')
    partial = target.with_name(target.name+'.partial')
    partial_meta = partial.with_name(partial.name+'.json')
    with lock(target.with_name(target.name+'.lock')):
        if record.is_symlink() or history.is_symlink():
            raise ValueError('Audit paths cannot be symlinks')
        previous = json.loads(record.read_text()) if record.is_file() else {}
        state = {'status': 'RUNNING', 'started_at': now(), 'spec': spec, 'config_sha256': config_hash,
                 'code_sha256': digest(__file__), 'scientific_status': 'NOT_REVIEWED'}
        atomic_json(record, state)
        try:
            if target.is_symlink():
                raise ValueError('Artifact target cannot be a symlink')
            if target.exists():
                if local and local != target.resolve():
                    verify(local, sha256, expected_bytes, compressed)
                size = verify(target, sha256, expected_bytes, compressed)
                same = (previous.get('config_sha256') == config_hash and previous.get('code_sha256') == state['code_sha256']
                        and previous.get('status') in {'SUCCESS_VERIFIED', 'VERIFIED_EXISTING', 'SKIPPED_EXISTING_VALID'})
                state.update(status='SKIPPED_EXISTING_VALID' if same else 'VERIFIED_EXISTING', size_bytes=size)
            else:
                if local is None and url is None:
                    raise ValueError('No local file or public URL available')
                if partial.is_symlink() or partial_meta.is_symlink():
                    raise ValueError('Partial paths cannot be symlinks')
                if partial.exists():
                    if not partial_meta.is_file() or json.loads(partial_meta.read_text()).get('config_sha256') != config_hash:
                        raise ValueError('Partial identity changed; choose a new target or review partial files')
                metadata = {'config_sha256': config_hash}
                if partial_meta.exists():
                    metadata = json.loads(partial_meta.read_text())
                atomic_json(partial_meta, metadata)
                if local is not None:
                    verify(local, sha256, expected_bytes, compressed)
                    with local.open('rb') as src, partial.open('wb') as dst:
                        shutil.copyfileobj(src, dst, CHUNK)
                        dst.flush(); os.fsync(dst.fileno())
                else:
                    complete = partial.exists() and partial.stat().st_size > 0 and digest(partial) == sha256
                    if not complete:
                        for attempt in range(retries+1):
                            try:
                                transfer(url, partial, metadata, expected_bytes, timeout, opener or urllib.request.urlopen)
                                break
                            except urllib.error.HTTPError as error:
                                if error.code not in {408, 429, 500, 502, 503, 504} or attempt == retries:
                                    raise
                                time.sleep(min(attempt+1, 3))
                            except (OSError, EOFError) as error:
                                if attempt == retries:
                                    raise
                                time.sleep(min(attempt+1, 3))
                try:
                    size = verify(partial, sha256, expected_bytes, compressed)
                except (ValueError, OSError, EOFError):
                    partial.unlink(missing_ok=True); partial_meta.unlink(missing_ok=True)
                    raise
                # A concurrent non-cooperating writer must not be overwritten either.
                os.link(partial, target)
                partial.unlink(); partial_meta.unlink(missing_ok=True)
                state.update(status='SUCCESS_VERIFIED', size_bytes=size)
        except Exception as error:
            state.update(status='FAILED', error_type=type(error).__name__, error=str(error)[:1000])
            raise
        finally:
            state['finished_at'] = now()
            atomic_json(record, state)
            with history.open('a') as stream:
                stream.write(json.dumps(state, ensure_ascii=False)+'\n')
                stream.flush(); os.fsync(stream.fileno())
    return state


def collect(catalog, ids, bindings, root, **transfer_options):
    from .resources import plan
    assessment = plan(catalog, ids, bindings)
    if any(item['blocks'] for item in assessment['items']):
        return {'status': 'BLOCKED', 'plan': assessment, 'artifacts': []}
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    with lock(root/'.collection.lock'):
        result = {'status': 'RUNNING', 'started_at': now(), 'artifacts': [],
                  'catalog_sha256': hashlib.sha256(json.dumps(catalog, sort_keys=True).encode()).hexdigest(),
                  'scientific_status': 'NOT_REVIEWED', 'selected_ids': ids}
        atomic_json(root/'collection_summary.json', result)
        for item in assessment['items']:
            key = item['dataset_id']; artifact = item['artifact']
            target = root/'raw'/key/artifact['file_name']
            identity = {'dataset_id': key, 'access_type': item['access_type'],
                        'release': artifact['release'], 'genome_build': artifact['genome_build'],
                        'ancestry': artifact['ancestry'], 'authorization_ref': item.get('authorization_ref')}
            try:
                if not target.resolve().is_relative_to(root):
                    raise ValueError('Artifact path escapes collection root')
                record = stage_file(target, item['sha256'],
                                    url=artifact['artifact_url'] if item['source_mode']=='remote' else None,
                                    local_path=item.get('local_path'), identity=identity,
                                    expected_bytes=bindings[key].get('expected_bytes'), **transfer_options)
                result['artifacts'].append({'dataset_id': key, 'status': record['status'],
                                            'path': str(target), 'sha256': item['sha256'], 'identity': identity,
                                            'size_bytes': record['size_bytes']})
            except Exception as error:
                result['artifacts'].append({'dataset_id': key, 'status': 'FAILED', 'error_type': type(error).__name__, 'error': str(error)[:1000]})
            atomic_json(root/'collection_summary.json', result)
        result['status'] = 'FAILED' if any(i['status']=='FAILED' for i in result['artifacts']) else 'SUCCESS_VERIFIED'
        result['finished_at'] = now()
        atomic_json(root/'collection_summary.json', result)
    return result
