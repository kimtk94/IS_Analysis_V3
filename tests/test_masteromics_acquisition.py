"""Transfer tests use tiny local fixtures and simulated HTTP, never public DBs."""
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import urllib.error
from unittest.mock import patch

from masteromics.acquisition import collect, stage_file
from masteromics.resources import DEFAULT, plan

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT/'tests/fixtures/gigastroke.tsv'
PAYLOAD = FIXTURE.read_bytes()
SHA = hashlib.sha256(PAYLOAD).hexdigest()
URL = 'https://example.invalid/fixture.tsv'


class Response(io.BytesIO):
    def __init__(self, body, status=200, headers=None):
        super().__init__(body)
        self.status = status
        self.headers = headers if headers is not None else {'Content-Length': str(len(body)), 'ETag': '"fixture-v1"'}
    def getcode(self):
        return self.status
    def geturl(self):
        return URL


def binding(**changes):
    b = {'artifact_url': URL, 'file_name': 'fixture.tsv', 'sha256': SHA,
         'release': 'fixture-v1', 'license_reviewed': True, 'genome_build': 'GRCh38', 'ancestry': 'EUR'}
    b.update(changes)
    return b


class AcquisitionTest(unittest.TestCase):
    def test_collect_local_controlled_audit_resume_and_no_overwrite(self):
        catalog = json.loads(DEFAULT.read_text())
        b = binding(local_path=str(FIXTURE), access_authorized=True, authorization_ref='test-approval-reference')
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)/'collection'
            result = collect(catalog, ['koges_cohort'], {'koges_cohort': b}, root)
            self.assertEqual(result['status'], 'SUCCESS_VERIFIED')
            self.assertEqual(result['artifacts'][0]['status'], 'SUCCESS_VERIFIED')
            target = root/'raw/koges_cohort/fixture.tsv'
            self.assertEqual(target.read_bytes(), PAYLOAD)
            record = json.loads(target.with_name(target.name+'.acquisition.json').read_text())
            self.assertEqual(record['spec']['identity']['access_type'], 'controlled')
            self.assertEqual(record['scientific_status'], 'NOT_REVIEWED')
            self.assertEqual(collect(catalog, ['koges_cohort'], {'koges_cohort': b}, root)['artifacts'][0]['status'], 'SKIPPED_EXISTING_VALID')
            target.write_bytes(b'corrupted')
            result = collect(catalog, ['koges_cohort'], {'koges_cohort': b}, root)
            self.assertEqual(result['status'], 'FAILED')
            self.assertEqual(target.read_bytes(), b'corrupted')
            self.assertEqual(len(target.with_name(target.name+'.history.jsonl').read_text().splitlines()), 3)

    def test_blocked_preflight_creates_no_files_or_requests(self):
        catalog = json.loads(DEFAULT.read_text())
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)/'not-created'
            for b in [binding(), binding(local_path=str(FIXTURE)), binding(local_path=str(FIXTURE), access_authorized=True)]:
                with patch('urllib.request.urlopen', side_effect=AssertionError('No external requests')):
                    self.assertEqual(collect(catalog, ['koges_cohort'], {'koges_cohort': b}, root)['status'], 'BLOCKED')
            self.assertFalse(root.exists())
            self.assertEqual(plan(catalog, ['finngen_gwas'], {'finngen_gwas': binding(access_authorized=True, authorization_ref='fixture')})['items'][0]['status'], 'BLOCKED')
            self.assertEqual(collect(catalog, ['ukb_ppp_pqtl'], {'ukb_ppp_pqtl': binding(sha256=None)}, root)['status'], 'BLOCKED')

    def test_interrupted_download_resumes_range_and_entity(self):
        requests = []
        offset = 7
        def opener(request, **options):
            requests.append(request)
            if len(requests) == 1:
                return Response(PAYLOAD[:offset], headers={'Content-Length': str(len(PAYLOAD)), 'ETag': '"fixture-v1"'})
            return Response(PAYLOAD[offset:], 206, {'Content-Range': f'bytes {offset}-{len(PAYLOAD)-1}/{len(PAYLOAD)}',
                                                  'Content-Length': str(len(PAYLOAD)-offset), 'ETag': '"fixture-v1"'})
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)/'fixture.tsv'
            with patch('masteromics.acquisition.time.sleep'):
                result = stage_file(target, SHA, url=URL, opener=opener, expected_bytes=len(PAYLOAD))
            self.assertEqual(result['status'], 'SUCCESS_VERIFIED')
            self.assertEqual(requests[1].get_header('Range'), f'bytes={offset}-')
            self.assertEqual(requests[1].get_header('If-range'), '"fixture-v1"')
            self.assertEqual(target.read_bytes(), PAYLOAD)
            self.assertFalse(target.with_name(target.name+'.partial').exists())

    def test_cross_run_resume_ignored_range_restart_and_identity_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)/'fixture.tsv'
            with self.assertRaises(OSError):
                stage_file(target, SHA, url=URL, retries=0, opener=lambda *a, **k: Response(PAYLOAD[:7], headers={'Content-Length': str(len(PAYLOAD))}))
            self.assertEqual(target.with_name(target.name+'.partial').read_bytes(), PAYLOAD[:7])
            with self.assertRaisesRegex(ValueError, 'identity changed'):
                stage_file(target, SHA, url=URL, identity={'release': 'changed'}, opener=lambda *a, **k: self.fail('Must not request'))
            requests = []
            def opener(request, **options):
                requests.append(request)
                return Response(PAYLOAD)  # Server ignores Range and sends 200.
            stage_file(target, SHA, url=URL, opener=opener)
            self.assertEqual(requests[0].get_header('Range'), 'bytes=7-')
            self.assertEqual(target.read_bytes(), PAYLOAD)

    def test_bad_range_checksum_and_gzip_never_promoted(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
                stage_file(p/'wrong.tsv', SHA, url=URL, opener=lambda *a, **k: Response(b'wrong'), retries=0)
            self.assertFalse((p/'wrong.tsv').exists())
            self.assertFalse((p/'wrong.tsv.partial').exists())
            # A matching byte checksum is insufficient if an artifact claims gzip.
            with self.assertRaises(OSError):
                stage_file(p/'fake.gz', SHA, url=URL, opener=lambda *a, **k: Response(PAYLOAD), retries=0)
            self.assertFalse((p/'fake.gz').exists())
            target = p/'range.tsv'
            with self.assertRaises(OSError):
                stage_file(target, SHA, url=URL, retries=0, opener=lambda *a, **k: Response(PAYLOAD[:7], headers={'Content-Length': str(len(PAYLOAD))}))
            with self.assertRaisesRegex(ValueError, 'Content-Range'):
                stage_file(target, SHA, url=URL, retries=0, opener=lambda *a, **k: Response(PAYLOAD[7:], 206, {'Content-Range': f'bytes 8-{len(PAYLOAD)-1}/{len(PAYLOAD)}'}))
            self.assertFalse(target.exists())

    def test_genuine_gzip_and_empty_artifact(self):
        import gzip
        body = gzip.compress(PAYLOAD, mtime=0)
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)/'fixture.gz'
            stage_file(target, hashlib.sha256(body).hexdigest(), url=URL,
                       opener=lambda *a, **k: Response(body), expected_bytes=len(body))
            self.assertEqual(gzip.decompress(target.read_bytes()), PAYLOAD)
            empty = Path(tmp)/'empty.tsv'
            with self.assertRaisesRegex(ValueError, 'Empty artifact'):
                stage_file(empty, hashlib.sha256(b'').hexdigest(), url=URL, opener=lambda *a, **k: Response(b''))
            self.assertFalse(empty.exists())

    def test_lock_and_symlink_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp); target = p/'fixture.tsv'; held = p/'fixture.tsv.lock'
            held.write_text('another process')
            with self.assertRaises(FileExistsError):
                stage_file(target, SHA, local_path=FIXTURE)
            self.assertEqual(held.read_text(), 'another process')
            held.unlink(); target.symlink_to(FIXTURE)
            with self.assertRaisesRegex(ValueError, 'symlink'):
                stage_file(target, SHA, local_path=FIXTURE)
            self.assertEqual(FIXTURE.read_bytes(), PAYLOAD)

    def test_real_loopback_http_interruption_and_resume(self):
        from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
        import threading
        import urllib.request
        requests = []
        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass
            def do_GET(self):
                offset = self.headers.get('Range')
                requests.append(offset)
                self.send_response(206 if offset else 200)
                if offset:
                    self.send_header('Content-Range', f'bytes 7-{len(PAYLOAD)-1}/{len(PAYLOAD)}')
                self.send_header('Content-Length', str(len(PAYLOAD)-7 if offset else len(PAYLOAD)))
                self.send_header('ETag', '"fixture-v1"')
                self.end_headers()
                self.wfile.write(PAYLOAD[7:] if offset else PAYLOAD[:7])
        server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True); thread.start()
        def opener(request, **options):
            # Loopback transport fixture substitutes only the connection destination.
            req = urllib.request.Request(f'http://127.0.0.1:{server.server_port}/fixture', headers=dict(request.header_items()))
            response = urllib.request.urlopen(req, **options)
            response.geturl = lambda: URL
            return response
        try:
            with tempfile.TemporaryDirectory() as tmp, patch('masteromics.acquisition.time.sleep'):
                target = Path(tmp)/'fixture.tsv'
                stage_file(target, SHA, url=URL, opener=opener)
                self.assertEqual(target.read_bytes(), PAYLOAD)
                self.assertEqual(requests, [None, 'bytes=7-'])
        finally:
            server.shutdown(); server.server_close(); thread.join()

    def test_collect_cli_without_scientific_dependencies(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp); b = p/'bindings.json'
            b.write_text(json.dumps({'ukb_ppp_pqtl': binding(local_path=str(FIXTURE))}))
            command = [sys.executable, '-S', '-m', 'masteromics', 'resources', 'collect', '--ids', 'ukb_ppp_pqtl',
                       '--bindings', str(b), '--root', str(p/'collection')]
            result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(json.loads(result.stdout)['status'], 'SUCCESS_VERIFIED')


if __name__ == '__main__':
    unittest.main()
