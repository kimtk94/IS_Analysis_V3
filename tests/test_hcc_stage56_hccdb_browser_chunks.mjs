import test from 'node:test';
import assert from 'node:assert/strict';
import { isAllowedPublicResource, sourceSafePath, sourceLiteralPath,
  auditJavaScript, postflight } from '../scripts/hcc/audit_stage56_hccdb_browser_chunks.mjs';

test('source allows only official static JS and exact public origin', () => {
  const good = 'http://lifeome.net:809/static/js/chunk-123456.abc.js';
  assert.equal(isAllowedPublicResource(good, 'script'), true);
  assert.equal(sourceSafePath(good), '/static/js/chunk-123456.abc.js');
  assert.equal(sourceSafePath('https://evil.invalid/static/js/chunk-123.js'), null);
  assert.equal(sourceSafePath('http://lifeome.net:809/static/js/app.js?token=x'), null);
  assert.equal(sourceSafePath('http://lifeome.net:809/download/HCCDB25.patient.txt'), null);
  assert.equal(sourceSafePath('http://lifeome.net:809/static/js/../secret.js'), null);
});
test('other resource categories are blocked before any network access', () => {
  const root = 'http://lifeome.net:809/';
  assert.equal(isAllowedPublicResource(root, 'document'), true);
  assert.equal(isAllowedPublicResource(root, 'fetch'), false);
  assert.equal(isAllowedPublicResource(root, 'xhr'), false);
  assert.equal(isAllowedPublicResource('http://lifeome.net:809/api/private', 'document'), false);
  assert.equal(isAllowedPublicResource('http://lifeome.net:809/static/css/app.css', 'stylesheet'), true);
  assert.equal(isAllowedPublicResource('http://lifeome.net:809/data', 'script'), false);
});
test('source-derived literals exclude private IDs, query strings and tokens', () => {
  assert.equal(sourceLiteralPath('/api/repository/download'), '/api/repository/download');
  assert.equal(sourceLiteralPath('/api/download/OED000123'), null);
  assert.equal(sourceLiteralPath('/api/private/file'), null);
  assert.equal(sourceLiteralPath('/api/file?secret=abc'), null);
  assert.equal(sourceLiteralPath('/api/download/123456'), null);
});
test('triplet mentions in public JS are not downloaded data', () => {
  const js = Buffer.from('const x="/api/repository/download"; ' +
    'const y="HCCDB25_mRNA_level3.txt"; const z="HCCDB25.sample.txt";' +
    'const a="HCCDB25.patient.txt";');
  const a = auditJavaScript(js, '/static/js/public-chunk.js');
  assert.equal(a.safe_literal_endpoint_candidates_n, 1);
  assert.equal(a.triplet_filename_occurrences['HCCDB25.patient.txt'], 1);
  assert.equal(a.source_is_official_file_download, false);
  assert.equal(a.source_is_patient_to_PDC_mapping, false);
  assert.equal(a.sha256.length, 64);
});
test('source-only hints never open clinical and biology GO', () => {
  const base = { rna_matrix_downloaded: false,
    official_biopsy_crosswalk_verified: false,
    route_checks: [{ route: 'download',
      status: 'PUBLIC_BROWSER_ROUTE_JS_OBSERVED',
      new_chunks_vs_home_n: 1, blocked_backend_api_requests_n: 12 }] };
  const verdict = postflight(base);
  assert.match(verdict.G30_C, /^HOLD/);
  assert.match(verdict.G30_D, /^STOP/);
  assert.equal(verdict.biological_claims_allowed, false);
  assert.throws(() => postflight({ ...base, rna_matrix_downloaded: true }));
});
