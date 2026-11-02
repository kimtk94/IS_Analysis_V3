#!/usr/bin/env node
/**
 * Stage56: observe HCCDB2 route-loaded JS with a clean, restricted browser.
 *
 * Only first-party public documents, styles and JS are permitted. All XHR,
 * fetch, fonts, images, websocket and third-party traffic are blocked.
 * Navigates the published SPA URLs; captures script content loaded by the
 * browser itself, never guesses Webpack chunk filenames or calls data APIs.
 * Reports public static-JS path/hash and anonymous structural aggregate only.
 */
import { createHash } from 'node:crypto';
import { mkdir, writeFile } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

const SITE = 'http://lifeome.net:809';
const EXPECTED = [
  'HCCDB25_mRNA_level3.txt',
  'HCCDB25.sample.txt',
  'HCCDB25.patient.txt'
];
const ROUTES = [
  { key: 'home', hash: '#/' },
  { key: 'repository_buke', hash: '#/repository?type=buke' },
  { key: 'download', hash: '#/download' }
];
const MAX_JS_BYTES = 4 * 1024 * 1024;
const MAX_SCRIPTS_PER_ROUTE = 35;
const MAX_TOTAL_JS_PER_ROUTE = 18 * 1024 * 1024;

export function isAllowedPublicResource(value, type) {
  try {
    const u = new URL(value);
    if (u.origin !== SITE || u.username || u.password) return false;
    if (u.search || (u.hash && type !== 'document')) return false;
    if (type === 'document') return u.pathname === '/' &&
      (!u.hash || ['#/', '#/repository?type=buke', '#/download'].includes(u.hash));
    if (type === 'stylesheet') return /^\/(?:static\/)?css\/[A-Za-z0-9_.\-/]+\.css$/.test(u.pathname);
    if (type === 'script') {
      const valid = /^\/(?:static\/)?js\/[A-Za-z0-9_.\-/]+\.js$/.test(u.pathname);
      const parts = decodeURIComponent(u.pathname).split('/');
      return valid && !parts.includes('.') && !parts.includes('..');
    }
    return false;
  } catch {
    return false;
  }
}

export function sourceSafePath(value) {
  if (!isAllowedPublicResource(value, 'script')) return null;
  return new URL(value).pathname;
}

export function sourceLiteralPath(value) {
  if (typeof value !== 'string' || value.length > 120 ||
      /[?#]/.test(value) || /(?:OED|OES|OEZ|OEP)\d{5,}/i.test(value) ||
      /(?:token|key|secret|private|patient)/i.test(value) ||
      /\d{5,}/.test(value)) return null;
  if (!/^\/(?:api|download|repository|data)\/[A-Za-z0-9_./:$-]+$/.test(value)) return null;
  if (value.includes('..')) return null;
  return value;
}

export function auditJavaScript(body, sourcePath) {
  if (!Buffer.isBuffer(body)) throw new TypeError('expected public static JS Buffer');
  const text = body.toString('utf8');
  const routePatterns = {
    repository: /\/repository/gi,
    download: /\/download/gi,
    buke: /buke/gi,
    HCCDB25: /HCCDB25/gi,
    dynamic_loader: /(?:__webpack_require__|webpackJsonp|webpackChunkName|\.e\(|\.u\(|import\()/g
  };
  const count = Object.fromEntries(Object.entries(routePatterns).map(([key, regex]) =>
    [key, [...text.matchAll(regex)].length]));
  const strings = [...text.matchAll(/["'](\/(?:api|download|repository|data)\/[A-Za-z0-9_./:$-]{2,120})["']/g)]
    .map(m => sourceLiteralPath(m[1])).filter(Boolean);
  const fileMentions = Object.fromEntries(EXPECTED.map(filename => [
    filename, text.split(filename).length - 1
  ]));
  return {
    public_js_path: sourcePath,
    inspected_bytes: body.length,
    sha256: createHash('sha256').update(body).digest('hex'),
    marker_counts: count,
    safe_literal_endpoint_paths: [...new Set(strings)].sort().slice(0, 30),
    safe_literal_endpoint_candidates_n: new Set(strings).size,
    triplet_filename_occurrences: fileMentions,
    source_is_official_file_download: false,
    source_is_patient_to_PDC_mapping: false
  };
}

export function postflight(result) {
  if (result.rna_matrix_downloaded !== false ||
      result.official_biopsy_crosswalk_verified !== false) {
    throw new Error('UNAUTHORIZED_BIOLOGICAL_GO');
  }
  return {
    schema: 'hcc-stage56-browser-observed-chunks-postflight-v1',
    route_load_validation: result.route_checks.map(r => ({
      route: r.route,
      status: r.status,
      new_public_static_chunk_n: r.new_chunks_vs_home_n ?? null,
      blocked_backend_api_requests_n: r.blocked_backend_api_requests_n
    })),
    original_HCCDB25_expression_sample_patient_triplet_verified: false,
    G30_B: 'HOLD_OFFICIAL_PDC_ALIQUOT_JOIN',
    G30_C: 'HOLD_ORIGINAL_HCCDB25_RNA_TRIPLET_AND_IDS',
    G30_D: 'STOP_NO_SOURCE_AUTHORIZED_PATIENT_MATCH',
    biological_claims_allowed: false
  };
}

export async function inspectRoute(browser, def) {
  const context = await browser.newContext({ storageState: { cookies: [], origins: [] },
    serviceWorkers: 'block', acceptDownloads: false });
  const page = await context.newPage();
  const tasks = [];
  const assets = new Map();
  const result = {
    route: def.key,
    status: 'BROWSER_ROUTE_SOURCE_HOLD',
    blocked_backend_api_requests_n: 0,
    blocked_third_party_or_nonessential_requests_n: 0,
    blocked_resource_type_counts: {},
    attempted_downloads_n: 0,
    scripts_observed_n: 0,
    script_source_errors_n: 0
  };
  await context.route('**/*', async route => {
    const req = route.request();
    const rt = req.resourceType();
    if (rt === 'xhr' || rt === 'fetch' || rt === 'websocket') {
      result.blocked_resource_type_counts[rt] = (result.blocked_resource_type_counts[rt] || 0) + 1;
      result.blocked_backend_api_requests_n += 1;
      return route.abort('blockedbyclient');
    }
    if (!isAllowedPublicResource(req.url(), rt)) {
      result.blocked_resource_type_counts[rt] = (result.blocked_resource_type_counts[rt] || 0) + 1;
      result.blocked_third_party_or_nonessential_requests_n += 1;
      return route.abort('blockedbyclient');
    }
    return route.continue();
  });
  page.on('download', async download => {
    result.attempted_downloads_n += 1;
    await download.cancel().catch(() => {});
  });
  page.on('response', response => {
    const url = response.url();
    const safe = sourceSafePath(url);
    if (!safe || response.request().resourceType() !== 'script') return;
    const task = (async () => {
      result.scripts_observed_n += 1;
      if (result.scripts_observed_n > MAX_SCRIPTS_PER_ROUTE ||
          response.status() !== 200) return;
      try {
        const headers = await response.allHeaders();
        const declared = Number(headers['content-length'] ?? '0');
        if (declared > MAX_JS_BYTES) {
          result.script_source_errors_n += 1;
          return;
        }
        const content = await response.body();
        if (content.length > MAX_JS_BYTES) {
          result.script_source_errors_n += 1;
          return;
        }
        assets.set(safe, auditJavaScript(content, safe));
      } catch {
        result.script_source_errors_n += 1;
      }
    })();
    tasks.push(task);
  });
  try {
    // DOMContentLoaded waits on every slow upstream script; commit captures the
    // official document response without upgrading a partial SPA render to GO.
    const nav = await page.goto(SITE + '/' + def.hash, { waitUntil: 'commit', timeout: 30000 });
    result.document_http_status = nav?.status() ?? null;
    await page.waitForTimeout(8500);
    result.document_ready_state = await page.evaluate(() => document.readyState).catch(() => 'unknown');
    await Promise.race([Promise.allSettled(tasks), page.waitForTimeout(6000)]);
    const sorted = [...assets.values()].sort((a, b) => a.public_js_path.localeCompare(b.public_js_path));
    const size = sorted.reduce((a, v) => a + v.inspected_bytes, 0);
    if (size > MAX_TOTAL_JS_PER_ROUTE) {
      result.status = 'SCRIPT_TOTAL_BYTES_OVER_LIMIT_HOLD';
    } else {
      result.status = sorted.length ? 'PUBLIC_BROWSER_ROUTE_JS_OBSERVED' : 'NO_PUBLIC_ROUTE_JS_OBSERVED';
    }
    result.assets = sorted;
    result.script_bytes_inspected_n = size;
  } catch (error) {
    result.error_class = error?.constructor?.name || 'BrowserNavigationError';
    result.assets = [...assets.values()];
  } finally {
    await context.close();
  }
  return result;
}

export async function runAudit(browser) {
  const routes = [];
  for (const def of ROUTES) routes.push(await inspectRoute(browser, def));
  const home = new Set(routes.find(r => r.route === 'home')?.assets?.map(x => x.public_js_path) || []);
  const output = {
    schema: 'hcc-stage56-hccdb25-browser-observed-static-routes-v1',
    checked_utc: new Date().toISOString(),
    official_website: SITE + '/',
    network_policy: 'first_party_static_documents_scripts_styles_only;block_xhr_fetch_websocket_downloads',
    rna_matrix_downloaded: false,
    official_biopsy_crosswalk_verified: false,
    route_checks: routes.map(r => ({
      ...r,
      new_chunks_vs_home: (r.assets || []).map(a => a.public_js_path).filter(p => !home.has(p)),
      new_chunks_vs_home_n: (r.assets || []).filter(a => !home.has(a.public_js_path)).length
    }))
  };
  output.review = postflight(output);
  return output;
}

async function main() {
  const args = process.argv;
  const ix = args.indexOf('--out');
  if (ix < 0 || !args[ix + 1]) throw new Error('Expected --out');
  const path = resolve(args[ix + 1]);
  let report;
  let browser;
  try {
    const { chromium } = await import('playwright');
    browser = await chromium.launch({ headless: true });
    report = await runAudit(browser);
  } catch (error) {
    report = {
      schema: 'hcc-stage56-hccdb25-browser-observed-static-routes-v1',
      checked_utc: new Date().toISOString(),
      status: 'BROWSER_LAUNCH_OR_SOURCE_HOLD',
      error_class: error?.constructor?.name || 'BrowserError',
      rna_matrix_downloaded: false,
      official_biopsy_crosswalk_verified: false,
      route_checks: [],
    };
    report.review = postflight(report);
  } finally {
    if (browser) await browser.close().catch(() => {});
  }
  await mkdir(dirname(path), { recursive: true });
  await writeFile(path, JSON.stringify(report, null, 2) + '\n');
  console.log('HCC_STAGE56_BROWSER', JSON.stringify(report));
}

if (process.argv[1] && import.meta.url === pathToFileURL(resolve(process.argv[1])).href) {
  main().catch(error => {
    console.error('HCC_STAGE56_FATAL', error?.constructor?.name || 'Error');
    process.exitCode = 1;
  });
}
