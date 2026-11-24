#!/usr/bin/env node
/**
 * Stage59: parse SHA256-pinned, browser-observed public HCCDB JavaScript AST.
 *
 * Route names and component props are not HTTP APIs. Report request-config
 * candidates ONLY where URL, method and call-site structure coexist in AST.
 * No backend requests, expression matrices, clinical/sample data, or patient
 * joins are performed. Literal source identifiers are never exported.
 */
import * as acorn from 'acorn';
import * as walk from 'acorn-walk';
import http from 'node:http';
import { createHash } from 'node:crypto';
import { mkdir, writeFile } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

export const PINNED = Object.freeze({
  download: Object.freeze({
    path: '/static/js/download.b22d03a9.js',
    bytes: 11335,
    sha256: '973cddbabef72e012a6be3dd2c418f71c5e4ac27c851bd4d70ebd7a28b42c777'
  }),
  repository: Object.freeze({
    path: '/static/js/repository.509a9848.js',
    bytes: 108379,
    sha256: 'cf02b87d3424e0f1605f0d7b772948df51f6e9c52820445a2850e13a3b568893'
  })
});
const HOST = 'lifeome.net';
const PORT = 809;
const MAX_RESPONSE = 512 * 1024;
const IMPORTANT = new Set([
  '/conf/download', '/conf/download2', '/conf/download3', '/datasetId',
  '/download-table', '/bulkDatasets', '/bulkrna', '/system/chart/getGPointData'
]);
const FILE_NAMES = ['HCCDB25_mRNA_level3.txt', 'HCCDB25.sample.txt', 'HCCDB25.patient.txt'];

export function verifyPinnedSource(buffer, status, headers, pin) {
  if (!Buffer.isBuffer(buffer) || buffer.length !== pin.bytes || buffer.length > MAX_RESPONSE)
    throw new Error('SOURCE_BYTE_LENGTH_MISMATCH');
  const digest = createHash('sha256').update(buffer).digest('hex');
  if (digest !== pin.sha256) throw new Error('SOURCE_SHA256_MISMATCH');
  if (status !== 200 && status !== 206) throw new Error('SOURCE_HTTP_UNVERIFIED');
  const mime = String(headers['content-type'] || '').toLowerCase().split(';')[0].trim();
  if (!['application/javascript','text/javascript','application/x-javascript',
        'text/plain','application/octet-stream'].includes(mime))
    throw new Error('SOURCE_JS_MIME_UNVERIFIED');
  let complete = false;
  if (status === 206) {
    const m = /^bytes (\d+)-(\d+)\/(\d+)$/.exec(String(headers['content-range'] || '').trim());
    complete = !!m && Number(m[1]) === 0 && Number(m[2]) + 1 === buffer.length &&
      Number(m[3]) === buffer.length;
  } else {
    complete = Number(headers['content-length']) === buffer.length;
  }
  if (!complete) throw new Error('SOURCE_RANGE_COMPLETENESS_UNVERIFIED');
  return {sha256:digest,bytes:buffer.length,complete:true,http_status:status};
}

export async function fetchPinnedSource(mode) {
  const pin = PINNED[mode];
  if (!pin) throw new Error('UNLISTED_PUBLIC_CHUNK');
  return await new Promise((resolvePromise, reject) => {
    const request = http.get({
      hostname: HOST, port: PORT, path: pin.path, protocol: 'http:',
      timeout: 30000, maxHeaderSize: 16384,
      headers: {'Accept':'application/javascript,*/*',
        'User-Agent':'HCC-Stage59-source-pinned-parse/1.0',
        'Range':'bytes=0-' + String(MAX_RESPONSE-1)}
    }, response => {
      if (response.headers.location) {
        response.resume(); reject(new Error('SOURCE_REDIRECT_REJECTED')); return;
      }
      const parts = [];
      let bytes = 0;
      response.on('data', chunk => {
        bytes += chunk.length;
        if (bytes > MAX_RESPONSE) {
          request.destroy(new Error('SOURCE_TOO_LARGE'));
          return;
        }
        parts.push(chunk);
      });
      response.on('end', () => {
        try {
          const full = Buffer.concat(parts, bytes);
          const source = verifyPinnedSource(full,response.statusCode,response.headers,pin);
          resolvePromise({buffer:full,source});
        } catch(error) {reject(error);}
      });
      response.on('error',reject);
    });
    request.on('timeout',()=>request.destroy(new Error('SOURCE_REQUEST_TIMEOUT')));
    request.on('error',reject);
  });
}

export function safeSourcePath(value) {
  if (typeof value !== 'string' || value.length > 150) return null;
  if (!value.startsWith('/') || value.includes('//') || value.includes('..') ||
      /[?#@]/.test(value) || /(?:OEP|OED|OES|OEZ)\d{4,}/i.test(value) ||
      /(?:patient|secret|token|password|private|auth)/i.test(value)) return null;
  if (!/^\/[A-Za-z0-9_./:$-]+$/.test(value)) return null;
  if (!/\/(?:api|conf|download|repository|bulk|data|dataset|system)\b/i.test(value)) return null;
  return value;
}

function fieldName(node,computed=false) {
  if (!node) return null;
  if (node.type==='Identifier' && !computed) return node.name;
  if (node.type==='Literal' && typeof node.value==='string') return node.value;
  return null;
}
function literalString(node) {
  return node?.type==='Literal' && typeof node.value==='string' ? node.value : null;
}
function memberMethod(callee) {
  if (!callee) return 'unknown';
  if (callee.type==='Identifier') {
    const name = callee.name.toLowerCase();
    return ['fetch','open','get','post','request','download'].includes(name) ? name : 'indirect';
  }
  if (callee.type==='MemberExpression') {
    const name=fieldName(callee.property,callee.computed);
    if (name && ['get','post','put','delete','patch','request','fetch','open'].includes(name.toLowerCase()))
      return name.toLowerCase();
    return 'indirect_member';
  }
  return 'indirect_expression';
}
function objectProperties(node) {
  if (!node || node.type !== 'ObjectExpression') return null;
  const found = {};
  for (const prop of node.properties) {
    if (prop.type !== 'Property' || prop.computed || prop.kind !== 'init') continue;
    const key = fieldName(prop.key);
    if (key && !Object.hasOwn(found,key)) found[key] = prop.value;
  }
  return found;
}
export function analyzeAst(input, mode='download') {
  const text = Buffer.isBuffer(input) ? input.toString('utf8') : input;
  if (typeof text !== 'string' || text.length > MAX_RESPONSE) throw new Error('INVALID_SOURCE_AST_BYTES');
  const root = acorn.parse(text,{sourceType:'script',ecmaVersion:'latest',allowReturnOutsideFunction:true});
  let totalStrings=0,matchedStrings=0,propertyKeys=0,directArgs=0,objectProps=0;
  const structuralClasses={};
  const watched=Object.fromEntries([...IMPORTANT].map(x=>[x,{
    ast_literal_occurrences:0,property_key:0,object_property_value:0,
    direct_call_argument:0,other_literal_context:0,
    verified_http_request_config:0
  }]));
  const requestShapes=[];
  const methods={};
  let callN=0,objectN=0;
  walk.fullAncestor(root,(node,ancestors)=>{
    if (node.type==='CallExpression') {
      callN++;
      const method=memberMethod(node.callee);
      methods[method]=(methods[method]||0)+1;
      for (const a of node.arguments) {
        const props=objectProperties(a);
        if (!props) continue;
        const url=literalString(props.url) || literalString(props.path);
        const verb=literalString(props.method);
        const safe=safeSourcePath(url);
        if (!safe) continue;
        const reqMethod=verb && /^(GET|POST|PUT|DELETE|PATCH|HEAD)$/i.test(verb)
          ? verb.toUpperCase() : null;
        // Transport method + explicit URL must occur in one invoked config.
        const eligible=!!reqMethod && !!props.url;
        if (eligible) {
          requestShapes.push({path:safe,method:reqMethod,
            structural_type:'SOURCE_DECLARED_INVOKED_REQUEST_CONFIG',
            invoked_via:method});
          if(watched[safe]) watched[safe].verified_http_request_config++;
        }
      }
      const direct=node.arguments.map(literalString).filter(Boolean).map(safeSourcePath).filter(Boolean);
      if (direct.length && ['get','post','fetch','request','open'].includes(method)) {
        for(const p of direct) requestShapes.push({path:p,method:method.toUpperCase(),
          structural_type:'SOURCE_DECLARED_DIRECT_CALL_ARGUMENT',invoked_via:method});
      }
    }
    if (node.type==='ObjectExpression') objectN++;
    if (node.type!=='Literal' || typeof node.value!=='string') return;
    totalStrings++;
    const safe=safeSourcePath(node.value);
    if(!safe) return;
    matchedStrings++;
    const parent=ancestors[ancestors.length-2];
    let kind='other_literal_context';
    if(parent?.type==='Property' && parent.key===node){kind='property_key';propertyKeys++;}
    else if(parent?.type==='Property' && parent.value===node){kind='object_property_value';objectProps++;}
    else if(parent?.type==='CallExpression' && parent.arguments.includes(node)){
      kind='direct_call_argument';directArgs++;
    }
    structuralClasses[kind]=(structuralClasses[kind]||0)+1;
    if(watched[safe]){
      watched[safe].ast_literal_occurrences++;
      watched[safe][kind]++;
    }
  });
  const shapes=[...new Map(requestShapes.map(x=>[JSON.stringify(x),x])).values()];
  return {
    module:mode,
    parse_status:'AST_PARSED',
    AST_node_source:'sha256_pinned_public_route_javascript',
    AST_literal_strings_n:totalStrings,
    safe_source_route_literal_occurrences_n:matchedStrings,
    AST_objects_n:objectN,
    AST_calls_n:callN,
    safe_literal_structural_classes:structuralClasses,
    observed_key_path_roles:watched,
    source_invocation_method_counts:methods,
    source_declared_request_shapes:shapes.slice(0,30),
    source_declared_request_shapes_n:shapes.length,
    exact_RNA_triplet_file_string_mentions:Object.fromEntries(
      FILE_NAMES.map(n=>[n,text.split(n).length-1])),
    source_declared_paths_are_not_verified_network_endpoints:true,
    actual_RNA_file_content_downloaded:false,
    official_PDC_RNA_identity_verified:false
  };
}
export function review(out) {
  if(out.original_rna_file_downloaded !== false || out.pdc_identity_matched !== false)
    throw new Error('UNJUSTIFIED_MULTIOMICS_GO');
  return {schema:'hcc-stage59-route-AST-review-v1',
    verified_source_ast:out.status==='VERIFIED_PINNED_AST_ANALYZED',
    verified_RNA_triplet:false,
    G30_B:'HOLD_OFFICIAL_PDC_QUANT_ALIQUOT_MAP',
    G30_C:'HOLD_NODE_HCCDB25_ORIGINAL_RNA_SOURCE',
    G30_D:'STOP_NO_SAME_BIOPSY_RNA_PROTEIN_JOIN',
    biological_claims_allowed:false};
}
async function main() {
  const i=process.argv.indexOf('--module'),j=process.argv.indexOf('--out');
  if(i<0 || !PINNED[process.argv[i+1]] || j<0 || !process.argv[j+1])
    throw new Error('INVALID_STAGE59_ARGS');
  const mode=process.argv[i+1],outPath=resolve(process.argv[j+1]);
  const out={schema:'hcc-stage59-public-AST-download-callsite-v1',module:mode,
    checked_utc:new Date().toISOString(),status:'SOURCE_AST_HOLD',
    original_rna_file_downloaded:false,pdc_identity_matched:false};
  try {
    const fetched=await fetchPinnedSource(mode);
    out.verified_source=fetched.source;
    out.analysis=analyzeAst(fetched.buffer,mode);
    out.status='VERIFIED_PINNED_AST_ANALYZED';
  } catch(error) {
    out.status='SOURCE_OR_AST_PARSE_HOLD';
    out.error_class=error?.constructor?.name || 'SourceError';
    out.error_code=String(error?.message||'SOURCE_HOLD').slice(0,65).replace(/[^A-Z0-9_]/g,'_');
  }
  out.review=review(out);
  await mkdir(dirname(outPath),{recursive:true});
  await writeFile(outPath,JSON.stringify(out,null,2)+'\n');
  console.log('HCC_STAGE59_AST',JSON.stringify(out));
}
if (process.argv[1] && import.meta.url===pathToFileURL(resolve(process.argv[1])).href) {
  main().catch(e=>{console.error('HCC_STAGE59_FATAL',e?.constructor?.name || 'Error');process.exitCode=1;});
}
