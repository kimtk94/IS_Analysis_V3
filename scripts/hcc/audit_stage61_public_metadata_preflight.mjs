#!/usr/bin/env node
/**
 * Stage61: source-backed public HEAD and restricted JSON metadata preflight.
 *
 * Positive scientific claims require ORIGINAL RNA expression bytes, authoritative
 * sample/patient links and a separate official PDC aliquot bridge. This probe
 * deliberately obtains none of those. No arbitrary URLs, credentials, redirects,
 * private APIs, downloadable binary content, expression values or patient IDs.
 *
 * All three paths are AST-declared literal GET request configs in TWO separately
 * SHA-pinned, browser-observed HCCDB25 route JS bundles. An actual runtime
 * binding to the same axios instance is NOT presumed from co-occurrence.
 */
import * as acorn from 'acorn';
import * as walk from 'acorn-walk';
import http from 'node:http';
import { createHash } from 'node:crypto';
import { mkdir, writeFile } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

export const SITE = 'http://lifeome.net:809';
export const BUNDLES = Object.freeze({
  repository: Object.freeze({path:'/static/js/repository.509a9848.js',
    size:108379,sha256:'cf02b87d3424e0f1605f0d7b772948df51f6e9c52820445a2850e13a3b568893'}),
  download: Object.freeze({path:'/static/js/download.b22d03a9.js',
    size:11335,sha256:'973cddbabef72e012a6be3dd2c418f71c5e4ac27c851bd4d70ebd7a28b42c777'})
});
export const ROUTES=Object.freeze({
  download: '/api/conf/download',
  download2: '/api/conf/download2',
  download3: '/api/conf/download3'
});
const TARGETS=['HCCDB25_mRNA_level3.txt','HCCDB25.sample.txt','HCCDB25.patient.txt'];
const MAX_JS=256*1024,MAX_JSON=64*1024;
const ORIGIN='lifeome.net',PORT=809;

export function normalizedEndpoint(key) {
  if(!Object.hasOwn(ROUTES,key))throw new Error('ENDPOINT_NOT_SOURCE_AUTHORIZED');
  return ROUTES[key];
}
export function isJsonMime(value) {
  const x=String(value||'').toLowerCase().split(';')[0].trim();
  return x==='application/json'||x==='application/problem+json'||x==='application/hal+json';
}
export function strictHeaders(status,headers,expectedLength,maxLength=MAX_JS) {
  if(![200,206].includes(status)||headers.location)
    throw new Error('UNVERIFIED_SOURCE_REDIRECT_OR_STATUS');
  const len=Number(headers['content-length']);
  if(status===200) {
    if(!Number.isInteger(len)||len!==expectedLength)throw new Error('SOURCE_CONTENT_LENGTH_MISMATCH');
  }else{
    const m=/^bytes (\d+)-(\d+)\/(\d+)$/.exec(String(headers['content-range']||'').trim());
    if(!m||Number(m[1])!==0||Number(m[2])+1!==expectedLength||Number(m[3])!==expectedLength)
      throw new Error('SOURCE_RANGE_NOT_COMPLETE');
  }
  if(expectedLength>maxLength)throw new Error('SOURCE_EXCEEDS_SIZE_LIMIT');
  const mime=String(headers['content-type']||'').toLowerCase().split(';')[0].trim();
  if(!['application/javascript','text/javascript','application/x-javascript',
       'text/plain','application/octet-stream'].includes(mime))
    throw new Error('UNVERIFIED_JS_MIME');
}
export function verifyBundle(data,pin,status,headers) {
  if(!Buffer.isBuffer(data)||data.length!==pin.size||data.length>MAX_JS)
    throw new Error('PUBLIC_SOURCE_BYTES_DIFFER');
  strictHeaders(status,headers,data.length);
  const actual=createHash('sha256').update(data).digest('hex');
  if(actual!==pin.sha256)throw new Error('PUBLIC_JS_SOURCE_HASH_DIFFER');
  return {sha256:actual,bytes:data.length,pinned_source_verified:true};
}
function sourceStaticUrl(key) {
  if(!Object.hasOwn(BUNDLES,key))throw new Error('UNOBSERVED_JAVASCRIPT_BUNDLE');
  return BUNDLES[key].path;
}
export async function sourceBytes(key,fetcher=http.request) {
  const path=sourceStaticUrl(key),pin=BUNDLES[key];
  return new Promise((ok,fail)=>{
    const req=fetcher({protocol:'http:',hostname:ORIGIN,port:PORT,path,method:'GET',
      timeout:25000,headers:{Accept:'application/javascript',
      Range:'bytes=0-'+String(MAX_JS-1),'User-Agent':'HCC-Stage61-public-metadata-audit/1.0'}},res=>{
      if(res.headers.location) {res.resume();fail(new Error('SOURCE_REDIRECT_REJECTED'));return;}
      const chunks=[];let total=0;
      res.on('data',part=>{
        total+=part.length;
        if(total>MAX_JS) {req.destroy(new Error('SOURCE_JS_OVERSIZE'));return;}
        chunks.push(part);
      });
      res.on('end',()=>{
        try{
          const data=Buffer.concat(chunks,total),meta=verifyBundle(data,pin,res.statusCode,res.headers);
          ok({data,meta});
        }catch(error){fail(error);}
      });
      res.on('error',fail);
    });
    req.on('timeout',()=>req.destroy(new Error('SOURCE_HTTP_TIMEOUT')));
    req.on('error',fail);req.end();
  });
}
function propName(key) {
  if(key?.type==='Identifier')return key.name;
  if(key?.type==='Literal'&&typeof key.value==='string')return key.value;
  return null;
}
function objProps(obj) {
  if(obj?.type!=='ObjectExpression')return null;
  const out={};
  for(const p of obj.properties){
    if(p.type!=='Property'||p.computed)continue;
    const k=propName(p.key);
    if(k&&!Object.hasOwn(out,k))out[k]=p.value;
  }
  return out;
}
function stringValue(node) {
  return node?.type==='Literal'&&typeof node.value==='string'?node.value:null;
}
function moduleContext(ancestors) {
  // An ancestor Webpack module FunctionExpression is evidence of co-location,
  // not proof of exact client binding. Return only its ordinal source offset.
  for(const a of ancestors) {
    if(a.type==='Property'&&(a.value?.type==='FunctionExpression'||a.value?.type==='ArrowFunctionExpression')) {
      return a.start;
    }
  }
  return null;
}
export function sourceCooccurrenceAnalysis(source,kind) {
  const ast=acorn.parse(Buffer.isBuffer(source)?source.toString('utf8'):source,
    {ecmaVersion:'latest',sourceType:'script',allowReturnOutsideFunction:true});
  const creates=[],calls=[];
  walk.fullAncestor(ast,(n,a)=>{
    if(n.type!=='CallExpression')return;
    const m=n.callee?.type==='MemberExpression'?propName(n.callee.property):null;
    for(const arg of n.arguments){
      const props=objProps(arg);
      if(!props)continue;
      const base=stringValue(props.baseURL);
      if(m==='create'&&base==='/api'){
        creates.push({context:moduleContext(a),verified_invoked_create_config:true});
      }
      const url=stringValue(props.url);
      const method=stringValue(props.method);
      if(method?.toUpperCase()==='GET' &&
         ['/conf/download','/conf/download2','/conf/download3'].includes(url)) {
        calls.push({path:'/api'+url,method:'GET',context:moduleContext(a),
          client_binding_proved:false});
      }
    }
  });
  const linkedCalls=calls.filter(c=>creates.some(x=>x.context!==null && x.context===c.context));
  const paths=[...new Set(calls.map(x=>x.path))].sort();
  const collocatedPaths=[...new Set(linkedCalls.map(x=>x.path))].sort();
  return {
    bundle:kind,
    invoked_create_baseurl_api_count:creates.length,
    exact_get_download_config_count:calls.length,
    declared_api_get_paths:paths,
    source_module_collocated_api_get_paths:collocatedPaths,
    same_source_module_collocation_n:collocatedPaths.length,
    exact_request_instance_binding_verified:false,
    actual_network_request_observed:false
  };
}
export function sourceGate(analyses) {
  if(analyses.length!==2)throw new Error('TWO_INDEPENDENT_BUNDLES_REQUIRED');
  for(const a of analyses){
    if(a.invoked_create_baseurl_api_count<1)throw new Error('MISSING_CALLED_BASEURL');
    if(a.exact_get_download_config_count<3)throw new Error('MISSING_THREE_GET_CONFIGS');
    if(!Object.values(ROUTES).every(x=>a.declared_api_get_paths.includes(x)))
      throw new Error('MISSING_SOURCE_DECLARED_ENDPOINT');
  }
  return {independent_bundles_verified_n:2,
    source_declared_bounded_public_head_permitted:true,
    actual_http_transport_instance_binding_verified:false,
    source_proves_api_server_returns_metadata:false};
}
export function classifyHead(status,headers) {
  const type=String(headers['content-type']||'').toLowerCase().split(';')[0].trim();
  const raw=headers['content-length'],length=raw!==undefined&&raw!==null?Number(raw):null;
  const json=isJsonMime(type);
  const allowBody= status===200&&!headers.location&&json&&Number.isInteger(length)&&
    length>=2&&length<=MAX_JSON;
  return {
    status_code:Number.isInteger(status)?status:null,
    response_mime_class:json?'JSON':type.startsWith('text/html')?'HTML':
      type.includes('octet-stream')?'BINARY':'OTHER_OR_UNKNOWN',
    has_redirect:Boolean(headers.location),
    reported_length_bytes:Number.isInteger(length)?length:null,
    metadata_get_allowlisted:allowBody,
    reason:allowBody?'BOUNDED_PUBLIC_JSON_HEAD_VERIFIED':
       status===405?'HEAD_METHOD_UNSUPPORTED_HOLD':
       status===401||status===403?'PUBLIC_ACCESS_RESTRICTED_HOLD':
       status===404?'PUBLIC_ENDPOINT_NOT_FOUND_HOLD':
       headers.location?'PUBLIC_REDIRECT_REJECTED':
       'RESPONSE_NOT_SAFE_FOR_JSON_FETCH'
  };
}
function httpOnce(path,method,cap=0,fetcher=http.request) {
  if(!Object.values(ROUTES).includes(path)||!['HEAD','GET'].includes(method))
    throw new Error('UNLISTED_HTTP_SOURCE_ENDPOINT');
  return new Promise((ok,fail)=>{
    const req=fetcher({protocol:'http:',hostname:ORIGIN,port:PORT,path,method,
      timeout:16000,headers:{Accept:'application/json',
        'User-Agent':'HCC-Stage61-unauthed-public-metadata/1.0'}},res=>{
      if(method==='HEAD') {
        // Headers only; never consume a response body.
        const info={status:res.statusCode,headers:res.headers};
        res.destroy();ok(info);return;
      }
      if(res.statusCode!==200||res.headers.location||!isJsonMime(res.headers['content-type'])){
        res.destroy();fail(new Error('PUBLIC_METADATA_JSON_MIME_OR_STATUS_CHANGED'));return;
      }
      const length=Number(res.headers['content-length']);
      if(!Number.isInteger(length)||length<2||length>cap) {
        res.destroy();fail(new Error('PUBLIC_METADATA_BODY_LENGTH_UNSAFE'));return;
      }
      const chunks=[];let count=0;
      res.on('data',data=>{
        count+=data.length;
        if(count>cap){req.destroy(new Error('PUBLIC_JSON_TOO_LARGE'));return;}
        chunks.push(data);
      });
      res.on('end',()=>{
        if(count!==length){fail(new Error('PUBLIC_JSON_LENGTH_MISMATCH'));return;}
        ok({status:res.statusCode,headers:res.headers,body:Buffer.concat(chunks,count)});
      });
      res.on('error',fail);
    });
    req.on('timeout',()=>req.destroy(new Error('PUBLIC_API_HEADER_TIMEOUT')));
    req.on('error',fail);req.end();
  });
}
export function anonymousJsonSummary(data) {
  // No raw URLs, samples, IDs, expression values, or non-target filenames.
  if(data.length>MAX_JSON)throw new Error('PUBLIC_METADATA_TOO_LARGE');
  const v=JSON.parse(data.toString('utf8'));
  const root=Array.isArray(v)?'array':v===null?'null':typeof v;
  let topKeys=[],topArrayLength=0;
  if(Array.isArray(v))topArrayLength=v.length;
  if(root==='object'){
    // Only known non-PHI field names may be reported.
    const keys=new Set(['code','status','message','data','result','rows','files',
                        'total','page','list','success','count','items']);
    topKeys=Object.keys(v).filter(x=>keys.has(x)).sort();
  }
  const raw=data.toString('utf8');
  return {root_type:root,top_level_whitelisted_metadata_keys:topKeys,
    top_level_array_n:topArrayLength,
    expected_triplet_exact_filename_mentions:Object.fromEntries(
      TARGETS.map(name=>[name,raw.split(name).length-1])),
    observed_response_bytes:data.length,
    response_sha256:createHash('sha256').update(data).digest('hex'),
    raw_metadata_or_source_url_exported:false,
    RNA_matrix_or_patient_file_downloaded:false};
}
export function review(out) {
  if(out.rna_original_downloaded!==false||out.pdc_sample_join_verified!==false)
    throw new Error('BIOLOGICAL_GO_PROHIBITED');
  return {schema:'hcc-stage61-public-api-metadata-gate-v1',
    public_metadata_index_verified:false,
    original_RNA_triplet_verified:false,
    G30_B:'HOLD_OFFICIAL_PDC_ALIQUOT_MAPPING',
    G30_C:'HOLD_SOURCE_CONFIRMED_HCCDB25_RNA_FILES',
    G30_D:'STOP_NO_PDC_RNA_SAME_BIOPSY_BRIDGE',
    biological_analysis_permitted:false};
}
export async function inspect(endpoint,load=sourceBytes,request=httpOnce) {
  const out={schema:'hcc-stage61-restricted-public-meta-preflight-v1',
    checked_utc:new Date().toISOString(),endpoint_key:endpoint,
    status:'SOURCE_IDENTITY_HOLD',rna_original_downloaded:false,
    pdc_sample_join_verified:false};
  try{
    const bundles=[];
    for(const type of ['repository','download']){
      const fetched=await load(type);
      bundles.push({source:type,meta:fetched.meta,
        analysis:sourceCooccurrenceAnalysis(fetched.data,type)});
    }
    out.sha256_pinned_source_bundles=bundles.map(x=>({source:x.source,...x.meta}));
    const gate=sourceGate(bundles.map(x=>x.analysis));
    out.source_provenance=gate;
    out.module_collocation_counts=bundles.map(x=>({
      source:x.source,collocated:x.analysis.same_source_module_collocation_n}));
    out.status='SOURCE_DOUBLE_CONFIRMED';
    const head=await request(normalizedEndpoint(endpoint),'HEAD');
    const classification=classifyHead(head.status,head.headers);
    out.public_head=classification;
    out.status='PUBLIC_HEAD_PREFLIGHT_COMPLETE';
    if(classification.metadata_get_allowlisted) {
      try{
        const response=await request(normalizedEndpoint(endpoint),'GET',MAX_JSON);
        out.public_json_metadata=anonymousJsonSummary(response.body);
        out.status='PUBLIC_BOUNDED_JSON_METADATA_READ';
      }catch(error){
        out.metadata_hold_error_class=error?.constructor?.name||'MetadataError';
      }
    }
  }catch(error){
    out.error_class=error?.constructor?.name||'PreflightError';
    out.error_code=String(error?.message||'SOURCE_HOLD').slice(0,70).replace(/[^A-Z0-9_]/g,'_');
  }
  out.review=review(out);
  return out;
}
async function main(){
  const ia=process.argv.indexOf('--endpoint'),io=process.argv.indexOf('--out');
  if(ia<0||io<0||!process.argv[io+1])throw new Error('INVALID_ARGUMENTS');
  const name=process.argv[ia+1],dest=resolve(process.argv[io+1]);
  if(!Object.hasOwn(ROUTES,name))throw new Error('UNLISTED_PUBLIC_PREFLIGHT');
  const out=await inspect(name);
  await mkdir(dirname(dest),{recursive:true});
  await writeFile(dest,JSON.stringify(out,null,2)+'\n');
  console.log('HCC_STAGE61_PREFLIGHT',JSON.stringify(out));
}
if(process.argv[1]&&import.meta.url===pathToFileURL(resolve(process.argv[1])).href)
  main().catch(error=>{console.error('HCC_STAGE61_FATAL',error?.constructor?.name||'Error');process.exitCode=1;});
