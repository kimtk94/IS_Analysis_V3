#!/usr/bin/env node
/**
 * HCC Stage60 — resolve the transport baseURL behind SHA256-pinned public SPA JS.
 *
 * A public URL string in a bundle is NOT a verified API host. Distinguish
 * invoked HTTP-client construction, default assignments, unused config literals
 * and library defaults. Never request metadata, binary files or patient data.
 */
import * as acorn from 'acorn';
import * as walk from 'acorn-walk';
import http from 'node:http';
import { createHash } from 'node:crypto';
import { mkdir, writeFile } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

export const SOURCES = Object.freeze({
  runtime:{path:'/static/js/app.3f5b30f0.js',size:8621,
    sha:'9c791fca32897684014bfbc856acb9c8aa8613aa0d9a5e6953c324a8f21fa781'},
  vendor:{path:'/static/js/chunk-libs.be0532e0.js',size:2074822,
    sha:'da6dd59c1414f7a8143433801837dcb120e32dfd1a48f37f5379b2c4499c7707'},
  helpers:{path:'/rna.js',size:97383,
    sha:'c77fb5bdd3266b5d6d7e926ba151b931280cec1effdeefb53b1b6670d65d2277'},
  download:{path:'/static/js/download.b22d03a9.js',size:11335,
    sha:'973cddbabef72e012a6be3dd2c418f71c5e4ac27c851bd4d70ebd7a28b42c777'},
  repository:{path:'/static/js/repository.509a9848.js',size:108379,
    sha:'cf02b87d3424e0f1605f0d7b772948df51f6e9c52820445a2850e13a3b568893'}
});
const SITE='http://lifeome.net:809';
const MAX=3*1024*1024;
const CONF=new Set(['/conf/download','/conf/download2','/conf/download3']);
export function sourceUrl(mode) {
  if(!Object.hasOwn(SOURCES,mode))throw new Error('UNOBSERVED_SOURCE_PATH');
  return SITE+SOURCES[mode].path;
}
export function validateSource(bytes,pin,status,headers) {
  if(!Buffer.isBuffer(bytes)||bytes.length!==pin.size||bytes.length>MAX)
    throw new Error('SOURCE_BYTE_SIZE_MISMATCH');
  const hash=createHash('sha256').update(bytes).digest('hex');
  if(hash!==pin.sha)throw new Error('SOURCE_SHA256_MISMATCH');
  const mime=String(headers['content-type']||'').split(';')[0].trim().toLowerCase();
  if(!['application/javascript','text/javascript','application/x-javascript',
       'text/plain','application/octet-stream'].includes(mime))
    throw new Error('PUBLIC_JS_MIME_MISMATCH');
  let complete=false;
  if(status===206) {
    const m=/^bytes (\d+)-(\d+)\/(\d+)$/.exec(String(headers['content-range']||'').trim());
    complete=!!m && Number(m[1])===0 && Number(m[2])+1===bytes.length &&
      Number(m[3])===bytes.length;
  } else if(status===200) {
    complete=Number(headers['content-length'])===bytes.length;
  }
  if(!complete)throw new Error('SOURCE_NOT_COMPLETELY_VERIFIED');
  return {source_bytes:bytes.length,sha256:hash,complete:true,http_status:status};
}
export async function fetchPublicMode(mode) {
  const pin=SOURCES[mode];
  if(!pin)throw new Error('UNOBSERVED_PUBLIC_JS');
  return await new Promise((ok,fail)=>{
    const req=http.get({hostname:'lifeome.net',port:809,path:pin.path,timeout:48000,
      headers:{'Range':'bytes=0-'+String(MAX-1),
        'Accept':'application/javascript,*/*',
        'User-Agent':'HCC-Stage60-public-request-config/1.0'}},r=>{
      if(r.headers.location){r.resume();fail(new Error('SOURCE_REDIRECT_REJECTED'));return;}
      const chunks=[];let count=0;
      r.on('data',data=>{count+=data.length;
        if(count>MAX){req.destroy(new Error('PUBLIC_JS_EXCEEDS_CAP'));return;}
        chunks.push(data);
      });
      r.on('end',()=>{
        try{
          const bytes=Buffer.concat(chunks,count);
          const source=validateSource(bytes,pin,r.statusCode,r.headers);
          ok({bytes,source});
        }catch(e){fail(e);}
      });
      r.on('error',fail);
    });
    req.on('timeout',()=>req.destroy(new Error('SOURCE_FETCH_TIMEOUT')));
    req.on('error',fail);
  });
}
function fieldName(node,computed=false) {
  if(!node)return null;
  if(node.type==='Identifier'&&!computed)return node.name;
  if(node.type==='Literal'&&typeof node.value==='string')return node.value;
  return null;
}
function literal(node) {
  return node?.type==='Literal'&&typeof node.value==='string' ? node.value : null;
}
function properties(node) {
  if(node?.type!=='ObjectExpression')return null;
  const result={};
  for(const p of node.properties){
    if(p.type!=='Property'||p.computed||p.kind!=='init')continue;
    const k=fieldName(p.key);
    if(k&&!Object.hasOwn(result,k))result[k]=p.value;
  }
  return result;
}
function methodName(callee) {
  if(callee?.type==='MemberExpression')return fieldName(callee.property,callee.computed)||'indirect';
  if(callee?.type==='Identifier')return callee.name;
  return 'indirect';
}
function isSensitive(s) {
  return /[?#@]|(?:OED|OEP|OES|OEZ)\d{4,}|(?:token|password|private|patient|bearer)/i.test(s);
}
export function classifyBaseURL(value) {
  if(typeof value!=='string'||value.length>200||isSensitive(value))
    return {class:'SENSITIVE_OR_NONLITERAL_HOLD'};
  if(value.startsWith('/')) {
    if(!/^\/[A-Za-z0-9_.\/-]*$/.test(value)||value.includes('..'))
      return {class:'NONCANONICAL_RELATIVE_HOLD'};
    return {class:'RELATIVE_SOURCE_LITERAL',path:value};
  }
  try {
    const u=new URL(value);
    if(!['http:','https:'].includes(u.protocol)||u.username||u.password||
       u.search||u.hash||u.pathname.includes('..'))
      return {class:'UNTRUSTED_URL_HOLD'};
    const port=u.port?':'+u.port:'';
    return {class:u.origin===SITE?'FIRST_PARTY_SOURCE_LITERAL':'OTHER_PUBLIC_ORIGIN_SOURCE_LITERAL',
      origin:u.protocol+'//'+u.hostname+port,path:u.pathname};
  }catch {return {class:'DYNAMIC_OR_INVALID_LITERAL_HOLD'};}
}
function classifyBaseNode(node) {
  const val=literal(node);
  return val===null ? {class:'DYNAMIC_EXPRESSION_UNRESOLVED',node_type:node?.type||'unknown'} :
    classifyBaseURL(val);
}
function targetPath(value) {
  if(typeof value!=='string'||isSensitive(value))return null;
  return CONF.has(value)?value:null;
}
export function analyzeTransportAst(source,mode) {
  const text=Buffer.isBuffer(source)?source.toString('utf8'):source;
  if(typeof text!=='string'||text.length>MAX)throw new Error('INVALID_PUBLIC_JS_SOURCE');
  const root=acorn.parse(text,{ecmaVersion:'latest',sourceType:'script',allowReturnOutsideFunction:true});
  const cfg=[];
  const assignments=[];
  const requests=[];
  const calleeMethods={};
  let objectCount=0,callCount=0;
  walk.full(root,node=>{
    if(node.type==='ObjectExpression')objectCount++;
    if(node.type==='CallExpression'){
      callCount++;
      const fn=methodName(node.callee);
      calleeMethods[fn]=(calleeMethods[fn]||0)+1;
      for(const [argIndex,arg] of node.arguments.entries()){
        const props=properties(arg);
        if(!props)continue;
        if(props.baseURL||props.baseUrl){
          const base=classifyBaseNode(props.baseURL||props.baseUrl);
          cfg.push({source:'INVOKED_OBJECT_CONFIG',method:fn,
            function_kind:node.callee.type,argument_index:argIndex,base});
        }
        const url=targetPath(literal(props.url));
        const method=literal(props.method);
        if(url && method && /^(GET|POST|HEAD)$/i.test(method)) {
          requests.push({source:'INVOKED_EXPLICIT_URL_METHOD',path:url,
            method:method.toUpperCase(),callee_kind:node.callee.type,
            transport_implementation_proved:false});
        }
      }
    }
    if(node.type==='AssignmentExpression'&&node.left?.type==='MemberExpression'){
      const name=fieldName(node.left.property,node.left.computed);
      if(name==='baseURL'||name==='baseUrl') {
        assignments.push({source:'BASEURL_MEMBER_ASSIGNMENT',
          left_object_kind:node.left.object.type,
          base:classifyBaseNode(node.right)});
      }
    }
    if(node.type==='Property'){
      const key=fieldName(node.key,node.computed);
      if(key==='baseURL'||key==='baseUrl') {
        cfg.push({source:'DECLARED_OBJECT_PROPERTY_NOT_PROVEN_INVOKED',
          key,base:classifyBaseNode(node.value)});
      }
    }
  });
  return {mode,parse_status:'SOURCE_AST_PARSED',AST_nodes:{objects:objectCount,calls:callCount},
    baseurl_evidence:cfg.slice(0,45),baseurl_evidence_count:cfg.length,
    baseurl_assignments:assignments.slice(0,25),baseurl_assignments_n:assignments.length,
    target_download_configs:requests.slice(0,20),target_download_configs_n:requests.length,
    call_method_counts:Object.fromEntries(Object.entries(calleeMethods)
      .filter(([k])=>/^(create|get|post|request|fetch|ajax)$/i.test(k))),
    source_only_not_live_http_origin:true,request_api_network_call_made:false,
    original_rna_triplet_obtained:false,pdc_sample_identity_verified:false};
}
export function review(out) {
  if(out.original_RNA_downloaded!==false||out.PDC_sample_join_verified!==false)
    throw new Error('BIOLOGICAL_RESULT_INVALID');
  return {schema:'hcc-stage60-transport-origin-evidence-gates-v1',
    actual_baseURL_and_public_metadata_API_verified:false,
    source_JS_verified:out.status==='SOURCE_AST_VERIFIED',
    original_RNA_matrix_verified:false,
    G30_B:'HOLD_PDC_OFFICIAL_ALIQUOT_MAP',
    G30_C:'HOLD_HCCDB25_ORIGINAL_RNA_SOURCE',
    G30_D:'STOP_NO_OFFICIAL_SAME_BIOPSY_JOIN',
    biological_claims_allowed:false};
}
async function main(){
  const im=process.argv.indexOf('--mode'),io=process.argv.indexOf('--out');
  if(im<0||io<0||!process.argv[io+1])throw new Error('INVALID_ARGUMENTS');
  const mode=process.argv[im+1];
  if(!Object.hasOwn(SOURCES,mode))throw new Error('MODE_NOT_SOURCE_VERIFIED');
  const out={schema:'hcc-stage60-source-transport-ast-v1',mode,
    checked_utc:new Date().toISOString(),status:'SOURCE_AST_HOLD',
    original_RNA_downloaded:false,PDC_sample_join_verified:false};
  try{
    const fetched=await fetchPublicMode(mode);
    out.source=fetched.source;
    out.analysis=analyzeTransportAst(fetched.bytes,mode);
    out.status='SOURCE_AST_VERIFIED';
  }catch(e){
    out.status='SOURCE_TRANSPORT_OR_AST_HOLD';
    out.error_class=e?.constructor?.name||'SourceError';
    out.error_code=String(e?.message||'UNKNOWN_ERROR').slice(0,60).replace(/[^A-Z0-9_]/g,'_');
  }
  out.review=review(out);
  const dest=resolve(process.argv[io+1]);
  await mkdir(dirname(dest),{recursive:true});
  await writeFile(dest,JSON.stringify(out,null,2)+'\n');
  console.log('HCC_STAGE60_TRANSPORT',JSON.stringify(out));
}
if(process.argv[1]&&import.meta.url===pathToFileURL(resolve(process.argv[1])).href)
  main().catch(e=>{console.error('HCC_STAGE60_FATAL',e?.constructor?.name||'Error');process.exitCode=1;});
