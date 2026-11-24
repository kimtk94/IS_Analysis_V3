import test from 'node:test';
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {SOURCES,sourceUrl,validateSource,classifyBaseURL,analyzeTransportAst,review}
 from '../scripts/hcc/audit_stage60_transport_origin_ast.mjs';

test('all sources are SHA-pinned to prior verified official public JS',()=>{
 assert.equal(SOURCES.runtime.size,8621);
 assert.equal(SOURCES.vendor.size,2074822);
 assert.equal(SOURCES.download.size,11335);
 assert.equal(SOURCES.repository.size,108379);
 assert.ok(Object.values(SOURCES).every(x=>/^[a-f0-9]{64}$/.test(x.sha)));
 assert.equal(sourceUrl('runtime'),'http://lifeome.net:809/static/js/app.3f5b30f0.js');
 assert.throws(()=>sourceUrl('api-guess'));
});
test('exact source SHA, bytes, MIME and full HTTP range required',()=>{
 const bytes=Buffer.from('abc');
 const pin={size:3,sha:createHash('sha256').update(bytes).digest('hex')};
 const h={'content-type':'application/javascript','content-range':'bytes 0-2/3'};
 assert.equal(validateSource(bytes,pin,206,h).complete,true);
 assert.throws(()=>validateSource(bytes,{...pin,sha:'0'.repeat(64)},206,h));
 assert.throws(()=>validateSource(bytes,pin,206,{...h,'content-range':'bytes 0-2/300'}));
 assert.throws(()=>validateSource(bytes,pin,206,{...h,'content-type':'text/html'}));
});
test('source URL strings without secure provenance are NOT trusted',()=>{
 assert.deepEqual(classifyBaseURL('http://lifeome.net:809/'),{
  class:'FIRST_PARTY_SOURCE_LITERAL',origin:'http://lifeome.net:809',path:'/'});
 assert.deepEqual(classifyBaseURL('/api'),{class:'RELATIVE_SOURCE_LITERAL',path:'/api'});
 assert.equal(classifyBaseURL('https://some-public-host.org/v1').class,
  'OTHER_PUBLIC_ORIGIN_SOURCE_LITERAL');
 assert.equal(classifyBaseURL('https://secret:test@other.org/api').class,
  'SENSITIVE_OR_NONLITERAL_HOLD');
 assert.equal(classifyBaseURL('/api?secret=1').class,
  'SENSITIVE_OR_NONLITERAL_HOLD');
});
test('baseURL in unused config object must not be promoted to invoked source',()=>{
 const data='var conf={baseURL:"https://api.example.org/v1"};';
 const result=analyzeTransportAst(data,'runtime');
 assert.equal(result.baseurl_evidence_count,1);
 assert.equal(result.baseurl_evidence[0].source,
  'DECLARED_OBJECT_PROPERTY_NOT_PROVEN_INVOKED');
 assert.equal(result.request_api_network_call_made,false);
});
test('actual create function with explicit baseURL is an invoked config only',()=>{
 const data='var x=client.create({baseURL:"https://api.example.org/v1"});';
 const result=analyzeTransportAst(data,'runtime');
 assert.ok(result.baseurl_evidence.some(x=>x.source==='INVOKED_OBJECT_CONFIG'
   &&x.method==='create'&&x.base.origin==='https://api.example.org'));
 assert.equal(result.source_only_not_live_http_origin,true);
});
test('dynamic runtime baseURL is HOLD, never resolved by guess',()=>{
 const result=analyzeTransportAst('let c=f({baseURL:window.env.host});','helpers');
 assert.ok(result.baseurl_evidence.some(x=>x.base.class==='DYNAMIC_EXPRESSION_UNRESOLVED'));
});
test('baseURL assignment receives distinct evidence class',()=>{
 const result=analyzeTransportAst('client.defaults.baseURL="http://lifeome.net:809";','runtime');
 assert.equal(result.baseurl_assignments_n,1);
 assert.equal(result.baseurl_assignments[0].base.class,'FIRST_PARTY_SOURCE_LITERAL');
});
test('request-config URL/method is source declaration, not network validation',()=>{
 const result=analyzeTransportAst('const z=req({url:"/conf/download",method:"GET"});','download');
 assert.equal(result.target_download_configs_n,1);
 assert.equal(result.target_download_configs[0].method,'GET');
 assert.equal(result.target_download_configs[0].transport_implementation_proved,false);
});
test('simple unrelated front-end strings are not request configs',()=>{
 const result=analyzeTransportAst('const x={path:"/conf/download",method:"GET"};','repository');
 assert.equal(result.target_download_configs_n,0);
});
test('private original IDs, user records and query strings not exported',()=>{
 const result=analyzeTransportAst('x({url:"/conf/download?token=abc",method:"GET"});','runtime');
 assert.equal(JSON.stringify(result).includes('token=abc'),false);
 assert.equal(result.target_download_configs_n,0);
});
test('invalid script syntax fails and does not report no APIs',()=>{
 assert.throws(()=>analyzeTransportAst('const = ;','runtime'));
});
test('biological claim upgrades remain blocked',()=>{
 const base={status:'SOURCE_AST_VERIFIED',original_RNA_downloaded:false,
 PDC_sample_join_verified:false};
 assert.equal(review(base).actual_baseURL_and_public_metadata_API_verified,false);
 assert.match(review(base).G30_D,/^STOP/);
 assert.throws(()=>review({...base,original_RNA_downloaded:true}));
 assert.throws(()=>review({...base,PDC_sample_join_verified:true}));
});
