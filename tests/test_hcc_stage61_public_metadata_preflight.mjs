import test from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import {BUNDLES,ROUTES,normalizedEndpoint,isJsonMime,strictHeaders,
  verifyBundle,sourceCooccurrenceAnalysis,sourceGate,classifyHead,
  anonymousJsonSummary,review,inspect}
 from '../scripts/hcc/audit_stage61_public_metadata_preflight.mjs';

const testCode='var client=x.create({baseURL:"/api"});'+
'const q=x({url:"/conf/download",method:"GET"});'+
'const r=x({url:"/conf/download2",method:"GET"});'+
'const s=x({url:"/conf/download3",method:"GET"});';

test('exact three public source-derived paths only',()=>{
 assert.deepEqual(Object.keys(ROUTES).sort(),['download','download2','download3']);
 assert.equal(normalizedEndpoint('download'),'/api/conf/download');
 assert.throws(()=>normalizedEndpoint('OED00001234'));
 assert.throws(()=>normalizedEndpoint('private'));
});
test('both route JS source hashes remain pinned',()=>{
 assert.equal(BUNDLES.download.size,11335);
 assert.equal(BUNDLES.repository.size,108379);
 assert.ok(Object.values(BUNDLES).every(x=>/^[a-f0-9]{64}$/.test(x.sha256)));
});
test('verify complete byte count hash MIME and range fail-closed',()=>{
 const body=Buffer.from('abc'),pin={size:3,
 sha256:createHash('sha256').update(body).digest('hex')};
 const headers={'content-type':'text/javascript','content-range':'bytes 0-2/3'};
 assert.equal(verifyBundle(body,pin,206,headers).pinned_source_verified,true);
 assert.throws(()=>verifyBundle(body,{...pin,sha256:'0'.repeat(64)},206,headers));
 assert.throws(()=>strictHeaders(206,{...headers,'content-range':'bytes 0-2/30'},3));
 assert.throws(()=>strictHeaders(200,{...headers,'content-length':'4'},3));
 assert.throws(()=>strictHeaders(206,{...headers,'content-type':'text/html'},3));
});
test('both bundles must declare actual create and three GET configs',()=>{
 const a=sourceCooccurrenceAnalysis(testCode,'repository');
 assert.equal(a.invoked_create_baseurl_api_count,1);
 assert.equal(a.exact_get_download_config_count,3);
 const gate=sourceGate([a,sourceCooccurrenceAnalysis(testCode,'download')]);
 assert.equal(gate.independent_bundles_verified_n,2);
 assert.equal(gate.actual_http_transport_instance_binding_verified,false);
});
test('co-occurrence alone cannot promote client-instance binding',()=>{
 const a=sourceCooccurrenceAnalysis(testCode,'repository');
 assert.equal(a.exact_request_instance_binding_verified,false);
 assert.equal(a.actual_network_request_observed,false);
});
test('missing download request and baseURL prevent any preflight',()=>{
 const missing=sourceCooccurrenceAnalysis('x({url:"/conf/download",method:"GET"});','repository');
 assert.throws(()=>sourceGate([missing,sourceCooccurrenceAnalysis(testCode,'download')]));
});
test('inert baseURL properties cannot authorize a GET',()=>{
 const src='const config={baseURL:"/api",url:"/conf/download",method:"GET"};';
 const data=sourceCooccurrenceAnalysis(src,'repository');
 assert.equal(data.invoked_create_baseurl_api_count,0);
 assert.equal(data.exact_get_download_config_count,0);
});
test('metadata JSON HEAD needs exact safe length MIME, status, no redirect',()=>{
 assert.equal(isJsonMime('application/json; charset=utf-8'),true);
 assert.equal(classifyHead(200,{'content-type':'application/json','content-length':'31'}).metadata_get_allowlisted,true);
 assert.equal(classifyHead(200,{'content-type':'text/html','content-length':'31'}).metadata_get_allowlisted,false);
 assert.equal(classifyHead(200,{'content-type':'application/json'}).metadata_get_allowlisted,false);
 assert.equal(classifyHead(405,{'content-type':'application/json','content-length':'31'}).metadata_get_allowlisted,false);
 assert.equal(classifyHead(302,{'location':'/auth','content-type':'application/json','content-length':'31'}).metadata_get_allowlisted,false);
 assert.equal(classifyHead(200,{'content-type':'application/json','content-length':'999999'}).metadata_get_allowlisted,false);
});
test('anonymous JSON metadata excludes arbitrary patient identifiers and links',()=>{
 const bytes=Buffer.from(JSON.stringify({code:200,message:'sample PATIENT_ID 888',
  data:[{SAMPLE_ID:'P00001',name:'HCCDB25.sample.txt',url:'https://example.invalid/private'}]}));
 const x=anonymousJsonSummary(bytes);
 assert.deepEqual(x.top_level_whitelisted_metadata_keys,['code','data','message']);
 assert.equal(x.expected_triplet_exact_filename_mentions['HCCDB25.sample.txt'],1);
 assert.equal(JSON.stringify(x).includes('P00001'),false);
 assert.equal(JSON.stringify(x).includes('example.invalid'),false);
 assert.equal(x.RNA_matrix_or_patient_file_downloaded,false);
});
test('preflight only sends HEAD after both source identities validated',async()=>{
 let calls=0;
 const load=async key=>({data:Buffer.from(testCode),meta:{source:key,pinned_source_verified:true}});
 const request=async(path,method)=>{calls++;
   assert.ok(Object.values(ROUTES).includes(path));
   assert.equal(method,'HEAD');
   return {status:405,headers:{'content-type':'text/html'}};
 };
 const x=await inspect('download',load,request);
 assert.equal(calls,1);
 assert.equal(x.status,'PUBLIC_HEAD_PREFLIGHT_COMPLETE');
 assert.equal(x.public_head.reason,'HEAD_METHOD_UNSUPPORTED_HOLD');
 assert.equal(x.review.biological_analysis_permitted,false);
});
test('HEAD verified JSON permits just one bounded metadata GET',async()=>{
 let calls=[];
 const body=Buffer.from(JSON.stringify({code:200,files:['HCCDB25_mRNA_level3.txt']}));
 const load=async key=>({data:Buffer.from(testCode),meta:{source:key,pinned_source_verified:true}});
 const request=async(path,method)=>{
  calls.push(method);
  assert.equal(path,'/api/conf/download');
  if(method==='HEAD')return {status:200,headers:{'content-type':'application/json',
    'content-length':String(body.length)}};
  return {status:200,headers:{'content-type':'application/json',
    'content-length':String(body.length)},body};
 };
 const x=await inspect('download',load,request);
 assert.deepEqual(calls,['HEAD','GET']);
 assert.equal(x.public_json_metadata.expected_triplet_exact_filename_mentions['HCCDB25_mRNA_level3.txt'],1);
 assert.equal(x.review.original_RNA_triplet_verified,false);
});
test('source failures never send backend requests',async()=>{
 let requested=false;
 const load=async key=>({data:Buffer.from('const x={baseURL:"/api"};'),meta:{}});
 const request=async()=>{requested=true;throw new Error('must not run');};
 const x=await inspect('download2',load,request);
 assert.equal(requested,false);
 assert.equal(x.status,'SOURCE_IDENTITY_HOLD');
 assert.match(x.review.G30_D,/^STOP/);
});
test('unverifiable endpoint response or file values do not promote GO',()=>{
 const s={rna_original_downloaded:false,pdc_sample_join_verified:false};
 assert.equal(review(s).public_metadata_index_verified,false);
 assert.throws(()=>review({...s,rna_original_downloaded:true}));
 assert.throws(()=>review({...s,pdc_sample_join_verified:true}));
});
