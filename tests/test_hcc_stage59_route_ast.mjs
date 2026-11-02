import test from 'node:test';
import assert from 'node:assert/strict';
import {PINNED,safeSourcePath,verifyPinnedSource,analyzeAst,review}
 from '../scripts/hcc/audit_stage59_route_ast.mjs';
import {createHash} from 'node:crypto';

test('only the two actual Stage57 public route scripts are SHA-pinned', () => {
  assert.deepEqual(Object.keys(PINNED).sort(),['download','repository']);
  assert.equal(PINNED.download.bytes,11335);
  assert.equal(PINNED.repository.bytes,108379);
  assert.ok(Object.values(PINNED).every(x=>/^[a-f0-9]{64}$/.test(x.sha256)));
});
test('sha, bytes, content range and MIME must all match', () => {
  const raw=Buffer.from('abc');
  const pin={bytes:3,sha256:createHash('sha256').update(raw).digest('hex')};
  const good={'content-type':'application/javascript','content-range':'bytes 0-2/3'};
  assert.equal(verifyPinnedSource(raw,206,good,pin).complete,true);
  assert.throws(()=>verifyPinnedSource(raw,206,good,{...pin,sha256:'0'.repeat(64)}));
  assert.throws(()=>verifyPinnedSource(raw,206,{...good,'content-range':'bytes 0-2/30'},pin));
  assert.throws(()=>verifyPinnedSource(raw,206,{...good,'content-type':'text/html'},pin));
});
test('source literal safe paths do not export private IDs or query URLs', () => {
  assert.equal(safeSourcePath('/conf/download'),'/conf/download');
  assert.equal(safeSourcePath('/api/download'),'/api/download');
  assert.equal(safeSourcePath('/system/chart/getGPointData'),'/system/chart/getGPointData');
  assert.equal(safeSourcePath('/api/download/OED012345'),null);
  assert.equal(safeSourcePath('/api/download?token=123'),null);
  assert.equal(safeSourcePath('https://evil.invalid/api/download'),null);
  assert.equal(safeSourcePath('/api/private/patient'),null);
});
test('vue component path and property do not imply HTTP API',()=>{
  const js='const a={path:"/conf/download",component:foo,route:"/repository"};';
  const v=analyzeAst(js,'download');
  assert.equal(v.parse_status,'AST_PARSED');
  assert.equal(v.observed_key_path_roles['/conf/download'].object_property_value,1);
  assert.equal(v.source_declared_request_shapes_n,0);
});
test('invoked {url,method} config is classified as source-declared only',()=>{
  const code='const x=client({url:"/api/download",method:"POST"});';
  const r=analyzeAst(code,'download');
  assert.equal(r.source_declared_request_shapes_n,1);
  assert.equal(r.source_declared_request_shapes[0].method,'POST');
  assert.equal(r.source_declared_request_shapes[0].path,'/api/download');
  assert.equal(r.actual_RNA_file_content_downloaded,false);
});
test('direct GET invocation is distinct from string in route definition',()=>{
  const code='x.get("/api/metadata"); const routes={path:"/api/metadata"};';
  const r=analyzeAst(code,'repository');
  assert.equal(r.source_declared_request_shapes_n,1);
  assert.equal(r.source_declared_request_shapes[0].structural_type,'SOURCE_DECLARED_DIRECT_CALL_ARGUMENT');
  assert.equal(r.safe_literal_structural_classes.object_property_value,1);
});
test('dynamic request config without literal method is NOT promoted',()=>{
  const code='client({url:"/api/download",method:variable});';
  const r=analyzeAst(code,'download');
  assert.equal(r.source_declared_request_shapes_n,0);
});
test('raw patient identifiers and clinical values are not copied to report',()=>{
  const js='const a={url:"/api/download/OED000111", method:"GET",patient:"P01"};';
  const r=analyzeAst(js,'repository');
  assert.equal(r.source_declared_request_shapes_n,0);
  assert.equal(JSON.stringify(r).includes('OED000111'),false);
  assert.equal(JSON.stringify(r).includes('P01'),false);
});
test('invalid source must fail parse, not falsely return no API',()=>{
  assert.throws(()=>analyzeAst('const = ;','download'));
});
test('biological promotion without expression and official PDC identity forbidden',()=>{
  const base={status:'VERIFIED_PINNED_AST_ANALYZED',original_rna_file_downloaded:false,pdc_identity_matched:false};
  assert.equal(review(base).verified_RNA_triplet,false);
  assert.match(review(base).G30_D,/^STOP/);
  assert.throws(()=>review({...base,original_rna_file_downloaded:true}));
  assert.throws(()=>review({...base,pdc_identity_matched:true}));
});
