import json
import tempfile
import unittest
from pathlib import Path

from masteromics.architecture import initialize, inspect
from masteromics.is_binding import bind_frozen_prefix, emit_frozen_stage, materialize_frozen_prefix


class IsBindingTest(unittest.TestCase):
    def make_migration(self, root: Path):
        migration = root / 'migration'
        (migration / 'stages').mkdir(parents=True)
        stage_ids = [
            'acquisition','source_qc','normalize','gwas_loci','finemap',
            'cross_ancestry','molecular_coloc','mechanism','celltype',
            'human_annotation','evidence','report',
        ]
        ready = set(stage_ids[:9])
        blocked = stage_ids[9:]
        source = root / 'source.txt'
        source.write_text('fixture\n')
        import hashlib
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        for sid in stage_ids:
            payload = {
                'schema_version':1,
                'project':'ischemic_stroke',
                'stage_id':sid,
                'migration_status':'MIGRATED_FROZEN' if sid in ready else ('BLOCKED_AUTHOR_ANNOTATION' if sid=='human_annotation' else 'BLOCKED_DEPENDENCY'),
                'scientific_status':'FROZEN_RESULTS_NOT_RECOMPUTED' if sid in ready else 'NOT_EXECUTED',
                'source_artifacts':[] if sid in blocked else [{
                    'path':str(source),
                    'size_bytes':source.stat().st_size,
                    'sha256':digest,
                }],
            }
            (migration / 'stages' / f'{sid}.json').write_text(json.dumps(payload)+'\n')
        (migration / 'IS_STAGE_MIGRATION.json').write_text(json.dumps({
            'schema_version':1,
            'project':'ischemic_stroke',
            'catalog_stage_count':12,
            'migrated_frozen_stages':9,
            'blocked_stages':blocked,
            'masteromics_binding_status':'MIGRATION_PROVENANCE_REGISTERED_NOT_EXECUTABLE',
        })+'\n')
        (migration / 'IS_STAGE_MIGRATION.tsv').write_text('stage_id\tmigration_status\n')
        return migration

    def test_bind_prefix_leaves_final_three_unbound(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hub = root / 'hub'
            initialize(hub)
            project = hub / 'ischemic_stroke' / 'project.json'
            migration = self.make_migration(root)
            out = hub / 'ischemic_stroke' / 'config' / 'project.frozen.json'
            result = bind_frozen_prefix(project, migration, out)
            self.assertEqual(result['bound_stage_count'], 9)
            self.assertEqual(result['unbound_stages'], ['human_annotation','evidence','report'])
            cfg = json.loads(out.read_text())
            assessment = inspect(cfg)
            self.assertEqual(assessment['unbound_stages'], ['human_annotation','evidence','report'])
            bindings = {s['id']:s['binding'] for s in cfg['stages']}
            self.assertIsNotNone(bindings['celltype'])
            self.assertIsNone(bindings['human_annotation'])
            self.assertEqual(cfg['data_review_status'], 'FROZEN_BASELINE_AUDITED')

    def test_emit_verifies_source_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            migration = self.make_migration(root)
            stage = migration / 'stages' / 'finemap.json'
            out = root / 'out.json'
            emitted = emit_frozen_stage(stage, out)
            self.assertEqual(emitted['parity_status'], 'SOURCE_IDENTITY_VERIFIED')
            self.assertEqual(emitted['execution_mode'], 'FROZEN_PASSTHROUGH')
            source = Path(emitted['source_artifacts'][0]['path'])
            source.write_text('changed\n')
            with self.assertRaisesRegex(ValueError, 'size changed|SHA256 changed'):
                emit_frozen_stage(stage, root / 'out2.json')

    def test_materialize_verifies_nine_and_keeps_three_blocked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            migration = self.make_migration(root)
            out = root / 'materialized'
            result = materialize_frozen_prefix(migration, out)
            self.assertEqual(result['verified_stage_count'], 9)
            self.assertEqual(result['blocked_stages'], ['human_annotation','evidence','report'])
            self.assertEqual(result['full_blueprint_status'], 'BLOCKED_UNBOUND_ADAPTERS')
            self.assertEqual(result['scientific_status'], 'FROZEN_RESULTS_NOT_RECOMPUTED')
            self.assertEqual(
                result['verified_stages'],
                ['acquisition','source_qc','normalize','gwas_loci','finemap',
                 'cross_ancestry','molecular_coloc','mechanism','celltype'],
            )
            for stage_id in result['verified_stages']:
                payload = json.loads((out / f'{stage_id}.json').read_text())
                self.assertEqual(payload['parity_status'], 'SOURCE_IDENTITY_VERIFIED')
                self.assertEqual(payload['execution_mode'], 'FROZEN_PASSTHROUGH')
            summary = json.loads((out / 'FROZEN_PREFIX_VERIFICATION.json').read_text())
            self.assertEqual(summary['verified_stage_count'], 9)
            with self.assertRaisesRegex(ValueError, 'refusing overwrite'):
                materialize_frozen_prefix(migration, out)

    def test_bind_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hub = root / 'hub'
            initialize(hub)
            project = hub / 'ischemic_stroke' / 'project.json'
            migration = self.make_migration(root)
            out = hub / 'ischemic_stroke' / 'config' / 'project.frozen.json'
            bind_frozen_prefix(project, migration, out)
            with self.assertRaisesRegex(ValueError, 'refusing overwrite'):
                bind_frozen_prefix(project, migration, out)


if __name__ == '__main__':
    unittest.main()
