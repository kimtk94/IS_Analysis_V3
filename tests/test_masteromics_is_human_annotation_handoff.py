import hashlib
import json
import unittest
from pathlib import Path


class IsHumanAnnotationHandoffTest(unittest.TestCase):
    def test_drive_handoff_manifest_and_notebook_are_pinned(self):
        repo = Path(__file__).resolve().parents[1]
        manifest_path = repo / 'projects/ischemic_stroke.human_annotation.json'
        notebook_path = repo / 'notebooks/is/IS_Phase9F_E_Human_Vascular_12GB_R3.ipynb'
        manifest = json.loads(manifest_path.read_text())
        self.assertEqual(manifest['schema_version'], 1)
        self.assertEqual(manifest['project'], 'ischemic_stroke')
        self.assertEqual(manifest['status'], 'READY_TO_RUN_NOT_EXECUTED')
        self.assertEqual(manifest['notebook_version'], 'R3_STRICT_AUTHOR_ANNOTATION')
        self.assertEqual(manifest['drive_file_id'], '1PljK4AexdM_e5H2Eew6xDbcWPIcgrCxs')
        self.assertTrue(manifest['colab_url'].endswith(manifest['drive_file_id']))
        self.assertEqual(manifest['input_size_bytes'], 3493366138)
        self.assertEqual(manifest['scientific_guardrail'], 'REFERENCE_LOCALIZATION_NOT_DISEASE_DGE')
        self.assertEqual(manifest['current_blocker'], 'EXECUTE_R3_AND_FREEZE_AUTHOR_ANNOTATION')
        self.assertIn('R3_STRICT_AUTHOR', manifest['output_drive_path'])

        digest = hashlib.sha256(notebook_path.read_bytes()).hexdigest()
        self.assertEqual(digest, manifest['notebook_sha256'])

        nb = json.loads(notebook_path.read_text())
        self.assertEqual(nb['nbformat'], 4)
        self.assertEqual(len(nb['cells']), 1)
        source = ''.join(nb['cells'][0]['source'])
        required = [
            'R3_STRICT_AUTHOR_ANNOTATION',
            'AUTHOR_ANNOTATION_COLUMN = ""',
            'AUTHOR_ANNOTATION_COLUMN_AMBIGUOUS',
            'AUTHOR_ANNOTATION_COLUMN_UNRESOLVED',
            'predominantly numeric cluster labels',
            'HUMAN_AUTHOR_ANNOTATION_FREEZE.tsv',
            'HUMAN_RUN_MANIFEST.json',
            'REFERENCE_LOCALIZATION_NOT_DISEASE_DGE',
        ]
        for token in required:
            self.assertIn(token, source)
        self.assertNotIn(
            'if (nrow(candidates)>=1) annotation_col <- candidates$column[[1]]',
            source,
        )


if __name__ == '__main__':
    unittest.main()
