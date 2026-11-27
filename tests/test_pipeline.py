import argparse
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import lt_metabolism as pipeline
from smoke_test import make_toy


class InputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def populate(self):
        for sample, taxa in pipeline.EXPECTED.items():
            for taxon in taxa:
                old = 'Nitratidesulfovibrio' if taxon == 'Cupidesulfovibrio' else taxon
                (self.root/f'{sample}--{old}.faa').write_text('>protein 1\nMKTLQEVLNQ\n')

    def args(self, **kwargs):
        defaults = dict(manifest=None, proteins=self.root)
        defaults.update(kwargs)
        return argparse.Namespace(**defaults)

    def test_expected_major_bins_and_old_taxon(self):
        self.populate()
        rows = pipeline.inventory(self.args())
        self.assertEqual(len(rows), 10)
        self.assertIn('LT5_Cupidesulfovibrio', [r['model_id'] for r in rows])

    def test_missing_and_duplicate_bins_rejected(self):
        self.populate()
        (self.root/'LT1--Geobacter.faa').unlink()
        with self.assertRaisesRegex(ValueError, 'Missing expected'):
            pipeline.inventory(self.args())
        (self.root/'LT1-Geobacter.faa').write_text('>a\nMKLQ\n')
        (self.root/'LT1-Geobacter-copy.faa').write_text('>a\nMKLQ\n')
        with self.assertRaisesRegex(ValueError, 'Multiple files'):
            pipeline.inventory(self.args())

    def test_no_cross_sample_communities(self):
        self.populate()
        groups = pipeline.communities(pipeline.inventory(self.args()), True)
        self.assertEqual(len(groups), 9)
        for group, members in groups.items():
            self.assertTrue(all(m.startswith(group[:3]+'_') for m in members))
        self.assertNotIn('LT1', groups)
        self.assertNotIn('LT2', groups)

    def test_proteins_and_header_mapping(self):
        source = self.root/'path with spaces.faa'
        source.write_text('>original|id one\nMKTL*\n>original|id two\nMQLE\n')
        out, mapping = self.root/'stage.faa', self.root/'mapping.tsv'
        self.assertEqual(pipeline.normalize_proteins(source, out, mapping), 2)
        self.assertIn('gene_0000001\nMKTL\n', out.read_text())
        self.assertEqual(pipeline.read_tsv(mapping)[1]['original_header'], 'original|id two')
        source.write_text('>contig\n'+'ACGT'*30+'\n')
        with self.assertRaisesRegex(ValueError, 'Looks like DNA'):
            pipeline.normalize_proteins(source, out, mapping)

    def test_medium_requires_explicit_choice_and_excludes_oxygen(self):
        db = self.root/'media.tsv'
        db.write_text('medium\tcompound\nanaerobic\tac\naerobic\to2\n')
        args = self.args(media_db=db, media='anaerobic', exploratory=False, models_only=False, plan_only=False)
        self.assertEqual(pipeline.media_settings(args), {'anaerobic': ['ac']})
        args.media = 'aerobic'
        with self.assertRaisesRegex(ValueError, 'contains o2'):
            pipeline.media_settings(args)
        args.media_db, args.media = None, None
        with self.assertRaisesRegex(ValueError, 'Supply --media-db'):
            pipeline.media_settings(args)

    def test_zero_exit_without_sbml_is_failure(self):
        for folder in ('inputs', 'models', 'work', 'logs'):
            (self.root/folder).mkdir()
        source = self.root/'proteins.faa'
        source.write_text('>a\nMKTLQ\n')
        args = self.args(outdir=self.root, solver='scip', resume=False, threads=1)
        with patch.object(pipeline, 'run_command'):
            with self.assertRaisesRegex(RuntimeError, 'produced no SBML'):
                pipeline.reconstruct({'model_id': 'LT1_Geobacter', 'fasta': str(source)}, args, {})
        self.assertFalse((self.root/'models'/'LT1_Geobacter.xml').exists())


class RealIntegrationTests(unittest.TestCase):
    def test_analysis_exports_crossfeeding_and_marks_infeasibility(self):
        from reframed import set_default_solver
        set_default_solver('scip')
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/'models').mkdir()
            paths = [root/'models'/'LT3_Geobacter.xml', root/'models'/'LT3_Azospira.xml']
            make_toy(paths[0], True)
            make_toy(paths[1], False)
            rows = [{'sample': 'LT3', 'model_id': p.stem} for p in paths]
            db = root/'media.tsv'
            db.write_text('medium\tcompound\nfeed\tac\nfeed\th2o\nstarve\th2o\n')
            args = argparse.Namespace(outdir=root, media_db=db, exploratory=False, solver='scip', pairs=False)
            with redirect_stdout(io.StringIO()):
                pipeline.analyze(rows, paths, args, {'feed': ['ac', 'h2o'], 'starve': ['h2o']})
            output = next((root/'analysis').iterdir())
            positive = pipeline.read_tsv(output/'positive_exchanges.tsv')
            self.assertTrue(any(r['donor'] == 'LT3_Geobacter' and r['receiver'] == 'LT3_Azospira'
                                and r['compound'] == 'M_for_e' for r in positive))
            status = pipeline.read_tsv(output/'run_status.tsv')
            self.assertEqual(len(status), 4)
            self.assertTrue(all(r['result'] == 'skipped_no_joint_growth' for r in status if r['medium'] == 'starve'))
            self.assertTrue(json.loads((output/'settings.json').read_text())['MUS_solver_isolation'])


if __name__ == '__main__':
    unittest.main()
