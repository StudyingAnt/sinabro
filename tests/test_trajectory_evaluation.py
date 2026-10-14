import io
import contextlib
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices
from sinabro import Trajectory, MutInfo, EvaluationResult, ProteinIdentity, BlosumDistance
from sinabro.core import event
from sinabro.evaluation.policies import alignment_distance

class EvaluationTests(unittest.TestCase):

    def test_synonymous_identity(self):
        self.assertFalse(ProteinIdentity()('GCT', 'GCC', None).stop)
        self.assertTrue(ProteinIdentity()('GCT', 'GTC', None).stop)

    def test_stop_changes_and_stop_to_stop(self):
        p = ProteinIdentity()
        self.assertFalse(p('TAA', 'TAG', None).stop)
        self.assertTrue(p('TAA', 'CAA', None).stop)
        self.assertTrue(p('TGG', 'TGA', None).stop)

    def test_blosum_threshold_boundary_with_independent_matrix_score(self):
        matrix = substitution_matrices.load('BLOSUM80')
        expected = float(matrix['A', 'A'] - matrix['A', 'V'])
        self.assertAlmostEqual(alignment_distance('GCT', 'GTC'), expected)
        self.assertTrue(BlosumDistance(expected)('GCT', 'GTC', None).stop)
        self.assertFalse(BlosumDistance(expected + 0.01)('GCT', 'GTC', None).stop)

    def test_blosum_synonymous_zero(self):
        r = BlosumDistance(1)('GCT', 'GCC', None)
        self.assertFalse(r.stop)
        self.assertEqual(r.score, 0)

    def test_blosum_original_sequence_reference(self):
        p = BlosumDistance(100)
        matrix = substitution_matrices.load('BLOSUM80')
        score = p('GCT', 'ATC', None).score
        self.assertEqual(score, float(matrix['A', 'A'] - matrix['A', 'I']))

    def test_stop_policy_is_explicit(self):
        self.assertFalse(BlosumDistance(1000, stop_policy='score')('TGG', 'TGA', None).stop)
        r = BlosumDistance(1000, stop_policy='terminate')('TGG', 'TGA', None)
        self.assertTrue(r.stop)
        self.assertEqual(r.reason, 'stop_change')

    def test_invalid_blosum_settings(self):
        for value in [0, -1, np.nan, np.inf, True]:
            with self.subTest(threshold=value), self.assertRaises(ValueError):
                BlosumDistance(value)
        for kwargs in [dict(stop_policy='bad'), dict(mode='bad'), dict(gap_open=1)]:
            with self.assertRaises(ValueError):
                BlosumDistance(10, **kwargs)

class TrajectoryTests(unittest.TestCase):

    def test_record_snapshot_cannot_corrupt_trajectory(self):
        t = Trajectory('x', 'AGCTT')
        records = t.get_records()
        records[0].sequence = 'invalid'
        self.assertEqual(str(t.get_records()[0].sequence), 'AGCTT')

    def test_mixed_coordinate_aliases_are_rejected(self):
        t = Trajectory('x', 'AGCTT')
        with self.assertRaises(ValueError):
            t.autofill(start=1, cds_start=0)

    def test_missing_motif_arguments_are_explicit(self):
        t = Trajectory('x', 'AGCTT')
        for method in ['mut_type', 'mut_types']:
            with self.assertRaises(ValueError):
                t.add_mutated_sequence(method)

    def test_already_non_neutral_manual_path_cannot_be_recounted(self):
        t = Trajectory('x', 'AGCTT')
        t.append('AGTTT')
        with self.assertRaises(ValueError):
            t.autofill(condition='nonsynonymous')

    def test_first_nonsynonymous_counts_zero(self):
        t = Trajectory('x', 'AATGT')
        with patch('sinabro.engine.random_mutation', return_value=event('AATGT', 3, 'C', '[G>C]', 1)):
            t.autofill(condition='nonsynonymous', method='random')
        self.assertEqual(t._length, 1)
        self.assertEqual(t.accepted_steps, 0)
        self.assertEqual(t.terminal_reason, 'nonsynonymous')

    def test_completed_trajectory_cannot_be_recounted(self):
        t = Trajectory('x', 'AATGT', seed=1)
        t.autofill(condition='nonsynonymous')
        with self.assertRaises(ValueError):
            t.autofill(condition='nonsynonymous')

    def test_reserved_custom_terminal_reason(self):
        t = Trajectory('x', 'AGCTT', seed=2)
        with self.assertRaises(ValueError):
            t.autofill(condition='evaluation', evaluator=lambda *args: EvaluationResult(True, 0, 'max_steps'))

    def test_synonymous_then_nonsynonymous(self):
        events = [event('AGCTT', 3, 'C', '[T>C]', 1), event('AGCCT', 2, 'T', '[C>T]', 1)]
        t = Trajectory('x', 'AGCTT')
        with patch('sinabro.engine.random_mutation', side_effect=events):
            t.autofill(condition='nonsynonymous', method='random')
        self.assertEqual(t.accepted_steps, 1)
        self.assertEqual(t.synonymous_steps, 1)

    def test_blosum_allows_missense_below_threshold(self):
        matrix = substitution_matrices.load('BLOSUM80')
        distance = float(matrix['A', 'A'] - matrix['A', 'V'])
        t = Trajectory('x', 'AGCTT')
        first = event('AGCTT', 2, 'T', '[C>T]', 1)
        with patch('sinabro.engine.random_mutation', return_value=first):
            t.autofill(condition='evaluation', threshold=distance + 0.1, max_steps=1)
        self.assertEqual(t._length, 1)
        self.assertEqual(t.synonymous_steps, 0)
        self.assertEqual(t.terminal_reason, 'max_steps')
        self.assertIsNone(t.accepted_steps)
        self.assertEqual(t.evaluation, 'blosum_distance')

    def test_blosum_hit_counts_criterion_neutral_steps(self):
        events = [event('AGCTT', 3, 'C', '[T>C]', 1), event('AGCCT', 2, 'T', '[C>T]', 1)]
        t = Trajectory('x', 'AGCTT')
        with patch('sinabro.engine.random_mutation', side_effect=events):
            t.autofill(condition='evaluation', threshold=1)
        self.assertEqual(t.terminal_reason, 'evaluation')
        self.assertEqual(t.accepted_steps, 1)
        self.assertEqual(t.score, alignment_distance('GCT', 'GTC'))

    def test_custom_evaluator(self):
        t = Trajectory('x', 'AGCTT', seed=2)
        seen = []

        def evaluator(original, current, info):
            seen.append((original, current, info.idx_target))
            return EvaluationResult(True, 3, 'custom_threshold')
        t.autofill(condition='evaluation', evaluator=evaluator)
        self.assertEqual(seen[0][0], 'GCT')
        self.assertEqual(len(seen[0][1]), 3)
        self.assertEqual(t.terminal_reason, 'custom_threshold')
        self.assertEqual(t.accepted_steps, 0)

    def test_invalid_custom_return(self):
        for value in ['bad', EvaluationResult(True, np.nan), EvaluationResult(1, 0)]:
            t = Trajectory('x', 'AGCTT', seed=2)
            with self.assertRaises((ValueError, TypeError)):
                t.autofill(condition='evaluation', evaluator=lambda *args: value)

    def test_no_motif_is_not_a_negative_count(self):
        t = Trajectory('x', 'AAAAA', seed=42)
        t.autofill(condition='nonsynonymous', method='mut_type', mut_type='T[C>T]A')
        self.assertEqual(t._length, 0)
        self.assertIsNone(t.accepted_steps)
        self.assertEqual(t.terminal_reason, 'no_available_event')

    def test_step_cap_is_censoring(self):
        t = Trajectory('x', 'AGCTT', seed=42)
        t.autofill(condition='evaluation', evaluator=lambda *args: False, max_steps=3)
        self.assertEqual(t._length, 3)
        self.assertEqual(t.terminal_reason, 'max_steps')
        self.assertIsNone(t.accepted_steps)

    def test_fixed_zero_and_exact_length(self):
        for n in [0, 1, 10]:
            t = Trajectory('x', 'AGCTT', seed=42)
            t.autofill(condition='max_length', max_length=n)
            self.assertEqual(t._length, n)
            self.assertEqual(t.terminal_reason, 'max_length')
            self.assertIsNone(t.accepted_steps)

    def test_recurrent_hit_and_back_mutation_are_allowed(self):
        events = [event('AGCTT', 3, 'C', '[T>C]', 1), event('AGCCT', 3, 'T', '[C>T]', 1)]
        t = Trajectory('x', 'AGCTT')
        with patch('sinabro.engine.random_mutation', side_effect=events):
            t.autofill(max_length=2)
        self.assertEqual(t._data[0], t._data[-1])
        self.assertEqual(t._length, 2)

    def test_trajectory_seed_reproduces_signature_path(self):
        paths = []
        for _ in range(2):
            t = Trajectory('x', 'TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT', seed=123)
            t.autofill(condition='nonsynonymous', method='signature', mutational_signature='SBS2')
            paths.append((t._data, t._hgvs_mrnas, t._hgvs_aa, t.terminal_reason))
        self.assertEqual(paths[0], paths[1])

    def test_seed_cannot_be_reinitialized_per_event(self):
        t = Trajectory('x', 'AGCTT', seed=42)
        with self.assertRaises(ValueError):
            t.add_mutated_sequence('random', seed=2)

    def test_no_debug_output(self):
        with contextlib.redirect_stdout(io.StringIO()) as output:
            t = Trajectory('x', 'AGCTT', seed=42)
            t.autofill(condition='nonsynonymous')
        self.assertEqual(output.getvalue(), '')

    def test_fasta_export(self):
        t = Trajectory('x', 'AGCTT', seed=42)
        t.autofill(max_length=2)
        with tempfile.TemporaryDirectory() as directory:
            path = t.save(Path(directory) / 'new', prefix='x')
            records = list(SeqIO.parse(path, 'fasta'))
            self.assertEqual(len(records), 3)
            self.assertEqual(records[-1].seq, t._data[-1])

    def test_unknown_modes_and_missing_settings_fail_explicitly(self):
        for kwargs in [dict(condition='bad'), dict(method='bad'), dict(condition='evaluation'), dict(condition='evaluation', measure='blastp', threshold=10), dict(max_steps=0)]:
            t = Trajectory('x', 'AGCTT')
            with self.assertRaises((ValueError, TypeError)):
                t.autofill(**kwargs)
if __name__ == '__main__':
    unittest.main()
