import unittest
import math
import numpy as np
import sinabro
import sinabro.sinabro as legacy
import sinabro_v2 as v2
from sinabro import RobustnessComputer
SEQ = 'TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT'

class RobustnessTests(unittest.TestCase):

    def test_original_import_exposes_trajectory(self):
        self.assertIs(sinabro.Trajectory, legacy.Trajectory)

    def test_legacy_signature_keyword_names(self):
        from sinabro.core import signature_path
        from sinabro._mutate._mutate_seq_with_mutational_signature import make_full_mutational_signature, compute_mut_prob_matrix
        matrix = make_full_mutational_signature(input_file_path=signature_path(), column='SBS2')
        p = compute_mut_prob_matrix(seq=SEQ, full_mutational_signature=matrix)
        self.assertAlmostEqual(p.sum(), 1)

    def test_v2_record_and_trajectory(self):
        t = v2.Trajectory(0, v2.MutationRecord(SEQ), seed=42)
        t.autofill('random', 'max_length', maxlen=3)
        self.assertEqual(t.get_length(), 3)
        self.assertEqual(len(t.get_records()), 4)

    def test_v2_motif_api(self):
        result = v2.mutate.mutate_seq_with_mut_type('ATCAT', mut_type='T[C>T]A', seed=1)
        self.assertEqual(result.new_seq, 'ATTAT')

    def test_v2_weighted_motif_api(self):
        result = v2.mutate.mutate_seq_with_mut_types('ATCAT', ['T[C>T]A'], [1], seed=1)
        self.assertEqual(result.new_seq, 'ATTAT')

    def test_stopping_time_is_mean_preterminal_length(self):
        r = RobustnessComputer('x', SEQ, seed=42)
        value = r.compute_l_robustness(20, 'signature', 'nonsynonymous', mutational_signature='SBS2')
        self.assertEqual(value, float(np.mean([t._length - 1 for t in r.trajs['lrho']['trajs']])))

    def test_blosum_stopping_time_is_labeled_separately(self):
        r = RobustnessComputer('x', SEQ, seed=42)
        value = r.compute_l_robustness(5, 'signature', 'blosum', mutational_signature='SBS2', threshold=20)
        self.assertGreaterEqual(value, 0)
        self.assertEqual(r.trajs['lrho']['evaluation'], 'blosum_distance')

    def test_relaxed_criterion_retains_same_mutation_prefix(self):
        r = RobustnessComputer('x', SEQ, seed=42)
        strict = r.generate_trajectories(10, 'signature', 'nonsynonymous', mutational_signature='SBS2')
        relaxed = r.generate_trajectories(10, 'signature', 'blosum', mutational_signature='SBS2', threshold=20)
        for a,b in zip(strict, relaxed):
            self.assertGreaterEqual(b.get_length(), a.get_length())
            self.assertEqual(a._data, b._data[:len(a._data)])

    def test_incomplete_stopping_time_is_not_silently_averaged(self):
        r = RobustnessComputer('x', 'AAAAA', seed=42)
        with self.assertRaises(ValueError):
            r.compute_l_robustness(3, 'mut_type', 'nonsynonymous', mut_type='T[C>T]A')
        self.assertEqual(r.trajs['lrho']['n_incomplete'], 3)
        self.assertIsNone(r.trajs['lrho']['robustness'])

    def test_fixed_n_mean_matches_endpoint_translation(self):
        r = RobustnessComputer('x', SEQ, seed=42)
        value = r.compute_n_robustness(30, 'random', n=2)
        trajs = r.trajs['nrho-2']['trajs']
        expected = sum((t._data[0][1:-1].translate() == t._data[-1][1:-1].translate() for t in trajs)) / 30
        self.assertEqual(value, expected)

    def test_fixed_zero_is_one(self):
        self.assertEqual(RobustnessComputer('x', SEQ, seed=42).compute_n_robustness(3, 'random', n=0), 1)

    def test_fixed_n_blosum_and_zero_distance_score(self):
        r = RobustnessComputer('x', SEQ, seed=42)
        self.assertEqual(r.compute_n_robustness(3, 'random', n=0, phi_eval='blosum', threshold=1), 1)
        self.assertEqual(r.compute_n_robustness(3, 'random', n=0, phi_eval='blosum', by_score=True), math.inf)

    def test_fixed_n_cannot_accept_short_no_event_path(self):
        with self.assertRaises(ValueError):
            RobustnessComputer('x', 'AAAAA', seed=42).compute_n_robustness(3, 'mut_type', n=2, mut_type='T[C>T]A')

    def test_serial_parallel_paths_identical(self):
        r = RobustnessComputer('x', SEQ, seed=42)
        a = r.generate_trajectories(8, 'signature', 'nonsynonymous', mutational_signature='SBS2')
        b = r.generate_trajectories(8, 'signature', 'nonsynonymous', mutational_signature='SBS2', multiprocessing=True, n_cpu=2)
        self.assertEqual([t._data for t in a], [t._data for t in b])
        self.assertEqual([t.terminal_reason for t in a], [t.terminal_reason for t in b])

    def test_fixed_n_parallel_runs(self):
        r = RobustnessComputer('x', SEQ, seed=42)
        a = r.compute_n_robustness(8, 'random', n=2)
        b = r.compute_n_robustness(8, 'random', n=2, multiprocessing=True, n_cpu=2)
        self.assertEqual(a, b)

    def test_invalid_counts(self):
        for n in [0, -1, 1.5, True]:
            with self.assertRaises(ValueError):
                RobustnessComputer('x', SEQ).generate_trajectories(n, 'random', 'max_length')
if __name__ == '__main__':
    unittest.main()
