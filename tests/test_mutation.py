import csv
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch, Mock
import numpy as np
from Bio.Seq import Seq, MutableSeq
from sinabro.core import TRIPLETS, BASES, probability_matrix, full_signature, signature_path, signature_mutation, motif_mutation, random_mutation, motif_indices, expand_motif, reverse_motif, _load_signature

def signature_fixture(path, active=None):
    active = active or {'T[C>T]A': 1}
    with open(path, 'w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['Type', 'custom'])
        for ref in 'CT':
            for alt in 'ACGT':
                if alt == ref:
                    continue
                for left in 'ACGT':
                    for right in 'ACGT':
                        label = f'{left}[{ref}>{alt}]{right}'
                        writer.writerow([label, active.get(label, 0)])

class MutationTests(unittest.TestCase):

    def test_hand_computed_joint_distribution(self):
        matrix = np.zeros((4, 64))
        j = TRIPLETS.index('TCT')
        matrix[BASES.index('T'), j] = 2
        matrix[BASES.index('A'), j] = 6
        p = probability_matrix('ATCTATCTG', matrix)
        expected = np.zeros((4, 9))
        expected[3, [2, 6]] = 1 / 8
        expected[0, [2, 6]] = 3 / 8
        np.testing.assert_allclose(p, expected)

    def test_zero_mass_is_finite_zero(self):
        p = probability_matrix('ATCTATCTG', np.zeros((4, 64)))
        self.assertTrue(np.isfinite(p).all())
        self.assertEqual(p.sum(), 0)

    def test_invalid_matrix(self):
        for matrix in [np.ones((4, 63)), np.full((4, 64), -1), np.full((4, 64), np.nan), np.ones((4, 64))]:
            with self.assertRaises(ValueError):
                probability_matrix('ATCTATCTG', matrix)

    def test_packaged_profiles(self):
        for reference in ['GRCh37', 'GRCh38', 'mm9', 'mm10', 'rn6']:
            matrix = full_signature(signature_path(genome_ref=reference), 'SBS2')
            self.assertEqual(matrix.shape, (4, 64))
            self.assertAlmostEqual(matrix.sum(), 1)

    def test_all_packaged_signatures_normalize(self):
        path = signature_path()
        with open(path) as stream:
            names = next(csv.reader(stream, delimiter='\t'))[1:]
        self.assertEqual(len(names), 79)
        for name in names:
            with self.subTest(signature=name):
                matrix = full_signature(path, name)
                self.assertAlmostEqual(matrix.sum(), 1)
                p = probability_matrix('TATGCTGACTCGGTCATCGATCGGTTCTCATTGAT', matrix)
                self.assertTrue(np.isfinite(p).all())
                self.assertAlmostEqual(p.sum(), 1)
                self.assertEqual(p[:, 0].sum() + p[:, -1].sum(), 0)

    def test_signature_cache_is_protected_from_caller_mutation(self):
        a = full_signature(signature_path(), 'SBS2')
        a[:] = 0
        self.assertAlmostEqual(full_signature(signature_path(), 'SBS2').sum(), 1)

    def test_signature_cache_hits(self):
        _load_signature.cache_clear()
        full_signature(signature_path(), 'SBS2')
        full_signature(signature_path(), 'SBS2')
        self.assertEqual(_load_signature.cache_info().hits, 1)

    def test_custom_strand_expansion(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'sig.tsv'
            signature_fixture(p)
            matrix = full_signature(p, 'custom', 0.8)
            self.assertEqual(matrix[3, TRIPLETS.index('TCA')], 0.8)
            self.assertAlmostEqual(matrix[0, TRIPLETS.index('TGA')], 0.2)

    def test_single_possible_signature_event(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'sig.tsv'
            signature_fixture(p)
            result = signature_mutation('ATCAT', mutational_signature='custom', custom_signature_path=p, column='custom', seed=42)
            self.assertEqual(result.new_seq, 'ATTAT')
            self.assertEqual(result.idx_target, 2)
            self.assertEqual(result.hgvs_mrna, 'c.2C>T')

    def test_no_signature_event(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'sig.tsv'
            signature_fixture(p)
            result = signature_mutation('AAAAA', mutational_signature='custom', custom_signature_path=p, column='custom', seed=42)
            self.assertEqual(result.e, 1)
            self.assertEqual(result.new_seq, 'AAAAA')
            self.assertEqual(result.reason, 'no_available_event')

    def test_invalid_signature_column_or_weights(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'sig.tsv'
            signature_fixture(p)
            with self.assertRaises(ValueError):
                full_signature(p, 'not_found')
            signature_fixture(p, {'T[C>T]A': -1})
            with self.assertRaises(ValueError):
                full_signature(p, 'custom')

    def test_invalid_custom_signature_layout(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'sig.tsv'
            p.write_text('Type\tcustom\nT[C>T]A\t1\n')
            with self.assertRaises(ValueError):
                full_signature(p)

    def test_signature_bias_changes_real_strand_probability(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'sig.tsv'
            signature_fixture(p)
            a = probability_matrix('ATCAGATGA', full_signature(p, 'custom', 0.8))
            b = probability_matrix('ATCAGATGA', full_signature(p, 'custom', 0.2))
            self.assertFalse(np.allclose(a, b))
            self.assertAlmostEqual(a[3, 2], 0.8)

    def test_overlapping_motifs(self):
        self.assertEqual(motif_indices('CCCCC', 'CC[C>T]'), [0, 1, 2])

    def test_iupac_expansion(self):
        self.assertEqual(set(expand_motif('N[C>T]')), {f'{b}[C>T]' for b in 'ACGT'})
        self.assertEqual(len(expand_motif('NN[C>T]')), 16)

    def test_reverse_motif_is_involution(self):
        for motif in ['T[C>T]A', 'NN[C>T]', '[A>G]']:
            self.assertEqual(reverse_motif(reverse_motif(motif)), motif)

    def test_weighted_opportunities_and_overlap(self):
        mock = Mock()
        mock.choice.return_value = 0
        with patch('sinabro.core.rng_from', return_value=mock):
            motif_mutation('CCCCC', ['CC[C>T]', '[C>A]'], [2, 1], strand='single', start=0, end=5)
        probs = mock.choice.call_args.kwargs['p']
        np.testing.assert_allclose(probs, [2 / 11] * 3 + [1 / 11] * 5)

    def test_target_filter_not_motif_start(self):
        result = motif_mutation('CCCT', ['CC[C>T]'], strand='single', start=2, end=3, seed=1)
        self.assertEqual(result.idx_target, 2)
        self.assertEqual(result.new_seq, 'CCTT')

    def test_invalid_motif_weights(self):
        for weights in [[-1], [np.nan], [0], [1, 2]]:
            with self.assertRaises(ValueError):
                motif_mutation('ATCAT', ['T[C>T]A'], weights)

    def test_invalid_mutation_types(self):
        for motif in ['[C>C]', 'invalid', 'T[CC>T]']:
            with self.assertRaises(ValueError):
                expand_motif(motif)

    def test_sequence_types_and_immutability(self):
        for constructor in [str, Seq, MutableSeq]:
            seq = constructor('ATCAT')
            before = str(seq)
            result = random_mutation(seq, seed=42)
            self.assertIsInstance(result.new_seq, constructor)
            self.assertEqual(str(seq), before)

    def test_random_never_changes_flanks_or_selects_original_base(self):
        rng = np.random.default_rng(42)
        for _ in range(100):
            seq = 'ATCAT'
            result = random_mutation(seq, rng=rng)
            self.assertEqual(result.new_seq[0], seq[0])
            self.assertEqual(result.new_seq[-1], seq[-1])
            self.assertNotEqual(result.new_seq[result.idx_target], seq[result.idx_target])

    def test_seed_repeatability(self):
        a = random_mutation('ATCAT', seed=42)
        b = random_mutation('ATCAT', seed=42)
        self.assertEqual(a, b)

    def test_invalid_rng_and_bounds(self):
        for kwargs in [dict(rng=42), dict(seed=42, rng=np.random.default_rng(1)), dict(start=3, end=2), dict(start=-1), dict(start=True)]:
            with self.assertRaises((ValueError, TypeError)):
                random_mutation('ATCAT', **kwargs)
if __name__ == '__main__':
    unittest.main()
