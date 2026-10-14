import unittest
import numpy as np
from sinabro import MutationRecord
from sinabro_v2.metrics.dist import hamming, compute_distance
from sinabro_v2.metrics.dnds import _has_inframe_stop, compute_dnds

class MetricsTests(unittest.TestCase):

    def test_hamming_is_amino_acid_distance(self):
        self.assertEqual(hamming('GCT', 'GCC'), 0)
        self.assertEqual(hamming('GCT', 'GTC'), 1)

    def test_invalid_translation(self):
        with self.assertRaises(ValueError):
            hamming('GC', 'GCT')

    def test_one_record_distance_is_explicit(self):
        with self.assertRaises(ValueError):
            compute_distance([MutationRecord('AGCTT')])
        self.assertEqual(compute_distance([MutationRecord('AGCTT')], threshold_hit=True), 0)

    def test_distance_includes_or_excludes_terminal(self):
        records = [MutationRecord('AGCTT'), MutationRecord('AGCCT'), MutationRecord('AGTCT')]
        self.assertEqual(compute_distance(records), 0)
        self.assertEqual(compute_distance(records, threshold_hit=True), 1)

    def test_inframe_stop_search_checks_all_codons(self):
        self.assertTrue(_has_inframe_stop('ATAACCTAA'))
        self.assertFalse(_has_inframe_stop('ATAACC'))

    def test_dnds_premature_stop_is_nan(self):
        records = [MutationRecord('AATGTGGTAA' + 'T'), MutationRecord('AATGTGATAA' + 'T')]
        (dn, ds, omega, premature, stoploss) = compute_dnds(records, threshold_hit=True)
        self.assertTrue(premature)
        self.assertFalse(stoploss)
        self.assertTrue(np.isnan(dn) and np.isnan(ds) and np.isnan(omega))

    def test_dnds_requires_terminal_stop_and_complete_records(self):
        with self.assertRaises(ValueError):
            compute_dnds([MutationRecord('AGCTT')], threshold_hit=True)
        with self.assertRaises(ValueError):
            compute_dnds([])

    def test_valid_pairwise_dnds_synonymous_change(self):
        records = [MutationRecord('A'+'ATG'+'GCT'*10+'TAA'+'T'),
                   MutationRecord('A'+'ATG'+'GCC'+'GCT'*9+'TAA'+'T')]
        (dn, ds, omega, premature, stoploss) = compute_dnds(records, threshold_hit=True)
        self.assertEqual(dn, 0)
        self.assertGreater(ds, 0)
        self.assertEqual(omega, 0)
        self.assertFalse(premature or stoploss)

    def test_saturated_short_dnds_is_undefined_not_negative(self):
        records = [MutationRecord('AATGGCTTAAT'), MutationRecord('AATGGCCTAAT')]
        dn,ds,omega,premature,stoploss = compute_dnds(records,threshold_hit=True)
        self.assertTrue(np.isnan(dn) and np.isnan(ds) and np.isnan(omega))
        self.assertFalse(premature or stoploss)
if __name__ == '__main__':
    unittest.main()
