import itertools
import unittest
import numpy as np
from Bio.Seq import Seq, MutableSeq
from sinabro import Trajectory, MutationRecord
from sinabro.core import get_codon, dna
from sinabro._helper._helper import _get_idx_from_hgvs_mrna

class CoordinateTests(unittest.TestCase):

    def test_all_codons_and_all_nine_substitutions_against_translation(self):
        for letters in itertools.product('ACGT', repeat=3):
            codon = ''.join(letters)
            for pos in range(3):
                for alt in 'ACGT':
                    if alt == codon[pos]:
                        continue
                    with self.subTest(codon=codon, pos=pos, alt=alt):
                        new = codon[:pos] + alt + codon[pos + 1:]
                        t = Trajectory('oracle', 'A' + codon + 'T')
                        t.append('A' + new + 'T')
                        same = str(Seq(codon).translate()) == str(Seq(new).translate())
                        self.assertEqual(t._hgvs_aa[-1] == '.', same)
                        self.assertEqual(t._hgvs_mrnas[-1], f'c.{pos + 1}{codon[pos]}>{alt}')
                        if not same:
                            self.assertTrue(t._hgvs_aa[-1].startswith('p.1'))

    def test_codon_positions_in_multiple_codons(self):
        seq = 'AATGGCCTAAC'
        for idx in range(1, 10):
            self.assertEqual(get_codon(seq, idx), ['ATG', 'GCC', 'TAA'][(idx - 1) // 3])

    def test_third_position_annotation_is_not_next_codon(self):
        t = Trajectory('x', 'AATGT')
        t.append('AATCT')
        self.assertEqual(t._hgvs_aa[-1], 'p.1Met>Ile')

    def test_plain_cds_coordinates(self):
        t = Trajectory('x', 'ATGGCCTAA', cds_start=0, cds_end=9)
        t.append('ATCGCCTAA')
        self.assertEqual(t._hgvs_mrnas[-1], 'c.3G>C')
        self.assertEqual(t._hgvs_aa[-1], 'p.1Met>Ile')

    def test_nonstandard_offset(self):
        t = Trajectory('x', 'CCATGGCCTAA', cds_start=2, cds_end=11)
        t.append('CCATCGCCTAA')
        self.assertEqual(t._hgvs_aa[-1], 'p.1Met>Ile')

    def test_flanks_are_not_mutation_targets(self):
        for idx in [0, 4]:
            t = Trajectory('x', 'AATGT')
            s = list('AATGT')
            s[idx] = 'C'
            with self.assertRaises(ValueError):
                t.append(''.join(s))

    def test_codon_out_of_bounds(self):
        for idx in [0, 10]:
            with self.assertRaises(ValueError):
                get_codon('AATGGCCTAAC', idx)

    def test_input_validation(self):
        for seq in ['', 'ATN', 'AAAT', 123, None]:
            with self.subTest(seq=seq), self.assertRaises((ValueError, TypeError)):
                Trajectory('x', seq)

    def test_lowercase_and_input_types(self):
        for seq in ['aatgt', Seq('aatgt'), MutableSeq('aatgt')]:
            self.assertEqual(str(Trajectory('x', seq)._original_seq), 'AATGT')

    def test_hgvs_parser(self):
        self.assertEqual(_get_idx_from_hgvs_mrna('c.3G>C'), 3)
        self.assertEqual(_get_idx_from_hgvs_mrna('c.3G>C', offset=0), 2)
        for invalid in ['c.0G>C', '.', 'c.3del']:
            with self.assertRaises(ValueError):
                _get_idx_from_hgvs_mrna(invalid)

    def test_append_pop_invariants(self):
        t = Trajectory('x', 'AGCTT')
        t.append('AGCCT')
        self.assertEqual(t._length, 1)
        for arr in [t._data, t._hgvs_mrnas, t._hgvs_aa, t._mut_types, t._notes]:
            self.assertEqual(len(arr), 2)
        t.pop()
        self.assertEqual(t._length, 0)
        self.assertEqual(t._data, [Seq('AGCTT')])
        with self.assertRaises(IndexError):
            t.pop()

    def test_invalid_manual_steps(self):
        t = Trajectory('x', 'AGCTT')
        for seq in ['AGCTT', 'AGTTT', 'AGCTTC', 'ATCAT']:
            if seq == 'AGTTT':
                continue
            with self.assertRaises(ValueError):
                t.append(seq)

    def test_add_record(self):
        t = Trajectory('x', 'AGCTT')
        t.add_record(MutationRecord('AGCCT'))
        self.assertEqual(t._length, 1)
if __name__ == '__main__':
    unittest.main()
