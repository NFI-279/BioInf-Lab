"""Optional automated checks; run: py -m unittest discover -v"""
import unittest

from dna_analysis import base_frequencies, clean_dna, count_kmers
from ex1 import EXAMPLE as SHORT_SEQUENCE
from ex2 import EXAMPLE as LONG_SEQUENCE


class Lab1Tests(unittest.TestCase):
    def test_normalization(self):
        self.assertEqual(clean_dna('ac g\nt'), 'ACGT')

    def test_invalid_dna(self):
        with self.assertRaises(ValueError):
            clean_dna('ACNX')

    def test_exercise_one(self):
        result = base_frequencies(SHORT_SEQUENCE)
        self.assertEqual(sum(result.values()), 20)
        self.assertEqual(result, {'A': 6, 'C': 3, 'G': 7, 'T': 4})

    def test_exercise_two_length(self):
        self.assertEqual(len(LONG_SEQUENCE), 40)

    def test_overlapping_words(self):
        self.assertEqual(count_kmers('AAAA', 2), {'AA': 3})
        self.assertEqual(count_kmers('AAAA', 3), {'AAA': 2})

    def test_exercise_two_window_counts(self):
        self.assertEqual(sum(count_kmers(LONG_SEQUENCE, 2).values()), 39)
        self.assertEqual(sum(count_kmers(LONG_SEQUENCE, 3).values()), 38)
        self.assertEqual(count_kmers(LONG_SEQUENCE, 2)['CC'], 20)
        self.assertEqual(count_kmers(LONG_SEQUENCE, 3)['CCC'], 18)


if __name__ == '__main__':
    unittest.main()
