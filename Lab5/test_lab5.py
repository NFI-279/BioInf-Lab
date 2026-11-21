"""Optional unit tests: python -m unittest -v test_lab5.py"""
import unittest
from pathlib import Path
from dna_io import load_dna
from ex1 import fragments_from_dna, weight_da, travel_distance
from ex2 import motif_frequencies, matching_motifs

ROOT = Path(__file__).resolve().parent

class Lab5Tests(unittest.TestCase):
    def test_both_real_fasta_inputs(self):
        for name in ['480540219_genomic.fna', 'GCF_000919155.1_ViralProj244477_genomic.fna.gz']:
            header, dna = load_dna(ROOT / 'data' / name)
            self.assertTrue(header)
            self.assertGreater(len(dna), 2000)

    def test_fragments(self):
        dna = 'ACGT' * 500
        a = fragments_from_dna(dna, seed=3)
        self.assertEqual(a, fragments_from_dna(dna, seed=3))
        self.assertEqual(len(a), 10)
        self.assertTrue(all(100 <= len(fragment) <= 2000 for _, fragment in a))

    def test_molecular_weight(self):
        self.assertAlmostEqual(weight_da(100), 61832.04)

    def test_migration(self):
        self.assertAlmostEqual(travel_distance(100), 10)
        self.assertAlmostEqual(travel_distance(20000), 0)
        self.assertGreater(travel_distance(200), travel_distance(2000))

    def test_overlapping_motifs(self):
        counts = motif_frequencies('AAAAAA', sizes=(3,))
        self.assertEqual(counts['AAA'], 4)
        self.assertEqual(matching_motifs(counts, 4, 4), [('AAA', 4)])

if __name__ == '__main__':
    unittest.main()
