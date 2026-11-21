import unittest
from pathlib import Path
from needleman_wunsch import needleman_wunsch, read_fasta


class AlignmentTests(unittest.TestCase):
    def test_example(self):
        r = needleman_wunsch('ACCGTGAAGCCAATAC', 'AGCGTGCAGCCAATAC')
        self.assertEqual(len(r.scores), 17)
        self.assertEqual(len(r.scores[0]), 17)
        self.assertEqual(len(r.aligned1), len(r.aligned2))
        self.assertEqual(r.path[0], (16, 16))
        self.assertEqual(r.path[-1], (0, 0))
        self.assertEqual(r.score, r.matches - r.mismatches - 2 * r.gaps)

    def test_gap(self):
        r = needleman_wunsch('ACGT', 'AGT')
        self.assertEqual(r.score, 1)
        self.assertEqual(r.aligned1, 'ACGT')
        self.assertEqual(r.aligned2, 'A-GT')

    def test_identical(self):
        r = needleman_wunsch('ACGT', 'ACGT')
        self.assertEqual((r.score, r.matches, r.mismatches, r.gaps), (4, 4, 0, 0))

    def test_uploaded_fasta(self):
        folder = Path(__file__).resolve().parent
        for name in ('480540219_genomic.fna', 'GCF_000919155.1_ViralProj244477_genomic.fna.gz'):
            header, sequence, count = read_fasta(str(folder / name))
            self.assertTrue(header)
            self.assertGreater(len(sequence), 100)
            self.assertGreaterEqual(count, 1)


if __name__ == '__main__':
    unittest.main()
