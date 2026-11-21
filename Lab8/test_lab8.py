import unittest
from ex1 import FIRST, SECOND, align, statistics
from ex2 import DEMO_DNA, calculate_pattern, gc_percent, information_content


class Lab8Tests(unittest.TestCase):
    def test_alignment(self):
        result = align(FIRST, SECOND)
        self.assertEqual(result.score, 12)
        self.assertEqual(statistics(result), {'matches': 14, 'mismatches': 2, 'gaps': 0})

    def test_initialization(self):
        result = align('A', 'G')
        self.assertEqual(result.grid[0], [0, -2])
        self.assertEqual(result.grid[1][0], -2)

    def test_window_count(self):
        self.assertEqual(len(DEMO_DNA), 73)
        self.assertEqual(len(calculate_pattern('demo', DEMO_DNA).gc), 44)

    def test_gc(self):
        self.assertAlmostEqual(gc_percent('ACGT'), 50.0)

    def test_ic_formula(self):
        self.assertAlmostEqual(information_content('AAAA'), 99.99)


if __name__ == '__main__':
    unittest.main()
