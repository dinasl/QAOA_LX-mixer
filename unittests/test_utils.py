import unittest
import sys
import os


# Add path to module
# Add path
from utils import *

class TestFindBestCost(unittest.TestCase):
    def setUp(self):
        self.test_cases = [
            {
                "Xs": [0b1111, 0b0001, 0b1101],
                "Zs": [(1, 0), (1, 12)],
                "expected_best_Xs_reduced": [12,14,1],
                "expected_best_cost": 10
            }
        ]

    def test_find_best_cost(self):
        for case in self.test_cases:
            with self.subTest(case=case):
                best_Xs_reduced, best_cost = find_best_cost(case["Xs"], case["Zs"])
                self.assertEqual(set(best_Xs_reduced), set(case["expected_best_Xs_reduced"]))
                self.assertEqual(best_cost, case["expected_best_cost"])


if __name__ == '__main__':
    unittest.main()
