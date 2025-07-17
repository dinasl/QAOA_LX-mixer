import unittest
import sys
import os


# Add path to module
sys.path.insert(0, os.path.abspath("C:/Users/sanne/QAOA_LX-mixer"))
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
    
class TestSplitIntoSuborbits(unittest.TestCase):
    def setUp(self):
        self.test_cases = [
            {
                "family_of_valid_graphs" : {0b0010 : [(0, 1), (2, 7), (3, 9), (4, 13), (5, 10), (6, 8), (11, 14)],
            0b0111 : [(0, 2), (1, 7), (3, 4), (5, 14), (6, 12), (9, 13), (10, 11)],
            0b1010 : [(0, 3), (1, 9), (2, 4), (6, 11), (7, 13), (8, 14), (10, 12)],
            0b1101 : [(0, 4), (1, 13), (2, 3), (5, 8), (6, 10), (7, 9), (11, 12)],
            0b1110 : [(0, 5), (1, 10), (2, 14), (4, 8), (6, 13), (7, 11), (9, 12)],
            0b0001 : [(0, 6), (1, 8), (2, 12), (3, 11), (4, 10), (5, 13), (9, 14)],
            0b0101 : [(0, 7), (1, 2), (3, 13), (4, 9), (5, 11), (8, 12), (10, 14)],
            0b0011 : [(0, 8), (1, 6), (3, 14), (4, 5), (7, 12), (9, 11), (10, 13)],
            0b1000 : [(0, 9), (1, 3), (2, 13), (4, 7), (5, 12), (6, 14), (8, 11)],
            0b1100 : [(0, 10), (1, 5), (2, 11), (3, 12), (4, 6), (7, 14), (8, 13)],
            0b1011 : [(0, 11), (1, 14), (2, 10), (3, 6), (4, 12), (5, 7), (8, 9)],
            0b0110 : [(0, 12), (2, 6), (3, 10), (4, 11), (5, 9), (7, 8), (13, 14)],
            0b1111 : [(0, 13), (1, 4), (2, 9), (3, 7), (5, 6), (8, 10), (12, 14)],
            0b1001 : [(0, 14), (1, 11), (2, 5), (3, 8), (6, 9), (7, 10), (12, 13)],
            0b0100 : [(1, 12), (2, 8), (3, 5), (4, 14), (6, 7), (9, 10), (11, 13)]},

                "operators" : [0b1001, 0b0110, 0b111],
                "nodes" : (0, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13, 14),
                "expected_suborbits" : [{0, 12, 13, 14}, {9, 2, 5, 6}, {8, 10, 3, 7}]
            }
        ]

    def test_split_into_suborbits(self):
        for case in self.test_cases:
            with self.subTest(case=case):
                suborbits = split_into_suborbits(case["family_of_valid_graphs"], case["operators"], case["nodes"])
                self.assertEqual(suborbits, case["expected_suborbits"])


if __name__ == '__main__':
    unittest.main()
