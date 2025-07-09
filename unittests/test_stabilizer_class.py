"""
import unittest
import sys
import os
sys.path.insert(0, os.path.abspath("C:/Users/sanne/QAOA_LX-mixer"))

from Stabilizer import *

#TODO need a better test such as one where we get stabilizers that have more than 1 element... not sure if compute_projector_stabilizers work in that case
class TestStabilizer(unittest.TestCase):
    def setUp(self):
        self.test_cases = [([[0b1011, 0b1100, 0b0111, 0b0000, 0b1110, 0b1001, 0b0010, 0b0101]], 4, [[11, 14, 7]], [[(1, 13)]], [[(1, 0), (1, 13)]]),
                           ([[0b11000, 0b00100, 0b0111, 0b01101, 0b1110, 0b10001]], 5, [[9, 28]], [[(1, 2), (-1, 20), (1, 25)]], [[(1, 0), (1, 25), (-1, 20), (-1, 13), (1, 2), (1, 27), (-1, 22), (-1, 15)]])]

    def test_check_if_orbit(self):
        for B, n, expected_orbit, expected_minimal_generating_set, expected_projector in self.test_cases:
            with self.subTest(B=B, n=n, expected_orbit=expected_orbit, expected_minimal_generating_set=expected_minimal_generating_set, expected_projector=expected_projector):
                stab = Stabilizer(familiy_of_graphs=None, B=B, n=n)
                stab.check_if_orbit()
                self.assertEqual(stab.orbits, expected_orbit)
    
    def test_compute_minimal_generating_sets(self): 
        for B, n, expected_orbit, expected_minimal_generating_set, expected_projector in self.test_cases:
            with self.subTest(B=B, n=n, expected_orbit=expected_orbit, expected_minimal_generating_set=expected_minimal_generating_set, expected_projector=expected_projector):
                stab = Stabilizer(familiy_of_graphs=None, B=B, n=n)
                stab.orbits = expected_orbit
                stab.compute_minimal_generating_sets()
                self.assertEqual(stab.minimal_generating_sets, expected_minimal_generating_set)
        
    
    def test_compute_projector_stabilizers(self):
        for B, n, expected_orbit, expected_minimal_generating_set, expected_projector in self.test_cases:
            with self.subTest(B=B, n=n, expected_orbit=expected_orbit, expected_minimal_generating_set=expected_minimal_generating_set, expected_projector=expected_projector):
                stab = Stabilizer(familiy_of_graphs=None, B=B, n=n)
                stab.orbits = expected_orbit
                stab.minimal_generating_sets = expected_minimal_generating_set
                stab.compute_projector_stabilizers()
                self.assertEqual(stab.projectors, expected_projector)



if __name__ == '__main__':
    unittest.main()
"""
import unittest
import sys
import os
import copy

# Add path to module
sys.path.insert(0, os.path.abspath("C:/Users/sanne/QAOA_LX-mixer"))
from Stabilizer import *

class TestStabilizer(unittest.TestCase):
    def setUp(self):
        self.test_cases = [
            {
                "B": [[0b1011, 0b1100, 0b0111, 0b0000, 0b1110, 0b1001, 0b0010, 0b0101]],
                "n": 4,
                "expected_orbits": [[11, 14, 7]],
                "expected_mgs": [[(1, 13)]],
                "expected_projectors": [[(1, 0), (1, 13)]]
            },
            {
                "B": [[0b11000, 0b00100, 0b01101, 0b10001]],
                "n": 5,
                "expected_orbits": [[9, 28]],
                "expected_mgs": [[(1, 2), (-1, 20), (1, 25)]],
                "expected_projectors": [[(1, 0), (1, 25), (-1, 20), (-1, 13), (1, 2), (1, 27), (-1, 22), (-1, 15)]]
            }
        ]

    def test_check_if_orbit(self):
        for case in self.test_cases:
            with self.subTest(case=case):
                case = copy.deepcopy(case)
                stab = Stabilizer(familiy_of_graphs=None, B=case["B"], n=case["n"])
                stab.check_if_orbit()
                self.assertEqual(stab.orbits, case["expected_orbits"])

    def test_compute_minimal_generating_sets(self):
        for case in self.test_cases:
            with self.subTest(case=case):
                case = copy.deepcopy(case)
                stab = Stabilizer(familiy_of_graphs=None, B=case["B"], n=case["n"])
                stab.orbits = case["expected_orbits"]
                stab.compute_minimal_generating_sets()
                self.assertEqual(stab.minimal_generating_sets, case["expected_mgs"])

    def test_compute_projector_stabilizers(self):
        for case in self.test_cases:
            with self.subTest(case=case):
                case = copy.deepcopy(case)
                stab = Stabilizer(familiy_of_graphs=None, B=case["B"], n=case["n"])
                stab.orbits = case["expected_orbits"]
                stab.minimal_generating_sets = case["expected_mgs"]
                stab.compute_projector_stabilizers()
                self.assertEqual(stab.projectors, case["expected_projectors"])

if __name__ == '__main__':
    unittest.main()
