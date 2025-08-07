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
# Add path
from Stabilizer import *
from Mixer import Orbit

class TestStabilizer(unittest.TestCase):
    def setUp(self):
        self.test_cases = [
            {
                "B": [0b1011, 0b1100, 0b0111, 0b0000, 0b1110, 0b1001, 0b0010, 0b0101],
                "n": 4,
                "orbit_dictionary": {(0,1,2,3,4,5,6,7): Orbit(Xs={11, 14, 7})},
                "expected_orbits": {(0,1,2,3,4,5,6,7) : [11, 14, 7]},
                "expected_mgs": {(0,1,2,3,4,5,6,7):[(1, 13)]},
                "expected_projectors": {(0,1,2,3,4,5,6,7): [(1, 0), (1, 13)]}
            },
            {
                "B": [0b11000, 0b00100, 0b01101, 0b10001],
                "n": 5,
                "orbit_dictionary": {(0,1,2,3): Orbit(Xs={9, 28})},
                "expected_orbits": {(0,1,2,3): [9, 28]},
                "expected_mgs": {(0,1,2,3) : [(1, 2), (-1, 20), (1, 25)]},
                "expected_projectors": {(0,1,2,3) : [(1, 0), (1, 25), (-1, 20), (-1, 13), (1, 2), (1, 27), (-1, 22), (-1, 15)]}
            }
            # }, #TODO the one below is not working....
            # {
            #     "B": [0b1110, 0b1100, 0b1001, 0b0100, 0b0011],
            #     "n": 4,
            #     "orbit_dictionary": {(0,2,3,4): Orbit(Xs={10, 13, 7}), (0,1): Orbit(Xs={2}), (1,2): Orbit(Xs={5}), (1,3): Orbit(Xs={8}), (1,4): Orbit(Xs={15})},
            #     "expected_orbits": {(0,1,2,3,4,5,6,7) : [11, 14, 7]},
            #     "expected_mgs": {(0,1,2,3,4,5,6,7):[(1, 13)]},
            #     "expected_projectors": {(0,2,3,4): [(1,0), (-1, 14), (-1, 5), (1, 11)], (0,1): [(1, 0), (1, 1), (-1, 4), (-1, 5), (-1, 8), (-1, 9), (1, 12), (1, 13)], (1,2): [(1, 0), (-1, 5), (1, 2), (-1, 7), (-1, 8), (1, 13), (-1, 10), (1, 15)], (1,3): [(1, 0), (1, 1), (1, 2), (1, 3), (-1, 4), (-1, 5), (-1, 6), (-1, 7)], (1,4): [(1, 0), (-1, 9), (-1, 10), (1, 3), (1, 12), (-1, 5), (-1, 6), (1, 15)]}
            # }
        ]

    def test_compute_minimal_generating_sets(self):
        for case in self.test_cases:
            with self.subTest(case=case):
                case = copy.deepcopy(case)
                stab = Stabilizer(B=case["B"], n=case["n"], orbit_dictionary=case["orbit_dictionary"])
                stab.compute_minimal_generating_sets()
                for orbit_key, orbit in stab.orbit_dictionary.items():
                    with self.subTest(orbit_key=orbit_key):
                        self.assertEqual(orbit.Zs, case["expected_mgs"][orbit_key])
    
    def test_compute_projector_stabilizers(self):
        for case in self.test_cases:
            with self.subTest(case=case):
                case = copy.deepcopy(case)
                stab = Stabilizer(B=case["B"], n=case["n"], orbit_dictionary=case["orbit_dictionary"])
                for orbit_key, orbit in stab.orbit_dictionary.items():
                    with self.subTest(orbit_key=orbit_key):
                        orbit.Zs = case["expected_mgs"][orbit_key]
                
                stab.compute_projector_stabilizers()
                for orbit_key, orbit in stab.orbit_dictionary.items():
                    with self.subTest(orbit_key=orbit_key):
                        self.assertEqual(orbit.Zs, case["expected_projectors"][orbit_key])

    
    """
    def test_compute_projector_stabilizers(self):
        for case in self.test_cases:
            with self.subTest(case=case):
                case = copy.deepcopy(case)
                stab = Stabilizer(B=case["B"], n=case["n"])
                stab.orbits = case["expected_orbits"]
                stab.minimal_generating_sets = case["expected_mgs"]
                stab.compute_projector_stabilizers()
                self.assertEqual(stab.projectors, case["expected_projectors"])
    """
if __name__ == '__main__':
    unittest.main()
