import unittest
import sys
import os
sys.path.insert(0, os.path.abspath("C:/Users/sanne/QAOA_LX-mixer"))

from Stabilizer import *

#TODO need a better test such as one where we get stabilizers that have more than 1 element... not sure if compute_projector_stabilizers work in that case
class TestStabilizer(unittest.TestCase):
    def __init__(self, methodname):
        super().__init__(methodname)
        self.B = [[0b1011, 0b1100, 0b0111, 0b0000, 0b1110, 0b1001, 0b0010, 0b0101]]
        self.n = 4
        
        self.stabilizer = Stabilizer(familiy_of_graphs=None, B=self.B, n=self.n)

    def test_check_if_orbit(self):
        self.stabilizer.check_if_orbit()
        self.assertEqual(self.stabilizer.orbits, [[11, 14, 7]])
    def test_compute_minimal_generating_sets(self): 
        self.stabilizer.check_if_orbit()
        self.stabilizer.compute_minimal_generating_sets()
        self.assertEqual(self.stabilizer.minimal_generating_sets, [[(1,13)]])
    def test_compute_projector_stabilizers(self):
        self.stabilizer.check_if_orbit()
        self.stabilizer.compute_minimal_generating_sets()
        self.stabilizer.compute_projector_stabilizers()
        self.assertEqual(self.stabilizer.projectors, [[(1, 0), (1, 13)]])



if __name__ == '__main__':
    unittest.main()
