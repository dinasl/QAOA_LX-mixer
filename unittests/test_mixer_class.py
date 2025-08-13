import unittest
from Mixer import LXMixer

class TestLXMixer(unittest.TestCase):
    def test_compute_family_of_graphs(self):
        correct_family_of_valid_graphs = {
            0b0010: [(0, 1)],
            0b0101: [(1, 2)],
            0b0111: [(0, 2), (3, 4)],
            0b1000: [(1, 3)],
            0b1010: [(0, 3), (2, 4)],
            0b1101: [(0, 4), (2, 3)],
            0b1111: [(1, 4)]
        }
        
        B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]
        lxmixer = LXMixer(B, 4)
        lxmixer.compute_family_of_valid_graphs()
        
        # Sort edges for comparison
        for key in correct_family_of_valid_graphs:
            correct_family_of_valid_graphs[key] = sorted(correct_family_of_valid_graphs[key])
        
        for key in lxmixer.family_of_valid_graphs:
            lxmixer.family_of_valid_graphs[key] = sorted(lxmixer.family_of_valid_graphs[key])
        
        self.assertEqual(lxmixer.family_of_valid_graphs, correct_family_of_valid_graphs)
    
    def test_compute_all_orbits(self):
        correct_orbits = {
            (0, 2, 3, 4) : [0b1101, 0b1010],
            (0, 1) : [0b0010],
            (1, 2) : [0b0101],
            (1, 3) : [0b1000],
            (1, 4) : [0b1111]
        }
        
        B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]
        lxmixer = LXMixer(B, 4)
        lxmixer.family_of_valid_graphs = {
            0b0010: [(0, 1)],
            0b0101: [(1, 2)],
            0b0111: [(0, 2), (3, 4)],
            0b1000: [(1, 3)],
            0b1010: [(0, 3), (2, 4)],
            0b1101: [(0, 4), (2, 3)],
            0b1111: [(1, 4)]
        }
        lxmixer.compute_all_orbits()
        
        orbits = {}
        for nodes, orbit in lxmixer.orbits.items():
            orbits[nodes] = orbit.Xs
        
        self.assertEqual(orbits, correct_orbits)

if __name__ == '__main__':
    unittest.main()