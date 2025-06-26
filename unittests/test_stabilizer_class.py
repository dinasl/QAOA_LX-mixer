import unittest

class TestStabilizer(unittest.TestCase):
    def __init__(self, methodname):
        super().__init__(methodname)

        self.B = [0b1011, 0b1100, 0b0111, 0b0000, 0b1110, 0b1001, 0b0010, 0b0101]
        self.n = 4
    
    def test_check_if_orbit(self):
        pass
    def test_compute_minimal_generating_set(self):
        pass
    def test_compute_projector_stabilzers(self):
        pass

