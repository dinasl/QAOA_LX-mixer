# --- import classes ---
from Stabilizer import Stabilizer
from Mixer import LXMixer

#B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]
B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011, 0b0000, 0b1111, 0b1011, 0b1101, 0b0110, 0b0010, 0b0101, 0b1000, 0b0001, 0b0111]
lxmixer = LXMixer(B, 4)
lxmixer.compute_family_of_valid_graphs()
lxmixer.compute_all_orbits()

stabilizer = Stabilizer(B=lxmixer.B, n=4, orbit_dictionary=lxmixer.orbits)
stabilizer.compute_minimal_generating_sets()
stabilizer.compute_projector_stabilizers()