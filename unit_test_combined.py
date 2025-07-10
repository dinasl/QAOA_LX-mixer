# --- import classes ---
from Stabilizer import Stabilizer
from Mixer import LXMixer
from Mixer import Orbit
 
B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]

lxmixer = LXMixer(B, 4)
lxmixer.compute_family_of_valid_graphs()
lxmixer.compute_all_orbits()

stabilizer = Stabilizer(B=lxmixer.B, n=4, orbit_dictionary=lxmixer.orbits)
stabilizer.compute_minimal_generating_sets()