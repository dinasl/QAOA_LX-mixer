import networkx as nx
import numpy as np
from dataclasses import dataclass, field
from itertools import combinations, product

from utils import *
from Stabilizer import *
from Mixer import *

def test_compute_family_of_graphs():
    B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]
    LX = LXMixer(B)
    family_of_valid_graphs = LX.compute_family_of_graphs()
    print("Computed family of valid graphs:\n")
    for graph in family_of_valid_graphs:
        X = graph.X
        print(pauli_int_to_str(X, "X") + " ")
    print("\nCorrect family of valid graphs:\n",
          " IIXI, IXXX, XIXI, XXIX, IXIX, XIII, XXXX")