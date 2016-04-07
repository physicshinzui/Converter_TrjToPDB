import MDAnalysis
import sys
import numpy as np
from MDAnalysis.lib.distances import apply_PBC

if __name__ == "__main__":

  universe = MDAnalysis.Universe("trj.pdb") 
  protein  = universe.select_atoms("protein")
#  print protein.n_atoms
  with MDAnalysis.Writer("protein.pdb", protein.n_atoms) as W:
      for ts in universe.trajectory:
          W.write(protein)
