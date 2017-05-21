from MDAnalysis import Universe
from MDAnalysis.analysis.align import *
from MDAnalysis.analysis.rms import rmsd

ref  = Universe("npt_1.pdb")
traj = Universe("allframes.pdb")

rms_fit_trj(traj, ref, select=("resid 95-108 and name CA","resid 95-108 and name CA"),filename=None,rmsdfile="rmsd.dat")
