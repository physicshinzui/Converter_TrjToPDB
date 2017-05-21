#!/usr/local/bin/python

from MDAnalysis import Universe
import sys

"""
This program calculates center of mass and geometry.
At this time, I wanted to confirm if the com of s100b was canceled.

Caution: this program is specialized for s100b-CTD system.

Usage: python conform_com_cancel.py [ PDB file name ]   
"""



file_name = sys.argv[1]
print "Input file name : ", file_name

u = Universe(file_name)
f_out = open(file_name+"_comTraj.dat", "w")
print "No of snapshots: ", len(u.trajectory)

for i, ts in enumerate(u.trajectory):

    #Select the all atoms constitute s100b
    selected_atoms = u.select_atoms("resid 1-94")

    print "atom ids: ", selected_atoms.ids

    com = selected_atoms.center_of_mass()
    cog = selected_atoms.center_of_geometry()

    f_out.write(str(com[0]) + " " + str(com[1]) + " " + str(com[2]) + " \n")



