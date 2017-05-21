#!/usr/local/bin/python
import argparse
import time
import subprocess

p = argparse.ArgumentParser()
p.add_argument("-l","--trajectory_file_list", required=True)
args = p.parse_args()

trj_name_list = args.trajectory_file_list

print trj_name_list
fin_list = open(trj_name_list, "r")

for i, iname in enumerate(fin_list, 1):
   print i, iname 
   fin_para = open("base_para.inp","r")
   name_para_inp = "para"+str(i)+".inp"
   fout     = open(name_para_inp, "w")
   for line in fin_para:
       line = line.replace("TRAJ_DATA",iname.rstrip())
       line = line.replace("TRAJ_OUT_NAME","frame"+str(i)+".pdb")
       fout.write(line)
   fout.close()
