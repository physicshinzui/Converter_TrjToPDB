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

for i, iname in enumerate(fin_list):
   print i, iname 
   fin_para = open("base_para.inp","r")
   name_para_inp = "para"+str(i)+".inp"
   fout     = open(name_para_inp, "w")
   for line in fin_para:
       line = line.replace("TRAJ_DATA",iname.rstrip())
       line = line.replace("TRAJ_OUT_NAME","frame"+str(i)+".pdb")
       fout.write(line)
   fout.close()

##execution of fortran
#   fexe_in = open("base_exe.sh","r")
#   fexe_out= open("exe"+str(i)+".sh","w")
#   for line in fexe_in:
#       line = line.replace("PARA_INP", name_para_inp)
#       fexe_out.write(line)
#
#   cmd = "bash exe"+str(i)+".sh &"
#   print cmd
#   subprocess.call(cmd, shell=True)
#   time.sleep(0.1)

