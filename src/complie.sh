#!/bin/bash
gfortran -o hoge.exe \
            apply_PBC.f03  \
            Convert_trj_PDB.f03 \
	    main.f03

#./aho.exe < para01.inp
