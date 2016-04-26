#!/bin/bash
gfortran -o hoge.exe \
            variables.f03 \
	    arg_parse.f03 \
            apply_PBC.f03  \
            Convert_trj_PDB.f03 \
	    main.f03

#./aho.exe < para01.inp
