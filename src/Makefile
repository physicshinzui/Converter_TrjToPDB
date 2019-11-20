FC=gfortran
CFLAGS = -Wall -O2 -fbounds-check

PROGRAM = converter
OBJS    = variables.o \
          apply_pbc.o \
          convert_trj_pdb.o \
          main.o 

#.SUFFIXES: %.o %.f90
#*.o: *.f03
#	$(FC) -c $<

#***Making program
$(PROGRAM): $(OBJS)
	$(FC) -o $@ $^ $(CFLAGS) 

variables.o: variables.f03 
	$(FC) -c $<
apply_pbc.o: apply_PBC.f03 
	$(FC) -c $<
convert_trj_pdb.o: Convert_trj_pdb.f03 
	$(FC) -c $< 
main.o: main.f03 
	$(FC) -c $<
#
.PHONY:clean
clean:
	rm $(OBJS) $(PROGRAM) *.mod

