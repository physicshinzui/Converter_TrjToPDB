FC=gfortran
CFLAGS = -Wall -O2 -fbounds-check

#***name of program and objects
PROGRAM = converter
OBJS    = variables.o \
          apply_pbc.o \
	  assign_prob.o \
          analysis.o \
          main.o 


#.SUFFIXES: %.o %.f90
#*.o: *.f03
#	$(FC) -c $<

#***Making program
$(PROGRAM): $(OBJS)
	$(FC) -o $@ $^ $(CFLAGS) 


variables.o: variables.f03 
	$(FC) -c $<
assign_prob.o: assign_prob.f03
	$(FC) -c $<
apply_pbc.o: apply_PBC.f03 
	$(FC) -c $<
analysis.o: analysis.f03
	$(FC) -c $< 
main.o: main.f03 
	$(FC) -c $<
#
.PHONY:clean
clean:
	rm $(OBJS) $(PROGRAM) *.mod

