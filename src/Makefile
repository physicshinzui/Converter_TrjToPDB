FC=gfortran
#CFLAGS = -Wall -O2

OBJS = variables.o \
       apply_pbc.o \
       convert_trj_pdb.o \
       main.o 

EXEC = converter

all: $(EXEC) 

variables.o:  variables.f03
	$(FC) -c variables.f03
apply_pbc.o:  apply_pbc.f03
	$(FC) -c apply_pbc.f03
convert_trj_pdb.o: convert_trj_pdb.f03
	$(FC) -c convert_trj_pdb.f03
main.o: main.f03
	$(FC) -c main.f03 

#.f03.o:
#	$(FC) -c $< $(OBJS)

$(EXEC):$(OBJS)
	$(FC) -o $@ $(OBJS)

.PHONY:clean
clean:
	rm *.o *.mod

