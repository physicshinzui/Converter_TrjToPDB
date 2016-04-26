FC=gfortran

OBJS= main.o

MODS=apply_pbc.mod convert_trj_pdb.mod

.SUFFIXES:.o .f03
.f03.o:
	$(FC) -c $< $(OBJS)
.SUFFIXES:.mod.f03
.f03.mod:
	$(FC) -c $< $(MODS)

$(TARGET):$(OBJS)
	$(FC) -o $@ $(OBJS)

