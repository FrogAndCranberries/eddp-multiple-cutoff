.DEFAULT:
	-touch $@
all: forge flock frank forge_new flock_new frank_new
nn.o: $(NN)/src/nn.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c $(NN)/src/nn.f90
forge.o: forge.f90 data.o constants.o nn.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c forge.f90
forge_new.o: forge_new.f90 data_new.o constants.o nn.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c forge_new.f90
flock.o: flock.f90 data.o constants.o nn.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c flock.f90
flock_new.o: flock_new.f90 data_new.o constants.o nn.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c flock_new.f90
frank.o: frank.f90 niggli_reduce.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c frank.f90
frank_new.o: frank_new.f90 niggli_reduce.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c frank_new.f90
niggli_reduce.o: niggli_reduce.f90 constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c niggli_reduce.f90
data.o: data.f90 constants.o nn.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c data.f90
data_new.o: data_new.f90 constants.o nn.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c data_new.f90
constants.o: constants.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c constants.f90
FORGE_OBJ = nn.o forge.o data.o constants.o 
FORGE_NEW_OBJ = nn.o forge_new.o data_new.o constants.o 
FLOCK_OBJ = nn.o flock.o data.o constants.o 
FLOCK_NEW_OBJ = nn.o flock_new.o data_new.o constants.o 
FRANK_OBJ = frank.o niggli_reduce.o constants.o 
FRANK_NEW_OBJ = frank_new.o niggli_reduce.o constants.o 
clean: neat
	-rm -f forge flock frank forge_new flock_new frank_new
neat:
	-rm -f $(FORGE_OBJ) $(FLOCK_OBJ) $(FRANK_OBJ) *.mod

forge: $(FORGE_OBJ) 
	$(LD) $(FORGE_OBJ) -o forge $(LDFLAGS)

forge_new: $(FORGE_NEW_OBJ) 
	$(LD) $(FORGE_NEW_OBJ) -o forge_new $(LDFLAGS)

flock: $(FLOCK_OBJ) 
	$(LD) $(FLOCK_OBJ) -o flock $(LDFLAGS)

flock_new: $(FLOCK_NEW_OBJ) 
	$(LD) $(FLOCK_NEW_OBJ) -o flock_new $(LDFLAGS)

frank: $(FRANK_OBJ) 
	$(LD) $(FRANK_OBJ) -o frank $(LDFLAGS)

frank_new: $(FRANK_NEW_OBJ) 
	$(LD) $(FRANK_NEW_OBJ) -o frank_new $(LDFLAGS)
