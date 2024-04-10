.DEFAULT:
	-touch $@
all: repose repose_new ramble wobble

repose.o: repose.f90 constants.o rng.o cell.o ddp.o opt.o 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c repose.f90
repose_new.o: repose_new.f90 constants.o rng.o cell.o ddp_new.o opt.o 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c repose_new.f90
ramble.o: ramble.f90 constants.o rng.o cell.o ddp.o md.o opt.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c ramble.f90
wobble.o: wobble.f90 constants.o rng.o cell.o ddp.o ld.o 
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c wobble.f90
cell.o: cell.f90 constants.o rng.o spglib_f08.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c cell.f90
ddp.o: ddp.f90 constants.o cell.o nn.o spglib_f08.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c ddp.f90
ddp_new.o: ddp_new.f90 constants.o cell.o nn.o spglib_f08.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c ddp_new.f90
opt.o: opt.f90 constants.o rng.o cell.o ddp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c opt.f90
md.o: md.f90 constants.o rng.o cell.o ddp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c md.f90
ld.o: ld.f90 constants.o rng.o cell.o ddp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c ld.f90
nn.o: $(NN)/src/nn.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c $(NN)/src/nn.f90
rng.o: rng.f90 constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c rng.f90
constants.o: constants.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c constants.f90
spglib_f08.o: spglib_f08.f90
	$(FC) $(FFLAGS) -c $(SPGLIB)/spglib-$(SPG_VERSION)/example/spglib_f08.f90
OBJ = nn.o opt.o md.o ld.o ddp.o cell.o constants.o rng.o spglib_f08.o ddp_new.o
clean: neat
	-rm -f repose ramble wobble repose_new
neat:
	-rm -f $(OBJ) repose.o repose_new.o ramble.o wobble.o *.mod
repose: $(OBJ) repose.o
	$(LD) $(OBJ) repose.o -o repose $(LDFLAGS)
repose_new: $(OBJ) repose_new.o
	$(LD) $(OBJ) repose_new.o -o repose_new $(LDFLAGS)
ramble: $(OBJ) ramble.o
	$(LD) $(OBJ) ramble.o -o ramble $(LDFLAGS)
wobble: $(OBJ) wobble.o
	$(LD) $(OBJ) wobble.o -o wobble $(LDFLAGS)
