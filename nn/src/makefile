.DEFAULT:
	-touch $@
all: nn_sine nn_fourier nn_poly
nn_sine.o: nn_sine.f90 nn.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c nn_sine.f90
nn_fourier.o: nn_fourier.f90 nn.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c nn_fourier.f90
nn_poly.o: nn_poly.f90 nn.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c nn_poly.f90
nn.o: nn.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c nn.f90
constants.o: constants.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c constants.f90
SRC = nn_fourier.f90 nn_sine.f90 nn_poly.f90 nn.f90 constants.f90
OBJ = nn.o constants.o
clean: neat
	-rm -f nn_test
neat:
	-rm -f $(OBJ) nn_sine.o nn_fourier.o nn_poly.o *.mod
nn_sine: $(OBJ) nn_sine.o
	$(LD) $(OBJ) nn_sine.o -o nn_sine $(LDFLAGS)
nn_fourier: $(OBJ) nn_fourier.o
	$(LD) $(OBJ) nn_fourier.o -o nn_fourier $(LDFLAGS)
nn_poly: $(OBJ) nn_poly.o
	$(LD) $(OBJ) nn_poly.o -o nn_poly $(LDFLAGS)
