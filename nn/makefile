FC=gfortran
MATHLIB=
BUILD=

ifeq ($(BUILD), fast)
	ifeq ($(FC), ifort)
		OPT= -O3
		DEBUG=
	else
		OPT= -Ofast -fexternal-blas
		DEBUG=
	endif
else ifeq ($(BUILD), debug)
	ifeq ($(FC), ifort)
		OPT= -O0 -g
		DEBUG= -traceback -check all -ftrapuv -debug all
	else
		OPT= -O0 -g
		DEBUG= -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising \
			-Waliasing -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic \
			-fbacktrace -fall-intrinsics -ffpe-trap=invalid,zero,overflow -fbounds-check \
			-Wuninitialized -Wfatal-errors
	endif
else
	OPT= -O2
	DEBUG=
endif

ifeq ($(shell uname -s), Darwin)
    ifeq ($(FC), ifort)
	FFLAGS= $(OPT) $(DEBUG) -I"${MKLROOT}/include" -heap-arrays -fpp -D IFORT -qopenmp
	LDFLAGS= ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a \
		${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl
    else
	ifeq ($(MATHLIB), mkl)
		FFLAGS= $(OPT) $(DEBUG) -m64  -I"${MKLROOT}/include" -fopenmp -cpp
		LDFLAGS= -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_rt -lpthread -lm \
		-ldl -fopenmp
	else
		FFLAGS= $(OPT) $(DEBUG) -fopenmp -cpp
		LDFLAGS= -llapack -lblas -fopenmp
	endif
    endif
else
    ifeq ($(FC), ifort)
	FFLAGS= $(OPT) $(DEBUG) -I"${MKLROOT}/include" -heap-arrays -fpp -D IFORT -qopenmp
	LDFLAGS= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a \
		${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl
    else
	ifeq ($(MATHLIB), mkl)
		FFLAGS= $(OPT) $(DEBUG) -m64  -I"${MKLROOT}/include" -fopenmp -cpp
		LDFLAGS= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a \
			${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a \
			-Wl,--end-group -lgomp -lpthread -lm -ldl
	else
		FFLAGS= $(OPT) $(DEBUG) -fopenmp -cpp
		LDFLAGS= -llapack -lblas -fopenmp
	endif
    endif
endif

LD=$(FC)
PREFIX=$(PWD)
LDFLAGS+= -DSINGLE
export

all: nn 

nn:
	(cd src/; make)

install:
	(mv src/nn_sine bin/)
	(mv src/nn_fourier bin/)
	(mv src/nn_poly bin/)

test:
	(cd examples && time ../bin/nn_sine)
	(cd examples && time ../bin/nn_fourier)
	(cd examples && time ../bin/nn_poly)

neat:
	(cd src; make clean)

clean: neat
	(rm -f bin/nn_sine)
	(rm -f bin/nn_fourier)
	(rm -f bin/nn_poly)

dist: clean
	tar -czf ../bp-`date "+%d%m%Y"`.tgz  --exclude=".*" -C .. nn
	echo 'For release use "hg archive"'
