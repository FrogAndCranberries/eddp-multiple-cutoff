FC=gfortran
MATHLIB=
BUILD=fast
AIRSS=../../airss
NN=../../nn


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

all: ddp

ddp:
	(cd src/; make)

install:
	(mv src/forge bin/)
	(mv src/flock bin/)
	(mv src/frank bin/)
	(mv src/forge_new bin/)
	(mv src/flock_new bin/)
	(mv src/frank_new bin/)

neat:
	(cd src; make clean)

clean: neat
	(rm -f bin/forge)
	(rm -f bin/flock)
	(rm -f bin/frank)
	(rm -f bin/forge_new)
	(rm -f bin/flock_new)
	(rm -f bin/frank_new)