FC=gfortran
MATHLIB=
BUILD=
AIRSS=../../airss
NN=../../nn
SPGLIB=$(AIRSS)/external/spglib
SPG_VERSION=1.14.1

ifeq ($(BUILD), fast)
	ifeq ($(FC), ifort)
		OPT= -O3
		DEBUG=
	else
		OPT= -O3 -fexternal-blas
		DEBUG=
	endif
else ifeq ($(BUILD), debug)
	ifeq ($(FC), ifort)
		OPT= -O0 -g
		DEBUG= -traceback -check all -ftrapuv -debug all
	else
		OPT= -Og -g3
		DEBUG= -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising \
			-Waliasing -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic \
			-fbacktrace -fall-intrinsics -ffpe-trap=invalid,zero,overflow -fbounds-check \
			-Wuninitialized -Wfatal-errors
	endif
else
	OPT= -O2 -g
	DEBUG=
endif

ifeq ($(shell uname -s), Darwin)
    ifeq ($(FC), ifort)
	FFLAGS= $(OPT) $(DEBUG) -I"${MKLROOT}/include" -heap-arrays -fpp -D IFORT -qopenmp
	LDFLAGS= ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a \
		${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl  -L$(AIRSS)/lib -lsymspg
    else
	ifeq ($(MATHLIB), mkl)
		FFLAGS= $(OPT) $(DEBUG) -m64  -I"${MKLROOT}/include" -fopenmp -cpp
		LDFLAGS= -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_rt -lpthread -lm \
		-ldl -fopenmp -L$(AIRSS)/lib -lsymspg
	else
		FFLAGS= $(OPT) $(DEBUG) -cpp -fopenmp
		LDFLAGS= -llapack -lblas -L$(AIRSS)/lib -lsymspg -fopenmp
	endif
    endif
else
    ifeq ($(FC), ifort)
	FFLAGS= $(OPT) $(DEBUG) -I"${MKLROOT}/include" -heap-arrays -fpp -D IFORT -qopenmp
	LDFLAGS= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a \
		${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl -L$(AIRSS)/lib -lsymspg
    else
	ifeq ($(MATHLIB), mkl)
		FFLAGS= $(OPT) $(DEBUG) -m64  -I"${MKLROOT}/include" -fopenmp -cpp
		LDFLAGS= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a \
			${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a \
			-Wl,--end-group -lgomp -lpthread -lm -ldl -L$(AIRSS)/lib -lsymspg
	else
		FFLAGS= $(OPT) $(DEBUG) -cpp -fopenmp
		LDFLAGS= -llapack -lblas -L$(AIRSS)/lib -lsymspg -fopenmp
	endif
    endif
endif


LD=$(FC)
PREFIX=$(PWD)

export

all: repose

repose:
	(cd src/; make)

install:
	(mv src/repose bin/)
	(mv src/repose_new bin/)
	(mv src/ramble bin/)
	(mv src/wobble bin/)

bench:
	(cd benchmarks/SiO2-1536-ramble && ./bench.sh 1 `getconf _NPROCESSORS_ONLN`)

neat:
	(cd src; make clean)

clean: neat
	(rm -f bin/repose)
	(rm -f bin/repose_new)
	(rm -f bin/ramble)
	(rm -f bin/wobble)
