FC=gfortran   # GCC family only
FFLAGS=-Og -g # Optimise at -Og to give consistent results

#FFLAGS= -Og -g -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wno-implicit-interface -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace -fall-intrinsics -ffpe-trap=invalid,zero,overflow -fbounds-check -Wuninitialized -Wno-unused-dummy-argument -Wno-c-binding-type#-Werror -Wfatal-errors

CC=gcc # GCC family only
CFLAGS=

# Define the version of spglib to be used
export SPG_VERSION=1.14.1

PREFIX=$(PWD)

LDFLAGS= -L$(PREFIX)/lib -llapack -lblas -lsymspg
LD=$(FC)

GCC_VER_GTE46 := $(shell echo `$(FC) -dumpfullversion -dumpversion | cut -f1-2 -d.` \>= 4.6 | bc )
ifeq ($(GCC_VER_GTE46),0)
DFLAGS=-DCOMPAT
endif	

export

all: external internal 

internal: pp3 cabal buildcell cryan

ifeq ( , $(shell which wget))
	@echo '----------------------------------------------------------------------'
	@echo '  The "wget" command, used for fetching the external packages, should'
	@echo '  be installed to enable full symmetry functionality'
	@echo '----------------------------------------------------------------------'
endif
ifeq ( , $(shell which bc))
	@echo '----------------------------------------------------------------------'
	@echo '  The "bc" command, used for floating point shell arithmetic, should'
	@echo '  be installed to enable full functionality'
	@echo '----------------------------------------------------------------------'
endif

external: symmol spglib

spglib:
	(cd external/spglib; make)

symmol:
	(cd external/symmol; make)

pp3:
	(cd src/pp3/src; make)

cabal:
	(cd src/cabal/src; make)

buildcell:
	(cd src/buildcell/src; make)

cryan:
	(cd src/cryan/src; make)

install: install_external install_internal 

install_internal:
	(rm -f bin/pp3 && cp src/pp3/src/pp3 bin/pp3)
	(rm -f bin/cabal && cp src/cabal/src/cabal bin/cabal)
	(rm -f bin/buildcell && cp src/buildcell/src/buildcell bin/buildcell)
	(rm -f bin/cryan && cp src/cryan/src/cryan bin/cryan)
	@echo
	@echo 'Add '$(PWD)'/bin to your path by placing this line in your ~/.bashrc:'
	@echo 'export PATH='$(PWD)'/bin:$$PATH'
	@echo 'To update your path "source ~/.bashrc"'
	@([ -d .git ] && (printf "echo '";head -1 VERSION.md | tr -d '\n';printf " commit "; git describe --abbrev=18 --always  --long --match v* master | tr -d '\n';printf "'") > bin/airss_version) || echo " "
	@chmod u+x bin/airss_version

install_external:
	(rm -f bin/symmol && cp external/symmol/symmol bin/symmol)

check:
	(bash bin/check_airss)

neat_internal:
	(cd src/pp3/src; make clean)
	(cd src/cabal/src; make clean)
	(cd src/buildcell/src; make clean)
	(cd src/cryan/src; make clean)

neat_external:
	(cd external/symmol; make clean)
	(cd external/spglib; make clean)
	(cd lib && rm -f *.a *.so)

neat: neat_internal neat_external

clean: neat
	(rm -f bin/pp3 bin/cabal bin/cryan bin/buildcell bin/symmol)

dist: clean
	rm -fr .check external/spglib/master.zip external/symmol/symmol.zip
	rm -fr include/* lib/*
	find examples/ -name "*-*.*" | xargs rm	-f
	find examples/ -name ".spawnpid*" | xargs rm -f
	#tar -czf ../airss-`date "+%d%m%Y"`.tgz  --exclude=".*" -C .. airss
	@echo 'For release use "hg archive -t tgz ../airss-x.y.z.tgz"'


