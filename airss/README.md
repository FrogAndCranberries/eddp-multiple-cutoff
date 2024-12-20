            .o.       ooooo ooooooooo.    .oooooo..o  .oooooo..o
           .888.      '888' '888   'Y88. d8P'    'Y8 d8P'    'Y8
          .8:888.      888   888   .d88' Y88bo.      Y88bo.
         .8' '888.     888   888ooo88P'   ':Y8888o.   ':Y8888o.
        .88ooo8888.    888   888'88b.         ':Y88b      ':Y88b
       .8'     '888.   888   888  '88b.  oo     .d8P oo     .d8P
      o88o     o8888o o888o o888o  o888o 8::88888P'  8::88888P'

                Ab Initio Random Structure Searching
                Chris J. Pickard   (cjp20@cam.ac.uk)

                      Copyright (c) 2005-2022


Table of Contents

[TOC]

# Introduction

Ab initio Random Structure Searching (`AIRSS`) is a very simple, yet powerful and highly
parallel, approach to structure prediction. The concept was introduced in 2006 [1] and
its philosophy more extensively discussed in 2011 [2].

Random structures - or more precisely, random "sensible" structures - are generated
and then relaxed to nearby local energy minima. Particular success has been found
using density functional theory (DFT) for the energies, hence the focus on "ab initio"
random structure searching. The sensible random structures are constructed so that they
have reasonable densities, and atomic separations. Additionally they may embody
crystallographic, chemical or prior experimental/computational knowledge. Beyond these
explicit constraints the emphasis is on a broad, uniform, sampling of structure space.

`AIRSS` has been used in a number of landmark studies in structure prediction, from the
structure of SiH4 under pressure [1] to providing the theoretical structures which are used
to understand dense hydrogen (and anticipating the mixed Phase IV) [3], incommensurate
phases in aluminium under terapascal pressures [4], and ionic phases of ammonia [5].

The approach naturally extends to the prediction of clusters/molecules, defects in solids [6],
interfaces and surfaces [7].

**References**

    [1] C.J. Pickard and R.J. Needs, Phys. Rev. Lett., 97, 045504 (2006)
    [2] C.J. Pickard and R.J. Needs, J. Phys.: Condens. Matter 23, 053201 (2011)
    [3] C.J. Pickard and R.J. Needs, Nature Physics, 3, 473 (2007)
    [4] C.J. Pickard and R.J. Needs, Nature Materials, 9, 624 (2010)
    [5] C.J. Pickard and R.J. Needs, Nature Materials, 7, 775 (2008)
    [6] A.J. Morris, C.J. Pickard and R.J. Needs, Phys. Rev. B, 78, 184102 (2008)
    [7] G. Schusteritsch and C.J. Pickard, Phys. Rev. B, 90, 035424 (2014)

# Licence and citation

The AIRSS package is released under the `GPL 2.0` licence. See the `LICENCE` file for more details. You are not required to, but you might consider citing references [1] and [2] above in any work that makes use of the `AIRSS` package.

# Installation

Execute the following to perform a default installation:

    make ; make install ; make neat

The executables will be placed in `airss/bin`, which you should add to your path.

The installation can be tested:

    make check

[It is strongly recommended that `gcc` and `gfortran` version 5 and above are used to build the airss utilities. Other compiler families (such as ifort, clang) are not supported.]

# Documentation

The `AIRSS` package is documented through a growing list of worked
examples. Chapter 1 of the examples uses the provided pair potential
code (`pp3`). Please head to the `README.md` in `airss/examples`

# castep and AIRSS

The `AIRSS` package is tightly integrated with the `castep` first principles total
energy code. However, it is relatively straightforward to modify the scripts to
use alternative codes to obtain the core functionality. `xxx_relax` scripts for `castep`, `vasp`, `Quantum Espresso`, `pp3`, `gulp`, `psi4`, and `lammps` are provided and integrated 
with the `airss.pl` script.

# External utilities

The `AIRSS` package makes use of the following external utilities. Some of them are
available from package managers, others should be built locally.

## `spglib`

[This package is fetched and installed automatically]

An excellent library for finding and handling crystal symmetries,
written in `C` by Atsushi Togo.

http://atztogo.github.io/spglib/

## `cellsymm`

[This package is fetched and installed automatically]

This is a front-end to `spglib`, written by Michael Rutter.

http://www.tcm.phy.cam.ac.uk/sw/check2xsf/cellsym.tgz

[To be replaced by `c2x` in due course ...]

## `symmol`

[This package is fetched and installed automatically]

This fortran code symmetrises a group of atoms

    Pilati, T. & Forni, A. (2000). J. Appl. Cryst. 33, 417.and J. Appl. Cryst. (1998). 31, 503-504 doi:10.1107/S0021889898002180

https://www.mtg.msm.cam.ac.uk/files/symmol.zip

Before compilation a patch is applied to `symmol.f`

## `castep`

A high performance plane wave pseudopotential total energy
code. It is written and maintained by the Castep Developers Group, and available
under a no cost license for academic use.

http://www.castep.org/

The executable should be given the name `castep`, i.e. copy the default `castep.mpi`
or `castep.serial` to `castep` (which should be placed in your path).

## `OPTADOS`

Calculates high quality theoretical DOS, Projected-DOS, Joint-DOS, Optics
and core-loss spectroscopy.

http://www.tcm.phy.cam.ac.uk/~ajm255/optados/index.html

## `Gulp` [optional]

Structure prediction may be performed using a variety of empirical force fields,
as implemented in Julian Gale's powerful `Gulp` code.

http://gulp.curtin.edu.au/gulp/

## `LAMMPS` [optional]

http://lammps.sandia.gov/

## `pspot`

This is a directory containing the default `castep` `xx_00PBE.usp(cc)`
pseudopotentials. In more recent versions of `castep` the `QC5` set of high
throughput potentials provide an alternative. These can be used for general
searching, but tailored `OTFG` potentials are recommended for accurate results,
and/or very high pressures. It is assumed that pspot is in your home
directory. If not, set `PSPOT_DIR` appropriately. 

## `qhull`

Calculates the convex hulls.

http://www.qhull.org/

## `hull` [optional]

Ken Clarkson's convex hull code.

http://www.netlib.org/voronoi/hull.html

## `antiprism` [optional]

http://www.antiprism.com/files/antiprism-0.24.1.tar.gz

## `R/Rscript`

The statistical package `R` is used to visualise ternary convex hulls. The `ternary.r`
scripts can be executed using `Rscript`. The ggtern packaged is required.

https://cran.r-project.org/
http://www.ggtern.com/

## `openbabel`

Open Babel is a set of tools designed to convert chemical data. 
Used for the injection of `SMILES` strings.

## `xmgrace`

http://plasma-gate.weizmann.ac.il/Grace/

`Xmgrace/grace` is useful for visualising results. `AGR` scripts are generated to facilitate this.

## `cif2cell`

http://cif2cell.sourceforge.net/
http://www.sciencedirect.com/science/article/pii/S0010465511000336

This handy python utility can convert from `cif` files to a variety of electronic structure codes, including `castep`.

# The core AIRSS Scripts and Codes

## `airss.pl`

This perl script performs Ab Initio Random Structure Searching. See the set of
Examples for instruction in its use. Random structures are repeatedly generated
by buildcell, and relaxed with the chosen code (the default is `castep`).

## `buildcell`

This fortan code reads annotated `castep` cell files, and generates
random "sensible" structures from them. The type of randomness introduced
can be controlled through hash-tagged directives (which are treated as
comments and ignored by `castep`).

## `cabal`

A structure conversion tool.

    Usage: cabal in out < seed.in > seed.out
      in==out : Niggli reduce
      supports castep+,cell,shx,res,gulp*,cif*,psi4*,xtl,xyz(e)
      *output only +input only

The following converts a `castep` `.cell` file to a `SHLX` results file.

    cabal cell res < input.cell > output.res

## `ca`

A bash wrapper for `cryan` (see below). Uses the same command line options as `cryan`.

## `cryan`

A general purpose fortran program to analyse large amounts of structure data.
The structures are read from `STDIN`, for example:

     cat *.res | cryan -s
     gunzip -c lots.res.gz | cryan -f H2O -r
     find . -name "*.res" | xargs cat | cryan -m
     cat H2O-P21c.res | cryan -g

Experience suggests that cryan is suitable for the analysis of up to about 100K structures.
Other techniques are required for larger data sets.

## `castep_relax`

This bash script performs a self consistent geometry optimisation of the
specified structure using `castep`.

## `repose_relax`

This bash script performs structural optimisation using EDDPs and the
the `repose` code.

## `python_relax`

This bash script enables structural optimisation using a python script.

## `vasp_relax`

This bash script performs a self consistent geometry optimisation of the
specified structure using VASP.

## `qe_relax`

This bash script performs a self consistent geometry optimisation of the
specified structure using Quantum Espresso.

## `gulp_relax`

This bash script performs a geometry optimisation of the
specified structure using gulp.

## `pp3_relax`

This bash script performs a geometry optimisation of the
specified structure using `pp3`, a very simple pair potential code.

## `lammps_relax`

This bash script performs a geometry optimisation of the
specified structure using lammps.

[not currently recommended due to issues with structural optimisation]

## `psi4_relax`

This bash script performs a geometry optimisation of the
specified structure using `psi4`.

[not currently recommended due to issues with structural optimisation]

## `gencell`

This bash script will generate a set of recommended `castep` `.cell` and `.param` files
from a supplied unit cell volume, and atoms contained in the cell. It is strongly recommended
that this is used as a starting point for most projects.

## `pp3`

A simple pair potential code, for testing. It is used in the early chapters of 
the examples.

## `run.pl`

This perl script runs a batch of `castep` jobs in a directory. It is
useful for "polishing" your results, and high-throughput computation
in general. Failed runs are placed into `bad_castep`.

## `crud.pl`

The `castep` run daemon. A perl script for high-throughput batch calculations.
The structures to be relaxed are placed in the 'hopper' directory. Successful
calculations are placed in `good_castep`, and those that fail into `bad_castep`.
The script can be run as a daemon (i.e. continues running once the hopper is empty,
and waits for more).

## `spawn`

This script can be used to submit multiple jobs to the selection of
machines listed in the `~/.spawn` file. For example:

    node1 slots=8 root=
    node2 slots=8 root=
    node3 slots=12 root=

Typing `spawn airss.pl -seed Carbon` on your root node (on which it is
not advisable to run large jobs) will start a total of 28 instances of
`airss.pl` using the Carbon.* input files, on your 3 remote nodes. Spawn
uses ssh to run the commands remotely. Passwordless access to the
resources in your `.spawn` file is convenient 

The alternative to spawn or mpirun is to use the queueing system of a
multiuser computer cluster to submit multiple jobs. This should be
discussed with your system administrators.

## `spawn-slow`

Similar to the spawn script, this more slowly requests remote jobs. This
is recommended when launching the `run.pl` and `crud.pl` scripts.

## `despawn`

The spawn scripts records the PIDs of the remotely spawned jobs. This script
can be used to halt calculations in a controlled manner.

## `stopairss`

A script to kill spawned jobs. It will kill all jobs owned by you on
the remote nodes - so use with care. Use despawn in preference.

## `symm`

Finds the space group of the structure.

## `tidy.pl`

Removes the output of uncompleted calculations.

## `niggli`

Performs a Niggli transformation on all the SHLX res files in the current
directory.

## `prim`

Converts all the SHLX res files in the current directory to primitive cells.

## `conv`

Converts all the SHLX res files in the current directory to conventional cells.

# Structure Prediction


The use of the AIRSS scripts for a variety of structure prediction problems 
is illustrated by the collection of examples. Please now head to airss/examples.

