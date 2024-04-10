Example 5.2
===========

## Quantum Espresso free search for 8 atoms of Hydrogen at 100 GPa, followed by molecular units

#### Approximate runtime: 60 minutes (free search) and 60 minutes (units search), parallelised over 4 processes

In this example we explore hydrogen at 100 GPa. This example is identical to Example 3.2, except that we
use Quantum Espresso for structural relaxations.

In the first part, we carry out a free search. The molecular units that emerge in low-enthalpy structures from
the free search are then extracted, and used to carry out a follow-up 'units-based' search.


## Prerequisites

As with Example 5.1, you will need to have the `pw.x` binary of Quantum Espresso installed
and on your path.

You will also need a pseudopotential for hydrogen. This specific example has been set up to
work with the Quantum Espresso recommended pseudopotential from PSlibrary (https://dalcorso.github.io/pslibrary/).
Download it using:

    wget https://www.quantum-espresso.org/upf_files/H.pbe-rrkjus_psl.1.0.0.UPF

so that you have the file H.pbe-rrkjus_psl.1.0.0.UPF in the current directory.


## Input files

Examine `H.cell` first:

    host:5.2 cjp10$ cat H.cell

    #VARVOL=2.5
    #SPECIES=H
    #NATOM=8
    
    #SLACK=0.25
    #OVERLAP=0.1
    #MINSEP=1 H-H=0.7
    #COMPACT
    
    KPOINTS_MP_SPACING 0.07


The AIRSS directive `#VARVOL=2.5` specifies a volume of 2.5 cubic angstroms
per H atom, plus or minus 5 percent. We request hydrogen (`#SPECIES=H`),
and 8 atoms in the unit cell (`#NATOM=8`).

`#MINSEP=0.7` sets a minimum H-H separation of 0.7 angstroms.
Remember to avoid core overlap to prevent poor convergence of the electronic structure.

Finally, the directive `KPOINTS_MP_SPACING` sets a k-point sampling density of 0.07,
with 0.05 more appropriate for metallic (or nearly metallic systems).

The file `H.qe` contains DFT/relaxation parameters used by Quantum Espresso. For information
on these, consult the Quantum Espresso documentation at https://www.quantum-espresso.org/Doc/INPUT_PW.html.


## Run the example

We request 10 structures (`-max 10`) to be relaxed at 100 GPa (`-press 100`) using Quantum Espresso (`-qe`).
You could also add e.g. `-mpinp 4` to run Quantum Espresso in parallel.

    host:5.2 cjp10$ airss.pl -press 100 -max 10 -seed H -qe


## Examine the output

Use cryan (`ca`) to obtain an enthalpy-ranked list. Your exact output will differ, due to the stochastic
nature of AIRSS:

    host:5.2 cjp10$ ca -r
    H-10203-8632-9         100.01     2.301     -13.694   8 H            Pca21      1
    H-10203-8632-1         100.00     2.302       0.000   8 H            Pca21      1
    H-10203-8632-4          99.99     2.299       0.001   8 H            P21/c      1
    H-10203-8632-2         100.00     2.298       0.001   8 H            P21        1
    H-10203-8632-6         100.02     2.298       0.001   8 H            P212121    1
    H-10203-8632-8         100.00     2.298       0.002   8 H            P212121    1
    H-10203-8632-3         100.03     2.302       0.003   8 H            C2         1
    H-10203-8632-5         100.03     2.302       0.003   8 H            C2         1
    H-10203-8632-10        100.03     2.303       0.003   8 H            C2         1
    H-10203-8632-7         100.02     2.320       0.014   8 H            P1         1


In this case, the lowest-enthalpy structure found was H-10203-8632-9, whose structure can be found
in SHX format as H-10203-8632-9.res.


## Extract molecular units

Using the cryan code, we can directly extract the molecular units present in our lowest-enthalpy structure:

    host:5.2 cjp10$ cryan -g -bs 1 < H-10203-8632-9.res > H2.cell
    
    Degrees of freedom:    14  0.737
    
    Number of units:        4
    
    host:5.2 cjp10$ cat H2.cell
    %BLOCK LATTICE_CART
       1.89131   0.00000   0.00000
      -0.00120   2.99187   0.00000
      -1.88967  -2.99300   3.25316
    %ENDBLOCK LATTICE_CART
    
    #TARGVOL=4.60
    
    %BLOCK POSITIONS_ABS
        H  -0.90374  -1.98969   2.69356 # 1-D(inf)h % NUM=1
        H  -1.35056  -1.63007   3.14979 # 1-D(inf)h
    %ENDBLOCK POSITIONS_ABS
    
    #SYMMOPS=1
    ##SGRANK=20
    #NFORM=4
    #SLACK=0.25
    #OVERLAP=0.1
    #COMPACT
    #CELLADAPT
    #MINSEP=1.0 H-H=1.55


## Units-based search

We can use `H2.cell` to carry out a units-based search. Add the line `KPOINTS_MP_SPACING 0.07` from H.cell to your 
new H2.cell, and copy `H.qe` to `H2.qe`.

    host:5.2 cjp10$ echo "KPOINTS_MP_SPACING 0.07" >> H2.cell
    host:5.2 cjp10$ cp H.qe H2.qe

Run the search as above, but specify `H2` as the seed, rather than H.

    host:5.2 cjp10$ airss.pl -press 100 -max 10 -seed H2 -qe

Obtain an enthalpy-ranked list. Entries prefixed with `H` are from our first search, and with `H2` from our second
units-based search.

    host:5.2 cjp10$ ca -r
    H-10203-8632-9         100.01     2.301     -13.694   8 H            Pca21      1
    H2-21648-3658-5        100.00     2.302       0.000   8 H            Pca21      1
    H2-21648-3658-2        100.02     2.302       0.000   8 H            Pca21      1
    H-10203-8632-1         100.00     2.302       0.000   8 H            Pca21      1
    H2-21648-3658-8        100.00     2.299       0.001   8 H            P21/c      1
    H-10203-8632-4          99.99     2.299       0.001   8 H            P21/c      1
    H-10203-8632-2         100.00     2.298       0.001   8 H            P21        1
    H-10203-8632-6         100.02     2.298       0.001   8 H            P212121    1
    H-10203-8632-8         100.00     2.298       0.002   8 H            P212121    1
    H-10203-8632-3         100.03     2.302       0.003   8 H            C2         1
    H-10203-8632-5         100.03     2.302       0.003   8 H            C2         1
    H-10203-8632-10        100.03     2.303       0.003   8 H            C2         1
    H-10203-8632-7         100.02     2.320       0.014   8 H            P1         1
    H2-21648-3658-9         99.99     2.320       0.014   8 H            C2/m       1
    H2-21648-3658-1        100.00     2.333       0.015   8 H            Cmce       1
    H2-21648-3658-4        100.01     2.333       0.015   8 H            Cmce       1
    H2-21648-3658-10        99.99     2.337       0.020   8 H            C2/c       1
    H2-21648-3658-7        100.01     2.344       0.028   8 H            I4/mmm     1
    H2-21648-3658-6        100.03     2.261       0.035   8 H            C2/m       1

