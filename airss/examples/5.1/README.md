Example 5.1
===========

## Quantum Espresso free search for 2 atoms of Carbon at 100 GPa

#### Approximate runtime: 120 minutes (serial), 40 minutes (in parallel over 4 processes)

In this example, we will explore the energy landscape of Carbon at 100 GPa (or 1 MBar).
For comparison, the pressure at the centre of the Earth is about 350 GPa.

This example is a DFT analogue of Example 2.1, and is identical to Examples 3.1 and 4.1,
except that we use Quantum Espresso to relax the structures (rather than CASTEP or VASP).

## Prerequisites

To run this example, you will need to have the `pw.x` binary of Quantum Espresso installed
and on your path. This code is _not_ automatically installed as part of the AIRSS package, 
but can be downloaded from https://www.quantum-espresso.org/download. Ubuntu users could
try `(sudo) apt install quantum-espresso`.

You will also need a pseudopotential for carbon. This specific example has been set up to
work with the Quantum Espresso recommended pseudopotential from PSlibrary (https://dalcorso.github.io/pslibrary/).
Download it using:

    wget https://www.quantum-espresso.org/upf_files/C.pbe-n-rrkjus_psl.1.0.0.UPF

so that you have the file C.pbe-n-rrkjus_psl.1.0.0.UPF in the current directory.


## Input files

Examine the seedfile `C2.cell`:

    host:5.1 cjp10$ cat C2.cell 
    
    %BLOCK LATTICE_CART
    1.709975 0 0
    0 1.709975 0
    0 0 1.709975
    %ENDBLOCK LATTICE_CART
    
    #VARVOL=5
    
    %BLOCK POSITIONS_FRAC
    C 0.0 0.0 0.0 # C1 % NUM=1
    C 0.0 0.0 0.0 # C2 % NUM=1
    %ENDBLOCK POSITIONS_FRAC
    
    #MINSEP=1.3
    
    KPOINTS_MP_SPACING 0.07


In this file, the AIRSS directive `#VARVOL=5` specifies a volume of 5 cubic angstroms
for the generated cell, plus or minus 5 percent. There are two carbon atoms present
in the unit cell, labelled `C1` and `C2`. 

The AIRSS directive `#MINSEP=1.3` sets a minimum atom-atom separation of 1.3 angstroms.
This is not essential for the light elements, but for transition metals it is important 
to avoid core overlap to prevent poor convergence of the electronic structure.

Finally, the directive `KPOINTS_MP_SPACING` sets a k-point sampling density of 0.07, 
with 0.05 more appropriate for metallic (or nearly metallic systems).

The file `C2.qe` contains DFT/relaxation parameters used by Quantum Espresso. For information
on these, consult the Quantum Espresso documentation at https://www.quantum-espresso.org/Doc/INPUT_PW.html.


## Run the example

To run, use:

    host:5.1 cjp10$ airss.pl -press 100 -max 10 -seed C2 -qe

The command-line options `-press 100` sets the pressure to 100 GPa, `-max 10` tells AIRSS 
to generate 10 structures, and `-seed C2` gives the seedname of the .cell input file (`C2.cell`).
The last option `-qe` tells AIRSS to use Quantum Espresso for relaxation.

As a further option, you could add e.g. `-mpinp 4` to the above to execute Quantum Espresso in parallel,
provided that you have `mpirun` and your `pw.x` binary supports it. 

During the run xmgrace may be used to plot the "conv" files in order to monitor progress: `xmgrace *.conv`.
This will produce plots of enthalpy vs. iteration number for each structure as it relaxes.
Enthalpy is given in Rydbergs.


## Examine the output

After running, we can obtain an enthalpy-ranked list of structures using the cryan tool (`ca`).
Note that, due to the stochastic nature of AIRSS, your exact output will differ from
that shown here:

    host:5.1 cjp10$ ca -r
    C2-31019-5918-4        100.27     4.794    -165.716   2 C            Fd-3m      1
    C2-31019-5918-1         99.97     4.796       0.000   2 C            Fd-3m      1
    C2-31019-5918-7         99.95     4.796       0.000   2 C            Fd-3m      1
    C2-31019-5918-6        100.00     4.795       0.000   2 C            Fd-3m      1
    C2-31019-5918-9         99.98     4.795       0.000   2 C            Fd-3m      1
    C2-31019-5918-2         99.94     4.796       0.000   2 C            Fd-3m      1
    C2-31019-5918-5         99.86     4.796       0.000   2 C            Fd-3m      1
    C2-31019-5918-8        100.08     5.014       2.190   2 C            Pmma       1
    C2-31019-5918-3        100.72     4.616       2.484   2 C            P4/mmm     1


In this very small cell, most of the structures found are the diamond structure (space group Fd-3m). 
If the search is repeated at lower pressures, for example 1 GPa, more graphitic structures will be found.

