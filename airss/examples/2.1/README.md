Example 2.1
===========

## Gulp free search for 2 atoms of carbon at 100 GPa, using a Tersoff potential

In this example, we will explore the energy landscape of Carbon at 100 GPa (or 1MBar).
For comparison, the pressure at the centre of the Earth is about 350GPa.

To proceed, we generate carbon structures, and use the Gulp code and a Tersoff potential
to relax them.


## Input files

First, examine the seedfile `C2.cell`, which reads as follows:

    host:2.1 cjp10$ cat C2.cell
    
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
    
    KPOINTS_MP_SPACING 0.1

In this file, the AIRSS directive `#VARVOL=5` specifies a volume of 5 cubic angstroms
for the generated cell, plus or minus 5 percent. There are two carbon atoms present
in the unit cell, labelled `C1` and `C2`. The AIRSS directive `#MINSEP=1.3` sets
a minimum atom-atom separation of 1.3 angstroms.

The file `C2.lib` gives the Tersoff potential parameters needed by Gulp.


## Run the example

To run, use:

    host:2.1 cjp10$ airss.pl -gulp -press 100 -max 10 -seed C2

The command-line options `-gulp` tell AIRSS to use the Gulp code, `-press 100` sets
the pressure to 100 GPa, `-max 10` tells AIRSS to generate 10 structures, and `-seed C2`
gives the seedname of the .cell input file (`C2.cell`).


## Examine the output

After running, we can obtain an energy-ranked list of structures using the cryan tool (`ca`).
Note that, due to the stochastic nature of AIRSS, your exact output will differ from
that shown here:

    host:2.1 cjp10$ ca -r
    C2-41554-5032-8      100.00     4.892      -4.214   2 C            C2/m       1
    C2-41554-5032-2      100.00     4.893       0.000   2 C            C2/m       1
    C2-41554-5032-3      100.00     4.892       0.000   2 C            C2/m       1
    C2-41554-5032-10     100.00     4.784       0.070   2 C            Fd-3m      1
    C2-41554-5032-1      100.00     4.784       0.070   2 C            Fd-3m      1
    C2-41554-5032-6      100.00     4.387       2.732   2 C            P-1        1
    C2-41554-5032-9      100.00     4.233       2.860   2 C            P-1        1
    C2-41554-5032-4      100.00     4.602       3.055   2 C            I4/mmm     1
    C2-41554-5032-7      100.00     4.643       3.070   2 C            I4/mmm     1
    C2-41554-5032-5      100.00     4.693       3.236   2 C            P-1        1
    host:2.1 cjp10$ 


## DFT versions

We revisit this same carbon example using density-functional theory in 3.1/4.1/5.1 
(using different DFT codes).


