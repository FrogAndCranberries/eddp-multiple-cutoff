Example 3.1
===========

## Castep free search for 2 atoms of Carbon at 100GPa

In this example, we will explore the energy landscape of Carbon at 100 GPa (or 1MBar).
For comparison, the pressure at the centre of the Earth is about 350GPa.

This example is a DFT analogue of Example 2.1. We generate carbon structures, and use the 
CASTEP code to relax them.

**Note:** To run this example, you will need to have the `castep` binary installed and on your path. 
This code is _not_ automatically installed as part of the AIRSS package.
Refer to `http://www.castep.org/CASTEP/GettingCASTEP` for information on obtaining CASTEP.


## Input files

Examine the seedfile `C2.cell`:

    host:3.1 cjp10$ cat C2.cell 
    
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

Finally, the CASTEP directive `KPOINTS_MP_SPACING` sets a k-point sampling density of 0.07, 
with 0.05 more appropriate for metallic (or nearly metallic systems).

The file `C2.param` contains DFT/relaxation parameters used by CASTEP. For information
on these, consult the CASTEP documentation, or use `castep --help <keyword>`. 


## Run the example

To run, use:

    host:3.1 cjp10$ airss.pl -press 100 -max 10 -seed C2

The command-line options `-press 100` sets the pressure to 100 GPa, `-max 10` tells AIRSS 
to generate 10 structures, and `-seed C2` gives the seedname of the .cell input file (`C2.cell`).

As no other options are set, AIRSS will default to using the CASTEP code.

During the run xmgrace may be used to plot the "conv" files in order to monitor progress: `xmgrace *.conv`.
This will produce plots of enthalpy vs. iteration number for each structure as it relaxes.
Enthalpy is given in electronvolts (eV).

## Examine the output

After running, we can obtain an energy-ranked list of structures using the cryan tool (`ca`).
Note that, due to the stochastic nature of AIRSS, your exact output will differ from
that shown here:

    host:3.1 cjp10$ ca -r
    C2-90568-5971-4      100.00     4.699  -151.660  2 C            Fd-3m      1
    C2-90568-5971-5      100.00     4.699     0.000  2 C            Fd-3m      1
    C2-90568-5971-7      100.01     4.698     0.002  2 C            Fd-3m      1
    C2-90568-5971-8       99.99     4.698     0.002  2 C            Fd-3m      1
    C2-90568-5971-2       99.99     4.697     0.003  2 C            Fd-3m      1
    C2-90568-5971-1       99.94     4.698     0.003  2 C            Fd-3m      1
    C2-90568-5971-6      100.00     4.695     0.005  2 C            Fd-3m      1
    C2-90568-5971-10     100.00     4.691     0.010  2 C            Fd-3m      1
    C2-90568-5971-3      100.05     5.379     0.871  2 C            C2/m       1
    C2-90568-5971-9      100.03     4.477     2.555  2 C            P1         1

In this very small cell, most of the structures found are the diamond structure (space group Fd-3m). 
If the search is repeated at lower pressures, for example 1 GPa, more graphitic structures will be found.

