In this example we will use random searching to find the ground state of a Lennard-Jones solid.

    host:1.1 cjp10$ ls
    Al.cell	Al.pp	README
    host:1.1 cjp10$ cat Al.cell 
    %BLOCK LATTICE_CART
    2 0 0
    0 2 0
    0 0 2
    %ENDBLOCK LATTICE_CART.   
 
    %BLOCK POSITIONS_FRAC
    Al 0.0 0.0 0.0 # Al1 % NUM=8
    %ENDBLOCK POSITIONS_FRAC
 
    #MINSEP=1.5
    host:1.1 cjp10$

Al.cell is the "seed" input file which describes how the random structures should be built (by buildcell). In this case, 8 atoms are placed randomly into a unit cell of random shape, and a volume within 5% of 8 Ang^3/atom. If a larger number of atoms were requested then the volume of the cell would be scaled appropriately. Any configurations in which atoms are closer than 1.5 Ang are rejected. 

We now perform the search:

    host:1.1 cjp10$ airss.pl -pp3 -max 20 -seed Al
    host:1.1 cjp10$ 

The AIRSS perl script is executed. The pp3 pair potential code is chosen for the structural optimisations (the default is CASTEP), and 20 relaxations are requested (the default being a very large number - which would require manual halting of the scripts by the user). The Al.pp file contains the parameters for the pair potential. The results can now be analysed.

    host:1.1 cjp10$ ca -r
    Al-91855-9500-1       -0.00     7.561    -6.659  8 Al           P63/mmc    1
    Al-91855-9500-6       -0.00     7.561     0.000  8 Al           P63/mmc    1
    Al-91855-9500-5       -0.00     7.561     0.000  8 Al           P63/mmc    1
    Al-91855-9500-10       0.00     7.564     0.005  8 Al           P-1        1
    Al-91855-9500-17      -0.00     7.564     0.005  8 Al           P-1        1
    Al-91855-9500-3       -0.00     7.564     0.005  8 Al           C2/m       1
    Al-91855-9500-18      -0.00     7.564     0.005  8 Al           Pmmm       1
    Al-91855-9500-20      -0.00     7.564     0.005  8 Al           C2/m       1
    Al-91855-9500-4       -0.00     7.564     0.005  8 Al           Pmmm       1
    Al-91855-9500-16      -0.00     7.564     0.005  8 Al           C2/m       1
    Al-91855-9500-2       -0.00     7.564     0.005  8 Al           Cmmm       1
    Al-91855-9500-11      -0.00     7.564     0.005  8 Al           C2/m       1
    Al-91855-9500-9       -0.00     7.564     0.005  8 Al           Cmmm       1
    Al-91855-9500-8        0.00     7.564     0.005  8 Al           Fmmm       1
    Al-91855-9500-19       0.00     7.784     0.260  8 Al           C2/m       1
    Al-91855-9500-12      -0.00     8.119     0.700  8 Al           R-3c       1
    Al-91855-9500-14      -0.00     8.446     0.705  8 Al           Cm         1
    Al-91855-9500-13       0.00     8.453     0.794  8 Al           C2/m       1
    Al-91855-9500-15       0.00     8.465     0.802  8 Al           C2/m       1
    Al-91855-9500-7        0.00     8.505     0.834  8 Al           P21/m      1
    host:1.1 cjp10$

The above is a "ranking" of the results according to energy (strictly, enthalpy). The first column is the unique structure name, the second - pressure, the third - volume per formula unit (fu), the fourth - enthalpy/fu for the first, and then relative to that, the fifth - number of fu, the sixth - chemical formula, the seventh - symmetry, and the eighth - number of repeats.

For such a simple system you are likely to find the HCP ground state (space group P63/mmc) within these 20 attempts. Of course, this is a stochastic approach, and you may be unlucky. If so - simply run the AIRSS script again.

In the above example, we found the HCP structure 3 times. It can be helpful to unify repeated structures.

    host:1.1 cjp10$ ca -u 0.01 -r
    Al-91855-9500-1       -0.00     7.561    -6.659  8 Al           P63/mmc    3
    Al-91855-9500-10       0.00     7.564     0.005  8 Al           P-1       11
    Al-91855-9500-19       0.00     7.784     0.260  8 Al           C2/m       1
    Al-91855-9500-12      -0.00     8.119     0.700  8 Al           R-3c       1
    Al-91855-9500-14      -0.00     8.446     0.705  8 Al           Cm         1
    Al-91855-9500-13       0.00     8.453     0.794  8 Al           C2/m       1
    Al-91855-9500-15       0.00     8.465     0.802  8 Al           C2/m       1
    Al-91855-9500-7        0.00     8.505     0.834  8 Al           P21/m      1
    host:1.1 cjp10$

The value "0.01" is the threshold of similarity. If you are happy with unification, you can make it permanent.

    host:1.1 cjp10$ ca -u 0.01 -r --delete
    Deleting files. To confirm type <Enter>
    Al-91855-9500-1       -0.00     7.561    -6.659  8 Al           P63/mmc    3
    Al-91855-9500-10       0.00     7.564     0.005  8 Al           P-1       11
    Al-91855-9500-19       0.00     7.784     0.260  8 Al           C2/m       1
    Al-91855-9500-12      -0.00     8.119     0.700  8 Al           R-3c       1
    Al-91855-9500-14      -0.00     8.446     0.705  8 Al           Cm         1
    Al-91855-9500-13       0.00     8.453     0.794  8 Al           C2/m       1
    Al-91855-9500-15       0.00     8.465     0.802  8 Al           C2/m       1
    Al-91855-9500-7        0.00     8.505     0.834  8 Al           P21/m      1
    host:1.1 cjp10$ ca -r
    Al-91855-9500-1       -0.00     7.561    -6.659  8 Al           P63/mmc    3
    Al-91855-9500-10       0.00     7.564     0.005  8 Al           P-1       11
    Al-91855-9500-19       0.00     7.784     0.260  8 Al           C2/m       1
    Al-91855-9500-12      -0.00     8.119     0.700  8 Al           R-3c       1
    Al-91855-9500-14      -0.00     8.446     0.705  8 Al           Cm         1
    Al-91855-9500-13       0.00     8.453     0.794  8 Al           C2/m       1
    Al-91855-9500-15       0.00     8.465     0.802  8 Al           C2/m       1
    Al-91855-9500-7        0.00     8.505     0.834  8 Al           P21/m      1
    host:1.1 cjp10$

