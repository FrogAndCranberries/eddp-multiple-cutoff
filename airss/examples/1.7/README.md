In this example we build a slab supercell of Lennard-Jones HCP crystal, and add another atom to the surface randomly.

host:1.7 cjp10$ cat Al.cell 
%BLOCK LATTICE_CART
       1.908194408661965      -1.101696555507125       0.000000000000000
      -0.000000000000000       2.203393111014249       0.000000000000000
       0.000000000000000       0.000000000000000      35.964580000000005
#FIX
%ENDBLOCK LATTICE_CART

%BLOCK POSITIONS_FRAC
 Al   0.3333333333333333  -0.3333333333333333  0.3499999999999998 # Al1 
 Al   0.3333333333333333  -0.3333333333333333  0.4499999999999998 # Al2
 Al   0.0000000000000000   0.0000000000000000  0.2999999999999998 # Al3 
 Al   0.0000000000000000   0.0000000000000000  0.3999999999999999 # Al4
 Al   0.0000000000000000   0.0000000000000000  0.6000000000000000 # Al5 % NUM=1 ADATOM ZAMP=2.0
%ENDBLOCK POSITIONS_FRAC

FIX_ALL_CELL : true

#SUPERCELL=2 2 1
#SLAB
#POSAMP=0
#MINSEP=2.2
host:1.7 cjp10$ 

The shape of the primitive slab (atoms Al1-4) is fixed by "#FIX" in the last line of the LATTICE_CART block. The atom Al5 is an ADATOM. This means that is it added after the slab supercell has been built. The supercell is generated using "#SUPERCELL=2 2 1". In this case a 2x2x1 supercell is generated. You might also use "#SUPERCELL=4", which generates a random supercell of area 4. If "#SLAB" is omitted then the random supercell has a random volume of 4 (i.e. 2x1x2 is permitted). The global POSAMP is 0.0, but the added atom has "#ZAMP=2", which means that the X and Y coordinates are completely random, but the Z amplitude is less than or equal to 2.0. Initial configurations with atomic separations of less than 2.2 are rejected. The shape of the unit cell is not optimised during the relaxation.

host:1.7 cjp10$ airss.pl -pp3 -max 10 -seed Al
host:1.7 cjp10$ ca -r
Al-87679-1265-4        0.00    35.580    -5.570 17 Al           P3m1       1
Al-87679-1265-2        0.00    35.580     0.223 17 Al           P3m1       1
Al-87679-1265-7        0.00    35.580     0.223 17 Al           P3m1       1
Al-87679-1265-9        0.00    35.580     0.223 17 Al           P3m1       1
Al-87679-1265-5        0.00    35.580     0.223 17 Al           P3m1       1
Al-87679-1265-6        0.00    35.580     0.223 17 Al           P3m1       1
Al-87679-1265-3        0.00    35.580     0.223 17 Al           P3m1       1
Al-87679-1265-8        0.00    35.580     0.223 17 Al           P3m1       1
Al-87679-1265-10       0.00    35.580     0.223 17 Al           P3m1       1
Al-87679-1265-1        0.00    35.580     0.223 17 Al           P3m1       1
host:1.7 cjp10$ 


