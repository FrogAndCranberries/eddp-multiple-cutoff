In this example we will use random searching to find the ground state of Lennard-Jones cluster with 38 atoms, and the relax and shake algorithm.

host:1.5 cjp10$ cat Al.cell 
%BLOCK LATTICE_CART
20 0 0
0 20 0
0 0 20
#FIX
%ENDBLOCK LATTICE_CART
 
%BLOCK POSITIONS_FRAC
Al 0.0 0.0 0.0 # Al1 % NUM=38
%ENDBLOCK POSITIONS_FRAC

FIX_ALL_CELL : true

#MINSEP=1.5
#CLUSTER
#POSAMP=3.75
host:1.5 cjp10$

38 atoms are placed uniformly in a sphere of radius 3.75, with any structures with nearest neighbour distances less than 1.5 rejected.

host:1.5 cjp10$ airss.pl -pp3 -cluster -max 3000 -seed Al &
host:1.5 cjp10$ ca -u 0.01 -s -cl -r -t 10
Al-66136-6501-1502    -173.928  38 Al           Oh        3  3000
Number of structures   :   3000
Number of compositions :      1
Al-66136-6501-1502     0.00   210.526    -4.577 38 Al           Oh         3
Al-66136-6501-1916     0.00   210.526     0.021 38 Al           Cs         4
Al-66136-6501-1043     0.00   210.526     0.028 38 Al           C1         1
Al-66136-6501-489      0.00   210.526     0.042 38 Al           Cs         1
Al-66136-6501-1476     0.00   210.526     0.054 38 Al           C1         3
Al-66136-6501-2801     0.00   210.526     0.055 38 Al           C1         1
Al-66136-6501-2204     0.00   210.526     0.062 38 Al           Cs         3
Al-66136-6501-1936     0.00   210.526     0.064 38 Al           C3v        1
Al-66136-6501-2012     0.00   210.526     0.066 38 Al           C1         1
Al-66136-6501-1839     0.00   210.526     0.066 38 Al           C1         1

The known Oh LJ38 cluster is located 3 times in 3000 attempts.

host:1.5 cjp10$ airss.pl -pp3 -cluster -max 3000 -num 100 -amp 1.0 -seed Al &
host:1.5 cjp10$

The relax and shake (RASH) algorithm is a mimimal parameter "learning" algorithm. The first step is a structural optimisation of a random structure - which becomes the current "best" structure. The resulting local minimum is "shaken" - a random displacement of amplitude "amp" is applied to all atoms (subject to the minimum separation constraint) and relaxed. This shaking of the best (or most stable) structure is repeated "num" times (after which a new random structure is generated) or until a better structure is found, and then itself shaken up to "num" times.

host:1.5 cjp10$ ca -u 0.01 -s -cl -r
Al-46696-5893-2640    -173.928  38 Al           Oh        1    15
Number of structures   :     15
Number of compositions :      1
Al-46696-5893-2640     0.00   210.526    -4.577 38 Al           Oh         1
Al-46696-5893-834      0.00   210.526     0.018 38 Al           C5v        1
Al-46696-5893-1036     0.00   210.526     0.021 38 Al           Cs         9
Al-46696-5893-2807     0.00   210.526     0.045 38 Al           C1         1
Al-46696-5893-1475     0.00   210.526     0.066 38 Al           C1         1
Al-46696-5893-1363     0.00   210.526     0.067 38 Al           C2v        2
host:1.5 cjp10$

By default the rejected structures are removed so we are only left with 15 distinct structures. The AIRSS script can be called with the "-track" flag to suppress this behaviour.

In the above example, in both random search, and RASH, the known ground state for LJ38 was located - with random search apparently more successful. The mean number of attempts until first encounter (Na) is around 10K for random search and 1K for RASH with the parameters used, and due to the exponential probability distribution of such success/fail computational experiments the standard deviation is also Na, and hence large. Many runs must be performed before a good estimate of Na is obtained and reasonable assessment of the relative performance of methods may be made. An advantage of the RASH algorithm is that the structural optimisations following a shake require fewer steps than the relaxation of a completely random structure. A disadvantage is that the values of "num" and "amp" must be determined - and if "num" is large (so that the search is stuck in unpromising basins) or "amp" is small (so the local minimum is never escaped on shaking) Na will tend to infinity. A good recommendation for "amp" is about half a typical bond length. For large systems RASH is a very effective way to obtain "quite a good solution quite quickly".

It is clear, for this problem, the use of symmetry explored in 1.4 is recommended.
