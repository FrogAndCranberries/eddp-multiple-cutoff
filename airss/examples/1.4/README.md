In this example we will use random searching to find the ground state of Lennard-Jones cluster with 38 atoms, using symmetry.

host:1.4 cjp10$ cat Al.cell 
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
#SYMMOPS=2-5
#NFORM=1
host:1.4 cjp10$

Random clusters will be built containing 38 atoms and with the symmetry of a randomly chosen point group with between 2 and 5 symmetry operations. If "#NFORM=1" is omitted the default behaviour is that the number of atoms in the cluster will be nops x 38, which is unlikely to be satisfied in a sphere of radius 3.75 and minimum separation of 1.5.

host:1.4 cjp10$ airss.pl -pp3 -cluster -max 200 -seed Al &
host:1.4 cjp10$ 

The searches can be run in the background and monitored during their progress using the "ca" command. 

host:1.4 cjp10$ ca -u 0.01 -s -cl -r -t 10
Al-65004-3113-183     -173.928  38 Al           Oh        2   200
Number of structures   :    200
Number of compositions :      1
Al-65004-3113-183      0.00   210.526    -4.577 38 Al           Oh         2
Al-65004-3113-100      0.00   210.526     0.018 38 Al           C5v        3
Al-65004-3113-14       0.00   210.526     0.026 38 Al           C5v        7
Al-65004-3113-128      0.00   210.526     0.078 38 Al           D4h        5
Al-65004-3113-111      0.00   210.526     0.078 38 Al           C3         1
Al-65004-3113-110      0.00   210.526     0.083 38 Al           C2h        1
Al-65004-3113-33       0.00   210.526     0.097 38 Al           C1         1
Al-65004-3113-50       0.00   210.526     0.100 38 Al           C2         1
Al-65004-3113-164      0.00   210.526     0.100 38 Al           D2         2
Al-65004-3113-65       0.00   210.526     0.105 38 Al           C1         1
host:1.4 cjp10$

The above is a summary of the completed run, followed by the top ten (-t) structures. The known Oh ground state of LJ38 is found twice in 200 attempts in this example, and the C5v icosahedral minimum is also located. The use of symmetry is highly recommended - compare to the results of example 1.5.

