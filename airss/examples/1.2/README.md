In this example we will use random searching to find the ground state of a Lennard-Jones cluster.

host:1.2 cjp10$ ls
Al.cell	Al.pp	README
host:1.2 cjp10$ cat Al.cell 
%BLOCK LATTICE_CART
20 0 0
0 20 0
0 0 20
#FIX
%ENDBLOCK LATTICE_CART
 
%BLOCK POSITIONS_FRAC
Al 0.0 0.0 0.0 # Al1 % NUM=13
%ENDBLOCK POSITIONS_FRAC

FIX_ALL_CELL : true

#MINSEP=1.5
#CLUSTER
#POSAMP=3.0
host:1.2 cjp10$ 

13 atoms are placed in a large fixed box. Centred on the origin, they are uniformly distributed over a sphere of radius 3.0 Ang. If any atoms are closer than 1.5 Ang the configuration is rejected.

We now perform the search:

host:1.2 cjp10$ airss.pl -pp3 -cluster -max 20 -seed Al
host:1.2 cjp10$ ca -r
Al-72120-6057-19         0.00   615.385      -3.410  13 Al           Ih         1
Al-72120-6057-4          0.00   615.385       0.000  13 Al           Ih         1
Al-72120-6057-7          0.00   615.385       0.220  13 Al           Cs         1
Al-72120-6057-8          0.00   615.385       0.220  13 Al           Cs         1
Al-72120-6057-3          0.00   615.385       0.220  13 Al           Cs         1
Al-72120-6057-14         0.00   615.385       0.274  13 Al           C2v        1
Al-72120-6057-5          0.00   615.385       0.285  13 Al           C1         1
Al-72120-6057-10         0.00   615.385       0.285  13 Al           C1         1
Al-72120-6057-16         0.00   615.385       0.323  13 Al           Cs         1
Al-72120-6057-15         0.00   615.385       0.327  13 Al           Cs         1
Al-72120-6057-9          0.00   615.385       0.348  13 Al           C1         1
Al-72120-6057-12         0.00   615.385       0.351  13 Al           C1         1
Al-72120-6057-13         0.00   615.385       0.354  13 Al           C1         1
Al-72120-6057-11         0.00   615.385       0.354  13 Al           C1         1
Al-72120-6057-2          0.00   615.385       0.356  13 Al           C1         1
Al-72120-6057-18         0.00   615.385       0.359  13 Al           C1         1
Al-72120-6057-20         0.00   615.385       0.390  13 Al           C1         1
Al-72120-6057-6          0.00   615.385       0.416  13 Al           C1         1
Al-72120-6057-1          0.00   615.385       0.439  13 Al           C1         1
Al-72120-6057-17         0.00   615.385       0.490  13 Al           C1         1
host:1.2 cjp10$ ca -u 0.01 -r
Al-72120-6057-19         0.00   615.385      -3.410  13 Al           Ih         2
Al-72120-6057-8          0.00   615.385       0.220  13 Al           Cs         3
Al-72120-6057-14         0.00   615.385       0.274  13 Al           C2v        1
Al-72120-6057-5          0.00   615.385       0.285  13 Al           C1         2
Al-72120-6057-16         0.00   615.385       0.323  13 Al           Cs         1
Al-72120-6057-15         0.00   615.385       0.327  13 Al           Cs         1
Al-72120-6057-9          0.00   615.385       0.348  13 Al           C1         1
Al-72120-6057-12         0.00   615.385       0.351  13 Al           C1         1
Al-72120-6057-13         0.00   615.385       0.354  13 Al           C1         1
Al-72120-6057-11         0.00   615.385       0.354  13 Al           C1         1
Al-72120-6057-2          0.00   615.385       0.356  13 Al           C1         1
Al-72120-6057-18         0.00   615.385       0.359  13 Al           C1         1
Al-72120-6057-20         0.00   615.385       0.390  13 Al           C1         1
Al-72120-6057-6          0.00   615.385       0.416  13 Al           C1         1
Al-72120-6057-1          0.00   615.385       0.439  13 Al           C1         1
Al-72120-6057-17         0.00   615.385       0.490  13 Al           C1         1
host:1.2 cjp10$

The known icosahedral ground state has been found 6 times in this test. The resulting structures are stored in SHLX format .res files and can be converted to other formats using the cabal utility.

host:1.2 cjp10$ cat Al-72120-6057-19.res
TITL Al-72120-6057-19 0.0000000000 8000.0000000000 -44.3248941562 0 0 13 (Ih) n - 1
REM
REM in /Users/user/examples/1.2
REM
REM
REM
REM
REM
REM buildcell < ./Al.cell (8f0fc5c9d882ed20c27978dd052da9d4)
REM AIRSS Version 0.9.1 July 2018 build 92bfd83db9d4+ Sat, 30 Jun 2018 13:13:37 +0100
REM compiler GCC version 5.5.0
REM options -fPIC -feliminate-unused-debug-symbols -mmacosx-version-min=10.13.6 -mtune=core2 -g -O0
REM seed -1471360510 667860809 1027640838 1038292373 -213802532 -1539206485 -1872642957 -340800072 -697171857 -761177086 -1542735574 -885338462
REM
CELL 1.54180   20.00000   20.00000   20.00000   90.00000   90.00000   90.00000
LATT -1
SFAC Al
Al     1  0.6038930453823  0.4836194826238  0.5253309492534 1.0
Al     1  0.5723861716435  0.5636803790552  0.4509205384933 1.0
Al     1  0.4999999997846  0.5000000000790  0.5000000075654 1.0
Al     1  0.4276138527499  0.4363195677348  0.5490794435604 1.0
Al     1  0.4701141676138  0.6024370215993  0.4821888755138 1.0
Al     1  0.5298858276893  0.3975629918030  0.5178112246213 1.0
Al     1  0.4526399878259  0.4244380397184  0.4387532754310 1.0
Al     1  0.5473600525560  0.5755619778738  0.5612466872006 1.0
Al     1  0.4789066108987  0.5271042987692  0.3974126364780 1.0
Al     1  0.3961069520357  0.5163804721092  0.4746690478913 1.0
Al     1  0.5210933793772  0.4728958126855  0.6025874103552 1.0
Al     1  0.4384134677060  0.5463289943787  0.5759240890051 1.0
Al     1  0.5615864847371  0.4536709615702  0.4240759146313 1.0
END
host:1.2 cjp10$

The first line (TITL) contains stored data in the following format:

TITL <name> <pressure> <volume> <enthalpy> <spin> <modspin> <#ions> <(symmetry)> n - <#copies>

host:1.2 cjp10$ cabal res xyz < Al-72120-6057-19.res
13

Al     1.1791750764873     1.7827258262017    -0.3360534169762
Al     0.7363476807396     1.0192632282492     1.7607929157278
Al    -0.0000000029634    -0.0000000036287     0.0000000023555
Al    -0.7363484783607    -1.0192621771574    -1.7607931876182
Al     0.5875415154852    -1.2489218391593     1.6662791865905
Al    -0.5875402622518     1.2489211530929    -1.6662801418123
Al    -2.1222590886862     0.1555499373659    -0.3915580422390
Al     2.1222589581482    -0.1555490740574     0.3915591037870
Al    -1.3040498658070     0.0136119341742     1.7264894398813
Al    -1.1791758452594    -1.7827253947705     0.3360529947887
Al     1.3040514186501    -0.0136131959079    -1.7264882570909
Al     0.9383998151967    -1.8872737156990    -0.4889793760508
Al    -0.9384009213786     1.8872733212961     0.4889787786564