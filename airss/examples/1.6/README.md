In this example we will use random searching to find the ground state of the Lennard-Jones 56 atom cluster, using a prebuilt LJ55 icosahedral core.

host:1.6 cjp10$ cat Al.cell 
%BLOCK LATTICE_CART
20 0 0
0 20 0
0 0 20
#FIX
%ENDBLOCK LATTICE_CART
 
%BLOCK POSITIONS_FRAC
Al       0.1311103     0.0067124     0.1719006 # Ih
Al       0.0306298     0.0572550     0.1731731 # Ih
Al      -0.0702839     0.1069271     0.1719000 # Ih
Al      -0.0638148    -0.0053629     0.1731728 # Ih
Al      -0.0563753    -0.1175927     0.1719000 # Ih
Al       0.0376361    -0.0558453     0.1731731 # Ih
Al       0.0667831     0.0027273    -0.1693897 # Ih
Al       0.0593436     0.1149571    -0.1681170 # Ih
Al      -0.0346679     0.0532098    -0.1693899 # Ih
Al      -0.1281420    -0.0093480    -0.1681174 # Ih
Al      -0.0276615    -0.0598906    -0.1693899 # Ih
Al       0.0732522    -0.1095627    -0.1681169 # Ih
Al      -0.0278654     0.1454824    -0.1039667 # Ih
Al      -0.1223100     0.0828646    -0.1039667 # Ih
Al      -0.1109734    -0.1001356    -0.1039667 # Ih
Al      -0.0095222    -0.1506181    -0.1039667 # Ih
Al       0.1432920    -0.0493000    -0.1039658 # Ih
Al       0.1362857     0.0638000    -0.1039658 # Ih
Al       0.0124908     0.1479825     0.1077493 # Ih
Al      -0.1403238     0.0466647     0.1077489 # Ih
Al      -0.1333174    -0.0664356     0.1077490 # Ih
Al       0.0308333    -0.1481180     0.1077494 # Ih
Al       0.1252778    -0.0855000     0.1077497 # Ih
Al       0.1139417     0.0975000     0.1077496 # Ih
Al      -0.2082562    -0.0143111     0.0420247 # Ih
Al      -0.0748214    -0.1546632     0.0673151 # Ih
Al       0.1176064    -0.1764616     0.0420255 # Ih
Al       0.1724377     0.0092725     0.0673157 # Ih
Al       0.0951018     0.1868191     0.0420254 # Ih
Al      -0.0931643     0.1414373     0.0673150 # Ih
Al      -0.1694694    -0.0119081    -0.0635325 # Ih
Al      -0.0921333    -0.1894547    -0.0382422 # Ih
Al       0.0961326    -0.1440729    -0.0635318 # Ih
Al       0.2112244     0.0116752    -0.0382415 # Ih
Al       0.0777897     0.1520276    -0.0635319 # Ih
Al      -0.1146382     0.1738260    -0.0382424 # Ih
Al       0.0651309     0.0026250     0.0853667 # Ih
Al       0.0474506     0.0910579     0.0215974 # Ih
Al      -0.0337541     0.0518307     0.0853667 # Ih
Al      -0.1014989    -0.0076974     0.0215970 # Ih
Al      -0.0269250    -0.0584091     0.0853667 # Ih
Al       0.0585000    -0.0873138     0.0215975 # Ih
Al       0.1044667     0.0050618    -0.0178138 # Ih
Al       0.0298932     0.0557735    -0.0815830 # Ih
Al      -0.0555322     0.0846778    -0.0178143 # Ih
Al      -0.0621626    -0.0052606    -0.0815832 # Ih
Al      -0.0444824    -0.0936934    -0.0178142 # Ih
Al       0.0367222    -0.0544667    -0.0815830 # Ih
Al      -0.1513304    -0.1026356     0.0018912 # Ih
Al       0.1542987     0.1000000     0.0018919 # Ih
Al       0.0128207    -0.1843180     0.0018917 # Ih
Al       0.1656352    -0.0830000     0.0018920 # Ih
Al      -0.0098524     0.1816824     0.0018915 # Ih
Al      -0.1626667     0.0803646     0.0018911 # Ih
Al       0.0014841    -0.0013178     0.0018916 # Ih
Al       0.0000000     0.0000000     0.0000000 # Al % NUM=1 POSAMP=6 MINAMP=5
%ENDBLOCK POSITIONS_FRAC

FIX_ALL_CELL : true

#MINSEP=1.5
#CLUSTER
#POSAMP=0
#ANGAMP=0
host:1.6 cjp10$ 

We start with a 55 atom icosahedral (point group Ih) cluster, centred on the origin. On generating the random structure this cluster is neither shifted (#POSAMP=0) nor rotated (#ANGAMP=0). The "# Ih" label is identical for all the 55 atoms of the cluster. The 56th atom is randomly located in a shell of outer radius 6 (#POSAMP=6) and inner radius 5 (#MINAMP=5). Structures in which any atoms are closer than 1.5 are rejected.

host:1.6 cjp10$ airss.pl -pp3 -cluster -max 10 -seed Al
host:1.6 cjp10$
host:1.6 cjp10$ ca -s -cl -r
Al-52722-9075-6       -283.643  56 Al           C3v       1    10
Number of structures   :     10
Number of compositions :      1
Al-52722-9075-6        0.00   142.857    -5.065 56 Al           C3v        1
Al-52722-9075-1        0.00   142.857     0.006 56 Al           Cs         1
Al-52722-9075-10       0.00   142.857     0.006 56 Al           Cs         1
Al-52722-9075-8        0.00   142.857     0.006 56 Al           Cs         1
Al-52722-9075-5        0.00   142.857     0.006 56 Al           Cs         1
Al-52722-9075-3        0.00   142.857     0.006 56 Al           Cs         1
Al-52722-9075-9        0.00   142.857     0.006 56 Al           Cs         1
Al-52722-9075-7        0.00   142.857     0.006 56 Al           Cs         1
Al-52722-9075-4        0.00   142.857     0.006 56 Al           Cs         1
Al-52722-9075-2        0.00   142.857     0.006 56 Al           Cs         1
host:1.6 cjp10$ 

In the above run the known ground state structure of LJ56 is located after 10 attempts.
