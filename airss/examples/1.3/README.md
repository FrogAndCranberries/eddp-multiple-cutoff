In this example we will use random searching to find the ground state of Lennard-Jones clusters of a range of sizes.

host:1.3 cjp10$ cat Al.cell 
%BLOCK LATTICE_CART
20 0 0
0 20 0
0 0 20
#FIX
%ENDBLOCK LATTICE_CART
 
%BLOCK POSITIONS_FRAC
Al 0.0 0.0 0.0 # Al1 % NUM=8-16
%ENDBLOCK POSITIONS_FRAC

FIX_ALL_CELL : true

#MINSEP=1.5
#CLUSTER
#POSAMP=3.0
host:1.3 cjp10$

This is the same seed file as for example 1.2, except this time a range of clusters sizes is explored.

host:1.3 cjp10$ airss.pl -pp3 -cluster -max 200 -seed Al
host:1.3 cjp10$ ca -s
Al-70852-9173-163        0.00   500.000      -3.551  16 Al           Cs        1   200      ~
Number of structures   :    200
Number of compositions :      1

host:1.3 cjp10$

The command "ca -s" produces a summary of the results indicating the most stable structure (judged by energy/fu), which in this case is obviously the largest cluster.

host:1.3 cjp10$ ca -u 0.01 -s -cl
Al-70852-9173-70         -16.505    7 Al           D5h       5    19      ~
Al-70852-9173-99         -19.821    8 Al           Cs       16    22    -0.976
Al-70852-9173-167        -24.112    9 Al           C2v       2    15    -0.017
Al-70852-9173-146        -28.421   10 Al           C3v       5    26    -0.034
Al-70852-9173-43         -32.765   11 Al           C2v       5    22    -0.858
Al-70852-9173-77         -37.966   12 Al           C5v       2    13    -1.158
Al-70852-9173-174        -44.325   13 Al           Ih        3    23     2.841
Al-70852-9173-19         -47.843   14 Al           C3v       5    21    -0.959
Al-70852-9173-144        -52.320   15 Al           C2v       6    22    -0.016
Al-70852-9173-163        -56.813   16 Al           Cs        4    17      ~
Number of structures   :    200
Number of compositions :     10
host:1.3 cjp10$

The above is another summary of the results, this time in cluster mode, indicating the most stable structures for each cluster size. The sixth and seventh columns give the number of times the most stable structure is found, and the total number of structures. You may check your results using this excellent resource:

http://doye.chem.ox.ac.uk/jon/structures/LJ/tables.150.html

Plotting the energy/atom as a function of cluster size is a good way to identify "magic number", or particularly stable, clusters.

host:1.3 cjp10$ ca -u 0.01 -s -cl | awk '{print $3,$2/$3}' | xmgrace -pipe

It should be clear from your results that the 13 atom Ih symmetry cluster is particularly stable. Indeed, the final column of the summary presents the second derivative of the total energy with cluster size. The curvature is strongly positive for the 13 atom cluster, indicating high stability as compared to its neighbours.



