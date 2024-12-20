We now perform a ternary search.

host:1.10 cjp10$ airss.pl -pp3 -max 1000 -seed ABC

The output can be analysed using cryan. "ca -s" will give a summary of the covered compositions, and "ca -m" computes the convex hull. Both will generate a lot of output, so here we only consider those compositions on the convex hull, using the "-de 0" flag, which sets the energy window to zero.

host:1.10 cjp10$ ca -de 0 -m | sort -n -k 6 -k 5
ABC-20883-7172-331       0.00     8.849    -21.276  -4.371   0.000 +  4 BA           Pm-3m      1
ABC-20883-7172-679      -0.00    15.261    -38.136  -4.102  -0.000 +  1 CB2A         Fm-3m      1
ABC-20883-7172-932       0.00    34.558    -75.445  -3.686   0.000 +  1 B5A3         R-3m       1
ABC-20883-7172-215      -0.00    22.223    -51.280  -3.394  -0.000 +  1 C2B3A        P-3m1      1
ABC-20883-7172-880      -0.00    12.906    -26.894  -3.394   0.000 +  1 B2A          I4/mmm     1
ABC-20883-7172-887       0.00    10.612    -19.048  -2.839  -0.000 +  1 CA           Pm-3m      1
ABC-20883-7172-497      -0.00    16.996    -32.096  -2.801   0.000 +  1 B3A          I4/mmm     1
ABC-20883-7172-570      -0.00    20.333    -31.061  -1.916   0.000 +  2 C3A          Pmmn       1
ABC-20883-7172-656       0.00     7.004    -12.950  -1.879   0.000 +  1 CB           Pm-3m      1
ABC-20883-7172-383       0.00     5.013     -4.178   0.000   0.000 +  4 B            P63/mmc    1
ABC-20883-7172-948      -0.00     7.357     -8.356   0.000   0.000 +  2 A            P63/mmc    1
ABC-20883-7172-960       0.00     5.363     -5.014   0.000   0.000 +  2 C            P63/mmc    1

If R is installed, along with the ggtern package, then the terary convex hull can be visualised by running "Rscript ternary.r", and opening the resulting PDF file.
