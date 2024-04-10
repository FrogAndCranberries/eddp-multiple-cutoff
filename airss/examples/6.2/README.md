Example 6.2
===========

## Gaussian Approximation Potential (GAP) search for 4 atoms of Silicon

### Prerequisites

QUIP is required, with the GAP extension, as well as quippy.
Full installation instructions are [here](https://github.com/libAtoms/QUIP),
but it should be sufficient to install quippy via pip:

    host:1.1 cjp10$ pip install quippy-ase

GAP data files are also required, namely the Si_PRX_GAP potential available
[here](https://www.repository.cam.ac.uk/handle/1810/317974)). Place the Si_PRX_GAP
directory in the directory with this README:

    host:1.1 cjp10$ ls Si_PRX_GAP
    gp_iter6_sparse9k.xml
    gp_iter6_sparse9k.xml.sparseX.GAP_2017_6_17_60_4_3_56_1651
    gp_iter6_sparse9k.xml.xyz
    repulsive_2b.xml
    sparse_point_environment_and_config_type_list

### Performing the search

We now perform a search

    host:1.1 cjp10$ airss.pl -python -max 20 -seed Si
    host:1.1 cjp10$

obtaining

    host:1.1 cjp10$ ca -u 0.01 -r
    Si-10921-3998-10         0.00    20.357    -163.178   4 Si           Fd-3m          1
    Si-10921-3998-8          0.00    20.336       0.009   4 Si           P63/mmc        2
    Si-10921-3998-6          0.00    17.413       0.319   4 Si           C2/m           2
    Si-10921-3998-16         0.00    15.386       0.321   4 Si           I41/amd        1
    Si-10921-3998-9          0.00    19.927       0.331   4 Si           P4122          1
    Si-10921-3998-20         0.00    18.469       0.333   4 Si           C2/m           1
    Si-10921-3998-4          0.00    19.095       0.334   4 Si           I212121        1
    Si-10921-3998-14         0.00    15.144       0.342   4 Si           P6/mmm         1
    Si-10921-3998-13         0.00    15.530       0.347   4 Si           Cmmm           1
    Si-10921-3998-1          0.00    15.563       0.362   4 Si           Pbcm           1
    Si-10921-3998-3          0.00    20.006       0.387   4 Si           C2             1
    Si-10921-3998-11         0.00    17.891       0.407   4 Si           C2/m           1
    Si-10921-3998-19         0.00    21.727       0.409   4 Si           Imma           1
    Si-10921-3998-12         0.00    21.554       0.424   4 Si           Cmmm           1
    Si-10921-3998-2          0.00    34.552       0.496   4 Si           P-1            1
    Si-10921-3998-18         0.00    42.878       0.574   4 Si           Pm             1
    Si-10921-3998-17         0.00    28.126       0.763   4 Si           Cm             1
    Si-10921-3998-15         0.00    29.692       0.843   4 Si           Cm             1
    host:1.1 cjp10$

This search did yield the expected cubic diamond ground state
(space group Fd-3m), but only once. In general, a larger search may
be required to find it.
