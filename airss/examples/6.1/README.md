Example 6.1
===========

## Lennard-Jones solid using the Atomic Simulation Environment

This example, like Example 1.1, uses random searching to find the ground state
of a Lennard-Jones solid. However, it uses the Lennard-Jones potential from
the Atomic Simulation Environment (ASE), with structural optimisations
controlled by the `relax.py` script. The seed file is identical to that from
Example 1.1.

We use the `-python` flag to perform the search

    host:1.1 cjp10$ airss.pl -python -max 20 -seed Al
    host:1.1 cjp10$ 

and the resulting structures are

    host:1.1 cjp10$ ca -u 0.01 -r
    Al-31158-3253-14         0.00     7.459      -7.472   8 Al           P63/mmc        2
    Al-31158-3253-12         0.00     7.450       0.004   8 Al           R-3m           2
    Al-31158-3253-19         0.00     7.442       0.009   8 Al           Fm-3m         12
    Al-31158-3253-10         0.00     7.672       0.286   8 Al           C2/m           1
    Al-31158-3253-20         0.00     8.338       0.800   8 Al           Cm             1
    Al-31158-3253-13         0.00     8.344       0.803   8 Al           C2/m           1
    Al-31158-3253-18         0.00     8.361       0.906   8 Al           Cm             1
    host:1.1 cjp10$

While the Lennard-Jones parameters set in `relax.py` are the same as those from
Example 1.1, the results here differ slightly (e.g., the volume and energy for
the HCP structure). To obtain identical results, one must redo Example 1.1 using
`pp3 -e`, where the `-e` flag disables force shifting employed by the pp3 code.
