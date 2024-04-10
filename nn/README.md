nn
==

Fortran implementation of neural networks

examples
========

    Chris-MBP:nn cjp20$ make && make install && make test
    (cd src/; make)
    make[1]: Entering directory '/Users/cjp20/Code/nn/src'
    gfortran -O3 -c nn.f90
    gfortran -O3 -c constants.f90
    gfortran -O3 -c nn_sine.f90
    gfortran nn.o constants.o nn_sine.o -o nn_sine -L/usr/lib -llapack -lblas
    gfortran -O3 -c nn_fourier.f90
    gfortran nn.o constants.o nn_fourier.o -o nn_fourier -L/usr/lib -llapack -lblas
    make[1]: Leaving directory '/Users/cjp20/Code/nn/src'
    (cp src/nn_sine bin/)
    (cp src/nn_fourier bin/)
    (cd examples && ../bin/nn_sine)
     nn_sine test
     Number of layers:             3
     Number of weights:          151
     Number of steps:          17750
     Cost for training set:    9.1484273430263039E-007
     Cost for testing set:     7.8453966073413421E-005
    (cd examples && ../bin/nn_fourier)
     nn_fourier test
     Number of inputs:             4
     Number of outputs:           21
     Number of layers:             3
     Number of weights:          151
     Number of steps:           7311
     Cost for training set:    5.6525023826135419E-005
     Cost for testing set:     4.3284923486231256E-004
    Chris-MBP:nn cjp20$

