# Michal's ddp fork implementing separate cutoff radii for features by body order

Whenever an existing script was modified, the changes were made in a copy with a '_new' suffix, e.g. 'frank.f90' -> 'frank_new.f90'
All makefiles were modified to compile these new script versions. The following scripts were changed or added:

## ddp/src/frank_new - Modified frank code, main changes include:

rmax variable changed from real to a real(5) array where:

* rmax(1) is the largest of the subsequent cutoff radii, used to identify all atoms in the structure relevant for feature construction
* rmax(2) is the cutoff radius for 2-body interaction features
* rmax(3) is the cutoff radius for 3-body interaction features
* rmax(4) is the cutoff radius for 4-body interaction features
* rmax(5) is the cutoff radius for 5-body interaction features

Throughout the code, only the rmax and rmax squared values chosen by interaction body order are used for feature identification, computation, and scaling in the generate_features()

All 5 cutoff radii are saved when training data is written to a file in generate_features()

The highest cutoff radius rmax(1) is used to identify relevant ions in generate_supercluster()

## ddp/src/forge_new - Modified forge code, main changes include:

Added compatibility with ddp and feature file metadata containing 5 rmax values instead of one

## ddp/src/flock_new - Modified flock code, main changes include:

Added compatibility with ddp metadata containing 5 rmax values instead of one

## ddp/src/data_new

* data_read() method reads five separate cutoff radii instead of one from a saved training data file

* data_compatible() method reads and compares sets of five separate cutoff radii instead of one

## ddp/bin/franks_new

franks has been slightly modified to pass through 5 cutoff radius values instead of one

## ddp/bin/multiple-cutoff-test

Script for hyperparameter grid search. 
Search parameters are most easily tuned in the script. It should be run in a folder with existing data.res file containing training structures, and produces a set of folders with ddps trained with each combination of hyperparameters from a chosen range and saved training logs.
These can be used to evaluate the performance of networks trained with each parameter combination.

# End of Michal's ddp fork changelog
-------------------------