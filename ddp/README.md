# ddp - generate data derived potentials

The `ddp` package contains a suite of tools to construct and test data
derived interatomic potentials. They are designed to be used with the
`airss` first principles structure prediction package. `airss` can be
used to generate data, and can exploit the generated ddp potentials to
potentially accelerate searches.

## package contents

* `frank` - generate features from `airss` .res file

* `forge` - creates a potential

* `flock` - construct an ensemble potential

* `farm` - dispatch multiple instances of `forge`

`franks` is a script to help the computation of features using
`frank` and combining them into training, validation and test datasets.

## prerequisites

The `airss` and `repose` packages should be installed, as well as the source for the
`nn` package. Various command line tools, such as `gnu parallel` are
used. The `grace` plotting package is supported. `gfortran 9.3` or above is
recommended, but `ifort` is supported.

## getting started

### compilation

The `ddp` package is compiled and installed using `make && make
install`, and if successful you should add `ddp/bin` to your path.

### example

An example raw dataset is provided in `ddp/examples/C`. It has been generated
using the `airss` code, by randomly placing 8 carbon atoms in a cubic
unit cell. The single point energies are computed using the `castep`
code. The energy range of the data is large (up to 800 eV/atom).

#### preparing datasets

The training, validation and test feature sets can be constructed:

`franks 3.75 3 8 2 10 1000 100 100 100 0`

There are about 10K structures in the raw dataset, and `franks`
command will randomly select 1000 for the training set, and 100 each
for the validation and testing sets. A cutoff radius of 3.75 Ang is
used, and features considering up to three-body interactions are
constructed. The powers to which the generating function is raised range
from 2 to 10, and 8 such powers are requested. An energy window of 100
eV is applied, at a pressure of 0GPa.

#### fitting the potential

To fit the potential, execute the `forge` command. By default a purely
linear potential will be generated, and saved to `forge.ddp`. You can
examine the quality of the fit using `xmgrace forge-testing.agr`.

#### ensemble potentials

It may be possible to improve the quality of the fit by using a neural
network to introduce nonlinearity. Rerun the fit using `forge -nn 5`,
for a neural network with a single hidden layer of 5 nodes. Non-linear 
fitting requires non-convex optimisation, and the quality of the fit can 
be variable. The overall quality and robustness of the potential can be 
improved by combining the results of several fits:

`for i in {1..10}; do echo ====$i==== && forge -nn 5 -ompnp 4 -q -np
-s C.$i ; done`

The arguments `-ompnp 4` request four `openmp` cores (note: not `MPI` cores),
`-q` suppresses output, and `-np` prevents the generation of output
data. The `-s C.$i` specifies the seedname, appended by the number of the
potential.

The `flock` command can be used to combine the ten generated
potentials:

`ls C.*.ddp | flock -v > C.eddp`

The `-v` option requests the optimisation based on the validation
dataset, rather than the training set.

#### using the potential

The `.ddp` and `.eddp` potentials can be used by the `repose` code to perform local structural optimisations, and the `airss` package can
perform searches using the `repose` code, for example:

`airss.pl -repose -max 10 -seed C && cat C-*.res | cryan -r`

Given the nature of the training dataset it is unlikely that a high quality potential will have been generated. However, it should be rather *robust* and permit searching.

#### generating datasets

The quality and the utility of a given data derived potential is strongly dependent on the training dataset. Possible sources of 
data are i) entirely random, unrelaxed, structures, ii) partly relaxed structures, iii) shaken relaxed stuctures, or iv) a collection of snapshots 
during structural optimisation. These snapshots can be generated using the `-harvest` argument to `airss.pl`.

#### iterative fitting

The provided `chain` script implements an iterative scheme for fitting potentials. The input is a `<seed>.param/.cell` pair which would be used for the AIRSS search of interest. The `<seed>.param` file should be edited to turn off geometry optimisation by `Castep`. The script uses the `spawn/farm` commands to parallelise the generation process, which manage the running of commands on the cpus listed in your `.spawn` file.

The command:

`nohup chain -s <seed> & disown`

will run the `chain` script in the background. Defaults have been chosen for the settings, and they can be changed using arguments to the `chain` script. For example:

`nohup chain -s <seed> -f 10000 -r 5 & disown`

generates 10K rully random structures before the iteration begins, and a cutoff radius of 5 Ang is used. If a directory `marks` exists, and it contains structures in the `.res` format, they will be shaken and their energies computed using `Castep` as marker structures.






