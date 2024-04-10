import numpy as np
import os
import sys

from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.constraints import ExpCellFilter
from ase.optimize import LBFGSLineSearch
import ase.io

#------------------------------------------------

# Collect command line input and load random structure
pressure = float(sys.argv[1])
seed = sys.argv[2]
atoms = ase.io.read(seed + '.xyz')

# Relax cell using Lennard-Jones potential
epsilon = 1
sigma = 2
cutoff = 2.5*sigma
atoms.calc = LennardJones(sigma=sigma, epsilon=epsilon, rc=cutoff)
ecf = ExpCellFilter(atoms, scalar_pressure=pressure)
opt = LBFGSLineSearch(ecf)	
converged = opt.run(steps=1000)

# If converged, output final structure and fake castep file
if converged:
    ase.io.write(seed + '-out.xyz', atoms, 'extxyz')
    volume = atoms.get_volume()
    enthalpy = atoms.get_potential_energy() + pressure*volume
    with open(seed + '.castep', 'w') as f:
        f.write("Current cell volume = {:25.15f} A**3\n".format(volume))
        f.write("*  Pressure:   {:25.15f}\n".format(pressure))
        f.write("Python: Final Enthalpy     = {:25.15f} eV\n".format(enthalpy))
