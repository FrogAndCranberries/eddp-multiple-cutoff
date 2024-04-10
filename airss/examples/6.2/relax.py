from ase import Atoms
from ase.constraints import ExpCellFilter
import ase.io
from ase.optimize import LBFGSLineSearch
import numpy as np
from quippy.potential import Potential
import sys
to_eV_per_A3 = 1.0 / 160.2177 # from GPa

#------------------------------------------------

# Collect command line input and load random structure
pressure = float(sys.argv[1])
seed = sys.argv[2]
atoms = ase.io.read(seed + '.xyz')

# Relax cell using GAP
gap_file = 'Si_PRX_GAP/gp_iter6_sparse9k.xml'
atoms.calc = Potential(param_filename=gap_file)
ecf = ExpCellFilter(atoms, scalar_pressure=pressure*to_eV_per_A3)
optimizer = LBFGSLineSearch(ecf)
optimizer.run(2e-4, 1000)

# If converged, output final structure and fake castep file
if optimizer.converged():    
    ase.io.write(seed + '-out.xyz', atoms, 'extxyz')
    volume = atoms.get_volume()
    pv = pressure*to_eV_per_A3 * volume
    enthalpy = atoms.get_potential_energy() + pv
    with open(seed + '.castep', 'w') as f:
        f.write("Current cell volume = {:25.15f} A**3\n".format(volume))
        f.write("*  Pressure:   {:25.15f}\n".format(pressure))
        f.write("Python: Final Enthalpy     = {:25.15f} eV\n".format(enthalpy))
