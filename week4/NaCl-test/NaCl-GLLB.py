from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import write, read
from ase.optimize import QuasiNewton
from ase.units import Bohr
from ase.visualize import view
from gpaw import GPAW, FermiDirac, PoissonSolver, Mixer, setup_paths, restart
from gpaw.xc.functional import XCFunctional
from gpaw.xc.hybridg import HybridXC
from gpaw.xc.sic import SIC
from gpaw.test import equal
from gpaw.xc.tools import vxc
from gpaw.wavefunctions.pw import PW
from math import sqrt
setup_paths.insert(0, '.')

'''
### PBE
pair = Atoms('NaCl', positions=[(0, 0, 0),(0, 0, 2.82)])
pair.set_pbc((False, False, False))
pair.center(vacuum=5)

calc = GPAW(xc='PBE',
            mode='fd',
            h=0.14,
            txt='NaCl-PBE.txt',
            basis='szp(dzp)')
pair.set_calculator(calc)
pair.get_potential_energy()
pair.calc.write('NaCl-PBE.gpw', mode='all')
density = calc.get_all_electron_density(gridrefinement=4) * Bohr**3
write('NaCl-PBE.cube', pair, data=density)

### PBE-SIC
pair, calc = restart('NaCl-PBE.gpw')

calc = GPAW(xc='PBE-PZ-SIC',
            mode='fd',
            h=0.14,
            eigensolver='rmm-diis',
            txt='NaCl-PBE-SIC.txt',
           )
pair.set_calculator(calc)
pair.get_potential_energy()
pair.calc.write('NaCl-SIC.gpw', mode='all')
density = calc.get_all_electron_density(gridrefinement=4) * Bohr**3
write('NaCl-SIC.cube', pair, data=density)
'''
'''

### PBEsol
#gpaw-setup Na -f PBEsol
#gpaw-setup Cl -f PBEsol

pair = Atoms('NaCl', positions=[(0, 0, 0),(0, 0, 2.82)])
pair.set_pbc((False, False, False))
pair.center(vacuum=7)

calc = GPAW(xc='PBEsol',
            mode='fd',
            h=0.14,
            txt='NaCl-PBEsol.txt',
            basis='szp(dzp)')
pair.set_calculator(calc)
pair.get_potential_energy()
pair.calc.write('NaCl-PBEsol.gpw', mode='all')
density = calc.get_all_electron_density(gridrefinement=4) * Bohr**3
write('NaCl-PBEsol.cube', pair, data=density)
'''
### GLLBC
pair, calc = restart('NaCl-PBEsol.gpw')

calc = GPAW(xc='GLLBSC',
            mode='fd',
            h=0.14,
            eigensolver='rmm-diis',
            txt='NaCl-GLLBSC.txt',
           )
pair.set_calculator(calc)
pair.get_potential_energy()
pair.calc.write('NaCl-GLLBSC.gpw', mode='all')
density = calc.get_all_electron_density(gridrefinement=4) * Bohr**3
write('NaCl-GLLBSC.cube', pair, data=density)
#'''
'''
### PBE0 OK
pair = Atoms('NaCl', positions=[(0, 0, 0),(0, 0, 2.82)])
pair.set_pbc((False, False, False))
pair.center(vacuum=3.5)

calc = GPAW(xc='PBE0',
            mode='fd',
            h=0.14,
            nbands=12,
            eigensolver='rmm-diis',
            convergence={'eigenstates': 1.0e-6},
            txt='NaCl-PBE0.txt',
            basis='szp(dzp)')
pair.set_calculator(calc)
pair.get_potential_energy()
pair.calc.write('NaCl-PBE0.gpw', mode='all')
density = calc.get_all_electron_density(gridrefinement=4) * Bohr**3
write('NaCl-PBE0.cube', pair, data=density)
'''
