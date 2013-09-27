import os
import sys
from collections import defaultdict
from pprint import pprint
import distutils
# Desmond python tools
DESMOND_PATH = distutils.spawn.find_executable('desmond')
if DESMOND_PATH is None:
    raise ValueError('Could not find desmond in your path')

sys.path.insert(0, os.path.join(os.path.dirname(DESMOND_PATH), '..', 'lib', 'python'))
import framesettools
import molfile

import numpy as np
import matplotlib.pyplot as pp
from simtk import unit
from simtk import openmm as mm
from simtk.openmm import app
MAGIC = -12345

if 'GMXDATA' not in os.environ:
    raise ValueError('Could not find gromacs env variable GMXDATA to locate ff include dir')
        
gro = app.GromacsGroFile('conf.gro')
top = app.GromacsTopFile('topol.top', unitCellDimensions=gro.getUnitCellDimensions(), includeDir=os.path.join(os.environ['GMXDATA'], 'top'),
                         defines={'FLEXIBLE': True})

n_atoms = len(list(top.topology.atoms()))
print 'Number of bonds in topol.top: %d' % len(list(top.topology.bonds()))
print 'Number of atoms in topol.top: %d' % n_atoms
print

# Rename residues
table = {'HB2': 'HB1', 'HB3': 'HB2', 'HA2': 'HA1', 'HA3': 'HA2', 'HG2': 'HG1', 'HG3': 'HG2'}
for atom in top.topology.atoms():
    if atom.name in table:
        atom.name = table[atom.name]

# Rename H to H1 in residue 1
for atom in list(top.topology.residues())[0].atoms():
    if atom.name == 'H':
        atom.name = 'H1'

system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=0.9*unit.nanometers, rigidWater=True, constraints=app.HBonds)
context = mm.Context(system, mm.VerletIntegrator(0))
print 'Number of constraints', system.getNumConstraints()

openmm_residues = {residue.index: [(atom.name, atom.index) for atom in residue.atoms()] for residue in top.topology.residues()}
print 'OpenMM Atoms', sum(len(v) for k, v in openmm_residues.items())
desmond_residues = defaultdict(lambda : [])
for i, atom in enumerate(molfile.dms.read('1PLX.dms').atoms):
    desmond_residues[atom.resid-1].append((atom.name, i))

print 'OpenMM Atoms', sum(len(v) for k, v in openmm_residues.items())
print 'Desmond Atoms', sum(len(v) for k, v in desmond_residues.items())

assert sorted(openmm_residues.keys()) == sorted(desmond_residues.keys())
assert n_atoms == len(molfile.dms.read('1PLX.dms').atoms)
desmond_to_mm = MAGIC*np.ones(n_atoms, dtype=np.int)

# Build a mapping from the atom indices in the desmond
# trajectory to OpenMM
for i in openmm_residues.keys():
    openmm_r = sorted(openmm_residues[i])
    desmond_r  = sorted(desmond_residues[i])
    print 'OpenMM ', [e[0] for e in openmm_r]
    print 'Desmond', [e[0] for e in desmond_r]
    match = [e[0] for e in openmm_r] == [e[0] for e in desmond_r]
    print 'Atom nams in residue %d  match: %s ' % (i, match)
    assert match

    for n, m in zip([e[1] for e in desmond_r], [e[1] for e in openmm_r]):
        desmond_to_mm[m] = n

# Make sure we got everything
assert MAGIC not in desmond_to_mm

data = {'force_openmm': [],
        'force_desmond': [],
        'element': [],
        }

for f1, f2 in zip(framesettools.FrameSet('trajectory.dtr'), framesettools.FrameSet('forces.dtr')):
    #print f1.__labels__
    #print f2.__labels__
    #print 

    assert n_atoms == len(f1.POSITION)/3
    positions = unit.Quantity(f1.POSITION.reshape(n_atoms, 3)[desmond_to_mm], unit.angstroms)
    forces = unit.Quantity(f2.FORCES.reshape(n_atoms, 3)[desmond_to_mm], unit.kilocalories_per_mole/unit.angstrom)

    with open('topology-match.pdb', 'w') as f:
        app.PDBFile.writeFile(top.topology, positions, f)
    
    #print positions
    context.setPositions(positions)
    state = context.getState(getForces=True, getEnergy=True)
    print 'Potential Enegy (OpenMM)   ', state.getPotentialEnergy()
    print 'Potential Energy (desmond) ', unit.Quantity(f1.POT_ENERGY[0], unit.kilocalories_per_mole).in_units_of(unit.kilojoule_per_mole)

    print 'Forces (OpenMM) ', state.getForces()[0]
    print 'Forces (Desmond)', forces[0].in_units_of(unit.kilojoule_per_mole/unit.nanometers)

    for i, atom in enumerate(top.topology.atoms()):
        for j in range(3):
            data['force_openmm'].append(state.getForces()[i][j].value_in_unit_system(unit.md_unit_system))
            data['force_desmond'].append(forces[i, j].value_in_unit_system(unit.md_unit_system))
            data['element'].append(atom.element.symbol)

for k in data.keys():
    data[k] = np.array(data[k])

print np.mean(np.abs(data['force_openmm']))

color = {'C': 'g', 'N': 'b', 'O': 'r', 'H': 'white', 'S': 'y'}
colorseq = [color[d] for d in data['element']]

print 'Plotting'
pp.scatter(np.abs(data['force_openmm']), np.abs(data['force_desmond']), c=colorseq)
max_xy = max(max(data['force_openmm']), max(data['force_desmond']))
pp.plot([0, max_xy], [0, max_xy], 'k-')  # plot x==y
pp.xlim(0, max_xy)
pp.ylim(0, max_xy)
pp.xlabel('OpenMM  Force (kJ/mol/nm)')
pp.xlabel('Desmond Force (kJ/mol/nm)')
pp.title('OpenMM/Desmond Forces: 1PLX, charmm22*')
pp.savefig('force_error.png')



pp.figure()
error = np.abs(data['force_openmm'] - data['force_desmond'])
pp.hist(error, bins=np.logspace(-3, np.log10(np.max(error)), 50))
pp.xscale('log')

pp.title('OpenMM/Desmond Force Error: 1PLX, charmm22*')
pp.xlabel('Unsigned Error (kJ/mol/nm)')
pp.ylabel('Freq.')
pp.savefig('force_error_hist.png')


