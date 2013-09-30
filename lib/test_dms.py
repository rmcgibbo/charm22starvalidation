import os
import unittest
import distutils
import subprocess
import numpy as np

from dmsfile import DesmondDMSFile
from simtk import openmm as mm
from simtk.openmm import app

HAVE_VIPARR = distutils.spawn.find_executable('viparr')


@unittest.skipIf(not HAVE_VIPARR, 'Need viparr')
def test_amber99_constraints():
    subprocess.check_output('viparr ala2.pdb ala2.dms -f amber99SB-ILDN', shell=True)
    dms = DesmondDMSFile('ala2.dms')
    system1 = dms.createSystem()
    context1 = mm.Context(system1, mm.VerletIntegrator(0))
    context1.setPositions(dms.positions)
    state1 = context1.getState(getForces=True)
    force1 = state1.getForces(asNumpy=True)

    pdb = app.PDBFile('ala2.pdb')
    forcefield = app.ForceField('amber99sbildn.xml')
    system2 = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    context2 = mm.Context(system2, mm.VerletIntegrator(0))
    context2.setPositions(pdb.getPositions())
    state2 = context2.getState(getForces=True)
    force2 = state2.getForces(asNumpy=True)

    np.testing.assert_array_almost_equal(np.array(force1), np.array(force2), decimal=1)


@unittest.skipIf(not HAVE_VIPARR, 'Need viparr')
def test_amber99_no_constraints():
    subprocess.check_output('viparr ala2.pdb ala2.dms -f amber99SB-ILDN --without-constraints', shell=True)
    dms = DesmondDMSFile('ala2.dms')
    system1 = dms.createSystem()
    context1 = mm.Context(system1, mm.VerletIntegrator(0))
    context1.setPositions(dms.positions)
    state1 = context1.getState(getForces=True)
    force1 = state1.getForces(asNumpy=True)

    pdb = app.PDBFile('ala2.pdb')
    forcefield = app.ForceField('amber99sbildn.xml')
    system2 = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    context2 = mm.Context(system2, mm.VerletIntegrator(0))
    context2.setPositions(pdb.getPositions())
    state2 = context2.getState(getForces=True)
    force2 = state2.getForces(asNumpy=True)

    np.testing.assert_array_almost_equal(np.array(force1), np.array(force2), decimal=1)


@unittest.skipIf(not HAVE_VIPARR, 'Need viparr')
def test_1plx():
    subprocess.check_output('pdb2gmx -f 1PLX.pdb -ff charmm22star -ignh -water tip3p', shell=True)
    gro = app.GromacsGroFile('conf.gro')
    top = app.GromacsTopFile('topol.top')
    app.PDBFile.writeFile(top.topology, gro.positions, open('temp.pdb', 'w'))
    subprocess.check_output('viparr temp.pdb temp.dms -f charmm22star_aminoacids --without-constraints', shell=True)
    dms = DesmondDMSFile('temp.dms')
    system1 = dms.createSystem()
    context1 = mm.Context(system1, mm.VerletIntegrator(0))

    context1.setPositions(dms.positions)
    state1 = context1.getState(getForces=True)
    force1 = np.array(state1.getForces(asNumpy=True))


    system2 = top.createSystem(nonbondedMethod=app.NoCutoff, constraints=None)
    context2 = mm.Context(system2, mm.VerletIntegrator(0))
    context2.setPositions(gro.getPositions())
    state2 = context2.getState(getForces=True)
    force2 = np.array(state2.getForces(asNumpy=True))

    np.testing.assert_array_almost_equal(np.array(force1), np.array(force2), decimal=1)





@unittest.skipIf(not HAVE_VIPARR, 'Need viparr')
def test_charmm22star():
    subprocess.check_output('viparr ala2.pdb ala2.dms -f charmm22star_aminoacids', shell=True)
    dms = DesmondDMSFile('ala2.dms')
    system1 = dms.createSystem()
    context1 = mm.Context(system1, mm.VerletIntegrator(0))
    context1.setPositions(dms.positions)
    state1 = context1.getState(getForces=True)
    force1 = state1.getForces(asNumpy=True)

    print force1

@unittest.skipIf(not HAVE_VIPARR, 'Need viparr')
def test_charmm22star_1lyd():
    subprocess.check_output('viparr 1LYD.pdb 1LYD.dms -f charmm22star_aminoacids --without-constraints', shell=True)
    dms = DesmondDMSFile('1LYD.dms')
    system1 = dms.createSystem()
    context1 = mm.Context(system1, mm.VerletIntegrator(0))
    context1.setPositions(dms.positions)
    state1 = context1.getState(getForces=True)
    force1 = state1.getForces(asNumpy=True)

    print force1



def clone_system_with(system, *forces):
    new_system = mm.System()
    for i in range(system.getNumParticles()):
        new_system.addParticle(system.getParticleMass(i))
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if any(isinstance(force, f) for f in forces) or len(forces) == 0:
            new_system.addForce(force.__class__(force))
        else:
            print 'Discarding %s' % force.__class__.__name__
    return new_system

def print_cmap(cmap):
    print 'NUMBER OF MAPS'
    print cmap.getNumMaps()
    for i in range(cmap.getNumMaps()):
        print 'ROW 0 of MAP %d' % i
        size, energy = cmap.getMapParameters(i)
        print np.array(energy).reshape(size, size).T[0]
