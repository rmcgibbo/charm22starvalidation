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
