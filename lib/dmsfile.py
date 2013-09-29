import sqlite3
import warnings
from pprint import pprint

import networkx as nx


from simtk.unit import *
from simtk import openmm as mm
from simtk.openmm import app

class DesmondDMSFile(object):
    dmstype = {'integer': int, 'text': str, 'float': float}

    def __init__(self, file):
        self._open = False
        self.c = sqlite3.connect(file)
        self._open = True

        tables = {}
        for table in self.c.execute("SELECT name FROM sqlite_master WHERE type='table'"):
            names = []
            types = []
            for e in self.c.execute('PRAGMA table_info(%s)' % table):
                names.append(str(e[1]))
                types.append(str(e[2]))
            tables[str(table[0])] = {'names': names, 'types': types}
        self._tables = tables
        pprint(self._tables)


    def createSystem(self):
        sys = mm.System()

        # Buld the box dimensions
        boxVectors = []
        for i, x, y, z in self.c.execute('select * FROM global_cell'):
            boxVectors.append(mm.Vec3(x, y, z)*angstroms)
        sys.setDefaultPeriodicBoxVectors(*boxVectors)

        # Create all of the particles
        for term in self._iterRows('particle'):
            sys.addParticle(term['mass']*dalton)

        self._addBondsToSystem(sys)
        self._addAnglesToSystem(sys)
        self._addPeriodicTorsionsToSystem(sys)
        self._addImproperTorsionsToSystem(sys)
        self._addCMAPToSystem(sys)
        self._addNonbondedForceToSystem(sys, None, None)

        return sys

    def _addBondsToSystem(self, sys):
        # Create the harmonic bonds
        bonds = mm.HarmonicBondForce()
        sys.addForce(bonds)
        for term in self._iterRows('stretch_harm_term'):
            r0, fc = self.c.execute('SELECT r0, fc FROM stretch_harm_param WHERE id=%s' % term['param']).fetchone()
            if not term['constrained']:
                # Desmond writes the harmonic bond force without 1/2
                # so we need to to double the force constant
                bonds.addBond(term['p0'], term['p1'], r0*angstroms,
                              2*fc*kilocalorie_per_mole/angstrom**2)
            else:
                sys.addConstraint(term['p0'], term['p1'], r0*angstoms)

    def _addAnglesToSystem(self, sys):
        # Create the harmonic angles
        angles = mm.HarmonicAngleForce()
        sys.addForce(angles)
        for term in self._iterRows('angle_harm_term'):
            theta0, fc = self.c.execute('SELECT theta0, fc FROM angle_harm_param WHERE id=%s' % term['param']).fetchone()
            if not term['constrained']:
                # Desmond writes the harmonic angle force without 1/2
                # so we need to to double the force constant
                angles.addAngle(term['p0'], term['p1'], term['p2'],
                                theta0*degrees, 2*fc*kilocalorie_per_mole/radians**2)
            else:
                raise NotImplementedError('Angle restraints not implemeneted yet')

    def _addPeriodicTorsionsToSystem(self, sys):
        # Create the torsion terms
        periodic = mm.PeriodicTorsionForce()
        sys.addForce(periodic)
        for term in self._iterRows('dihedral_trig_term'):
            phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 = self.c.execute('''
                SELECT phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6
                FROM dihedral_trig_param
                WHERE id=%s''' % term['param']).fetchone()

            for order, fc in enumerate([fc0, fc1, fc2, fc3, fc4, fc5, fc6]):
                if fc == 0:
                    continue

                periodic.addTorsion(term['p0'], term['p1'], term['p2'], term['p3'],
                                    order, phi0*degree, fc*kilocalorie_per_mole)


    def _addImproperTorsionsToSystem(self, sys):
        # Create the improper harmonic torsion terms
        if not self._hasTable('improper_harm_term'):
            return

        harmonicTorsion = mm.CustomTorsionForce('0.5*k*(theta-theta0)^2')
        harmonicTorsion.addPerTorsionParameter('theta0')
        harmonicTorsion.addPerTorsionParameter('k')
        sys.addForce(harmonicTorsion)
        for term in self._iterRows('improper_harm_term'):
            phi0, fc = self.c.execute('SELECT phi0, fc FROM improper_harm_param WHERE id=%s' % term['param']).fetchone()
            harmonicTorsion.addTorsion(term['p0'], term['p1'], term['p2'],
                                       term['p3'], [phi0*degree, fc*kilocalorie_per_mole])

    def _addCMAPToSystem(self, sys):
        if not self._hasTable('torsiontorsion_cmap_term'):
            return

        # Create CMAP torsion terms
        cmap = mm.CMAPTorsionForce()
        sys.addForce(cmap)
        cmap_indices = {}

        for name in [k for k in self._tables.keys() if k.startswith('cmap')]:
            size2 = self.c.execute('SELECT COUNT(*) FROM %s' % name).fetchone()[0]
            fsize = math.sqrt(size2)
            if fsize != int(fsize):
                raise ValueError('Non-square CMAPs are not supported')
            size = int(fsize)

            map = [0 for i in range(size2)]
            for term in self._iterRows(name):
                i = int((size*term['phi']/360.0 + size/2) % size)
                j = int((size*term['psi']/360.0 + size/2) % size)
                map[i+size*j] = term['energy']*kilocalories_per_mole
            index = cmap.addMap(size, map)
            cmap_indices[name] = index


        for term in self._iterRows('torsiontorsion_cmap_term'):
            cmapid = self.c.execute('''SELECT cmapid FROM torsiontorsion_cmap_param
                                       WHERE id=%s''' %  term['param']).fetchone()[0]
            cmap.addTorsion(cmap_indices[cmapid],
                            term['p0'], term['p1'], term['p2'], term['p3'],
                            term['p4'], term['p5'], term['p6'], term['p7'])

    def _addNonbondedForceToSystem(self, sys, nonbondedMethod, nonbondedCutoff):
        # Create the nonbonded force
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        sys.addForce(nb)

        for term in self._iterRows('particle'):
            sigma, epsilon =  self.c.execute('''
                SELECT sigma, epsilon
                FROM nonbonded_param
                WHERE id=%s''' % term['nbtype']).fetchone()
            nb.addParticle(term['charge'], sigma*angstroms, epsilon*kilocalories_per_mole)


        # Bond graph (for debugging)
        g = nx.from_edgelist(self.c.execute('SELECT p0, p1 from stretch_harm_term').fetchall())
        pair_12_6_es_terms = set()
        nbnames = {1: '1-2', 2:'1-3', 3:'1-4'}

        for term in self._iterRows('pair_12_6_es_term'):
            l = nx.algorithms.shortest_path_length(g, term['p0'], term['p1'])
            pair_12_6_es_terms.add((term['p0'], term['p1']))
            print 'Scaling interaction for a %d-%d (%s) interaction' % (term['p0'], term['p1'], nbnames[l])

            a_ij, b_ij, q_ij, memo = self.c.execute('''
                SELECT aij, bij, qij, type
                FROM pair_12_6_es_param
                WHERE id=%s''' % term['param']).fetchone()
            a_ij = (a_ij*kilocalorie_per_mole*(angstroms**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
            b_ij = (b_ij*kilocalorie_per_mole*(angstroms**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
            q_ij = q_ij*elementary_charge**2

            #atom0params = nb.getParticleParameters(term['p0'])
            #atom1params = nb.getParticleParameters(term['p1'])
            #sigma_ij_orig = 0.5*(atom0params[1] + atom1params[1])
            #epsilon_ij_orig = (atom0params[2] * atom1params[2]).sqrt()
            #q_ij_orig = atom0params[0] * atom1params[0]

            new_epsilon =  b_ij**2/(4*a_ij)
            new_sigma = (a_ij / b_ij)**(1.0/6.0)

            nb.addException(term['p0'], term['p1'], q_ij, new_sigma, new_epsilon)

            if self.c.execute('SELECT COUNT(*) FROM exclusion WHERE p0=%d AND p1=%d' % (term['p0'], term['p1'])).fetchone()[0] != 1:
                raise ValueError('I can only support pair_12_6_es_terms that correspond to an exclusion')



        for term in self._iterRows('exclusion'):
            if (term['p0'], term['p1']) in pair_12_6_es_terms:
                # Desmond puts scaled 1-4 interactions in the pair_12_6_es
                # table, and then adds a corresponding exception here. We are
                # using the exception part of NonbondedForce, so we're just
                # adding the 1-4 interaction as an exception when its
                # registered, and then NOT registering it as an exception here.
                continue

            print 'Creating exception for a %d-%d (%s) interaction' % (term['p0'], term['p1'], nbnames[l])
            nb.addException(term['p0'], term['p1'], 0.0, 1.0, 0.0)



    def getPositions(self):
        positions = []
        for term in self._iterRows('particle'):
            positions.append(mm.Vec3(term['x'], term['y'], term['z'])*angstrom)
        return positions

    def _iterRows(self, table_name):
        for row in self.c.execute('SELECT * FROM %s' % table_name):
            yield dict(zip(self._tables[table_name]['names'], row))

    def _hasTable(self, table_name):
        return table_name in self._tables

    def close(self):
        if self._open:
            self.c.close()

    def __del__(self):
        self.close()


if __name__ == '__main__':
    import os
    import numpy as np

    os.system('viparr ala2.pdb ala2.dms -f amber99SB-ILDN --without-constraints')
    #os.system('viparr ala2.pdb ala2.dms -f amber96 --without-constraints')
    dms = DesmondDMSFile('ala2.dms')
    system1 = dms.createSystem()
    context1 = mm.Context(system1, mm.VerletIntegrator(0))
    context1.setPositions(dms.getPositions())
    state1 = context1.getState(getForces=True, getEnergy=True)
    force1 = state1.getForces(asNumpy=True)
    print 'DMS'
    print force1
    print 'DMS Energy', state1.getPotentialEnergy()

    pdb = app.PDBFile('ala2.pdb')
    forcefield = app.ForceField('amber99sbildn.xml')
    #forcefield = app.ForceField('amber96.xml')
    system_99 = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)
    system2 = mm.System()
    for i in range(system_99.getNumParticles()):
        system2.addParticle(system_99.getParticleMass(i))



    for i in range(system_99.getNumForces()):
        force = system_99.getForce(i)
        #print force
        if isinstance(force, mm.HarmonicBondForce):
            system2.addForce(mm.HarmonicBondForce(force))
        if isinstance(force, mm.HarmonicAngleForce):
            system2.addForce(mm.HarmonicAngleForce(force))
        if isinstance(force, mm.PeriodicTorsionForce):
            system2.addForce(mm.PeriodicTorsionForce(force))
        if isinstance(force, mm.NonbondedForce):
            system2.addForce(mm.NonbondedForce(force))

    context2 = mm.Context(system2, mm.VerletIntegrator(0))
    context2.setPositions(pdb.getPositions())
    state2 = context2.getState(getForces=True, getEnergy=True)
    force2 = state2.getForces(asNumpy=True)
    print 'PDB ENergy', state2.getPotentialEnergy()


    #print 'PDB'
    #print np.array([e.value_in_unit(nanometers) for e in pdb.getPositions()])
    #print 'DMS'
    #print np.array([e.value_in_unit(nanometers) for e in dms.getPositions()])
    #print '\n'

    print 'PDB'
    print force2


    diff = np.array([np.linalg.norm(e) for e in force1-force2])
    print diff
    print 'Atoms Where There is an Error', np.where(diff > 0.1)
