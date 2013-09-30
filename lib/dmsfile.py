'''
dmsfile.py: Load Desmond dms files

Portions copyright (c) 2013 Stanford University and the Authors
Authors: Robert McGibbon
Contributors:


Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import math

from simtk import openmm as mm
from simtk.openmm.app import forcefield as ff
from simtk.openmm.app import element, Topology, PDBFile
from simtk.unit import (Quantity, nanometer, angstrom, dalton,
                        kilocalorie_per_mole, kilojoule_per_mole,
                        radian, degree, elementary_charge)

try:
    import networkx as nx
except ImportError:
    # nx is only used when createSystem is called with verbose=True
    pass


class DesmondDMSFile(object):
    '''DesmondDMSFile parses a Desmond DMS (desmond molecular system) and
    constructs a topology and (optionally) an OpenMM System from it
    '''

    def __init__(self, file):
        '''Load a DMS file

        Parameters:
        - file (string) the name of the file to load
        '''

        # sqlite3 is included in the standard lib, but at python
        # compile time, you can disable support (I think), so it's
        # not *guarenteed* to be available. Doing the import here
        # means we only raise an ImportError if people try to use
        # this class, so the module can be safely imported
        import sqlite3

        self._open = False
        self._tables = None
        self._conn = sqlite3.connect(file)
        self._open = True
        self._readSchemas()

        if 'nbtype' not in self._tables['particle']:
            raise ValueError('No nonbonded parameters associated with this '
                             'DMS file. You can add a forcefield with the '
                             'viparr command line tool distributed with desmond')

        # Build the topology
        self.topology, self.positions = self._createTopology()
        self._topologyAtoms = list(self.topology.atoms())
        self._atomBonds = [{} for x in range(len(self._topologyAtoms))]

    def getPositions(self):
        '''Get the positions of each atom in the system
        '''
        return self.positions

    def getTopology(self):
        '''Get the topology of the system
        '''
        return self.topology

    def _createTopology(self):
        '''Build the topology of the system
        '''
        top = Topology()
        positions = []

        boxVectors = []
        for x, y, z in self._conn.execute('SELECT x, y, z FROM global_cell'):
            boxVectors.append(mm.Vec3(x, y, z)*angstrom)
        unitCellDimensions = [boxVectors[0][0], boxVectors[1][1], boxVectors[2][2]]
        top.setUnitCellDimensions(unitCellDimensions)

        atoms = {}
        lastChain = None
        lastResId = None
        c = top.addChain()
        q = '''SELECT id, name, anum, resname, resid, chain, x, y, z
        FROM particle'''
        for (atomId, atomName, atomNumber, resName, resId, chain, x, y, z) in self._conn.execute(q):
            if chain != lastChain:
                lastChain = chain
                c = top.addChain()
            if resId != lastResId:
                lastResId = resId
                if resName in PDBFile._residueNameReplacements:
                    resName = PDBFile._residueNameReplacements[resName]
                r = top.addResidue(resName, c)
                if resName in PDBFile._atomNameReplacements:
                    atomReplacements = PDBFile._atomNameReplacements[resName]
                else:
                    atomReplacements = {}

            if atomName in atomReplacements:
                atomName = atomReplacements[atomName]

            elem = element.get_by_atomic_number(atomNumber)
            atoms[atomId] = top.addAtom(atomName, elem, r)
            positions.append(mm.Vec3(x, y, z)*angstrom)

        for p0, p1 in self._conn.execute('SELECT p0, p1 FROM bond'):
            top.addBond(atoms[p0], atoms[p1])

        return top, positions

    def createSystem(self, nonbondedMethod=ff.NoCutoff, nonbondedCutoff=Quantity(value=1.0, unit=nanometer),
                     ewaldErrorTolerance=0.0005, removeCMMotion=True, hydrogenMass=None, verbose=False):
        '''Construct an OpenMM System representing the topology described by this dms file

        Parameters:
        - nonbondedMethod (object=NoCutoff) The method to use for nonbonded interactions.  Allowed values are
          NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
        - nonbondedCutoff (distance=1*nanometer) The cutoff distance to use for nonbonded interactions
        - ewaldErrorTolerance (float=0.0005) The error tolerance to use if nonbondedMethod is Ewald or PME.
        - removeCMMotion (boolean=True) If true, a CMMotionRemover will be added to the System
        - hydrogenMass (mass=None) The mass to use for hydrogen atoms bound to heavy atoms.  Any mass added to a hydrogen is
          subtracted from the heavy atom to keep their total mass the same.
        '''
        self._checkForUnsupportedTerms()
        sys = mm.System()

        # Buld the box dimensions
        sys = mm.System()
        boxSize = self.topology.getUnitCellDimensions()
        if boxSize is not None:
            sys.setDefaultPeriodicBoxVectors((boxSize[0], 0, 0), (0, boxSize[1], 0), (0, 0, boxSize[2]))
        elif nonbondedMethod in (ff.CutoffPeriodic, ff.Ewald, ff.PME):
            raise ValueError('Illegal nonbonded method for a non-periodic system')

        # Create all of the particles
        for mass in self._conn.execute('SELECT mass from particle'):
            sys.addParticle(mass[0]*dalton)

        # Add all of the forces
        bonds = self._addBondsToSystem(sys)
        angles = self._addAnglesToSystem(sys)
        perioic = self._addPeriodicTorsionsToSystem(sys)
        improper = self._addImproperHarmonicTorsionsToSystem(sys)
        cmap = self._addCMAPToSystem(sys)
        nb = self._addNonbondedForceToSystem(sys, verbose)

        # Finish configuring the NonbondedForce.
        methodMap = {ff.NoCutoff:mm.NonbondedForce.NoCutoff,
                     ff.CutoffNonPeriodic:mm.NonbondedForce.CutoffNonPeriodic,
                     ff.CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic,
                     ff.Ewald:mm.NonbondedForce.Ewald,
                     ff.PME:mm.NonbondedForce.PME}
        nb.setNonbondedMethod(methodMap[nonbondedMethod])
        nb.setCutoffDistance(nonbondedCutoff)
        nb.setEwaldErrorTolerance(ewaldErrorTolerance)

        # Adjust masses.
        if hydrogenMass is not None:
            for atom1, atom2 in self.topology.bonds():
                if atom1.element == element.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == element.hydrogen and atom1.element not in (element.hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)

        # Add a CMMotionRemover.
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        return sys

    def _addBondsToSystem(self, sys):
        '''Create the harmonic bonds
        '''
        bonds = mm.HarmonicBondForce()
        sys.addForce(bonds)

        q = '''SELECT p0, p1, r0, fc, constrained
        FROM stretch_harm_term INNER JOIN stretch_harm_param
        ON stretch_harm_term.param=stretch_harm_param.id'''
        for p0, p1, r0, fc, constrained in self._conn.execute(q):
            if constrained:
                sys.addConstraint(p0, p1, r0*angstrom)
            else:
                # Desmond writes the harmonic bond force without 1/2
                # so we need to to double the force constant
                bonds.addBond(p0, p1, r0*angstrom, 2*fc*kilocalorie_per_mole/angstrom**2)

            # Record information that will be needed for constraining angles.
            self._atomBonds[p0][p1] = r0*angstrom
            self._atomBonds[p1][p0] = r0*angstrom

        return bonds

    def _addAnglesToSystem(self, sys):
        '''Create the harmonic angles
        '''
        angles = mm.HarmonicAngleForce()
        sys.addForce(angles)
        degToRad = math.pi/180

        q = '''SELECT p0, p1, p2, theta0, fc, constrained
        FROM angle_harm_term INNER JOIN angle_harm_param
        ON angle_harm_term.param=angle_harm_param.id'''
        for p0, p1, p2, theta0, fc, constrained in self._conn.execute(q):
            if constrained:
                l1 = self._atomBonds[p1][p0]
                l2 = self._atomBonds[p1][p2]
                length = (l1*l1 + l2*l2 - 2*l1*l2*math.cos(theta0*degToRad)).sqrt()
                sys.addConstraint(p0, p2, length)
            else:
                # Desmond writes the harmonic angle force without 1/2
                # so we need to to double the force constant
                angles.addAngle(p0, p1, p2, theta0*degToRad, 2*fc*kilocalorie_per_mole/radian**2)

        return angles

    def _addPeriodicTorsionsToSystem(self, sys):
        '''Create the torsion terms
        '''
        periodic = mm.PeriodicTorsionForce()
        sys.addForce(periodic)

        q = '''SELECT p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6
        FROM dihedral_trig_term INNER JOIN dihedral_trig_param
        ON dihedral_trig_term.param=dihedral_trig_param.id'''
        for p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 in self._conn.execute(q):
            for order, fc in enumerate([fc0, fc1, fc2, fc3, fc4, fc5, fc6]):
                if fc == 0:
                    continue
                periodic.addTorsion(p0, p1, p2, p3, order, phi0*degree, fc*kilocalorie_per_mole)


    def _addImproperHarmonicTorsionsToSystem(self, sys):
        '''Create the improper harmonic torsion terms
        '''
        if not self._hasTable('improper_harm_term'):
            return

        harmonicTorsion = mm.CustomTorsionForce('k*(theta-theta0)^2')
        harmonicTorsion.addPerTorsionParameter('theta0')
        harmonicTorsion.addPerTorsionParameter('k')
        sys.addForce(harmonicTorsion)

        q = '''SELECT p0, p1, p2, p3, phi0, fc
        FROM improper_harm_term INNER JOIN improper_harm_param
        ON improper_harm_term.param=improper_harm_param.id'''
        for p0, p1, p2, p3, phi0, fc in self._conn.execute(q):
            harmonicTorsion.addTorsion(p0, p1, p2, p3, [phi0*degree, fc*kilocalorie_per_mole])

    def _addCMAPToSystem(self, sys):
        '''Create the CMAP terms
        '''
        if not self._hasTable('torsiontorsion_cmap_term'):
            return

        # Create CMAP torsion terms
        cmap = mm.CMAPTorsionForce()
        sys.addForce(cmap)
        cmap_indices = {}

        for name in [k for k in self._tables.keys() if k.startswith('cmap')]:
            size2 = self._conn.execute('SELECT COUNT(*) FROM %s' % name).fetchone()[0]
            fsize = math.sqrt(size2)
            if fsize != int(fsize):
                raise ValueError('Non-square CMAPs are not supported')
            size = int(fsize)

            map = [0 for i in range(size2)]
            for phi, psi, energy in self._conn.execute("SELECT phi, psi, energy FROM %s" % name):
                i = int((phi % 360) / (360.0 / size))
                j = int((psi % 360) / (360.0 / size))
                map[i+size*j] = energy
            index = cmap.addMap(size, map*kilocalorie_per_mole)
            cmap_indices[name] = index

        q = '''SELECT p0, p1, p2, p3, p4, p5, p6, p7, cmapid
        FROM torsiontorsion_cmap_term INNER JOIN torsiontorsion_cmap_param
        ON torsiontorsion_cmap_term.param=torsiontorsion_cmap_param.id'''
        for p0, p1, p2, p3, p4, p5, p6, p7, cmapid in self._conn.execute(q):
            cmap.addTorsion(cmap_indices[cmapid], p0, p1, p2, p3, p4, p5, p6, p7)

    def _addNonbondedForceToSystem(self, sys, verbose):
        '''Create the nonbonded force
        '''
        nb = mm.NonbondedForce()
        sys.addForce(nb)

        q = '''SELECT charge, sigma, epsilon
        FROM particle INNER JOIN nonbonded_param
        ON particle.nbtype=nonbonded_param.id'''
        for charge, sigma, epsilon in self._conn.execute(q):
            nb.addParticle(charge, sigma*angstrom, epsilon*kilocalorie_per_mole)

        if verbose:
            # Bond graph (for debugging)
            g = nx.from_edgelist(self._conn.execute('SELECT p0, p1 FROM stretch_harm_term').fetchall())
            nbnames = {1: '1-2', 2:'1-3', 3:'1-4'}

        q = '''SELECT p0, p1, aij, bij, qij
        FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param
        ON pair_12_6_es_term.param=pair_12_6_es_param.id;'''
        for p0, p1, a_ij, b_ij, q_ij in self._conn.execute(q):
            if verbose:
                l = nx.algorithms.shortest_path_length(g, p0, p1)
                print 'Scaling interaction for a %d-%d (%s) interaction' % (p0, p1, nbnames[l])
            a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
            b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
            q_ij = q_ij*elementary_charge**2

            if (b_ij._value == 0.0) or (a_ij._value == 0.0):
                new_epsilon = 0
                new_sigma = 1
            else:
                new_epsilon =  b_ij**2/(4*a_ij)
                new_sigma = (a_ij / b_ij)**(1.0/6.0)
            nb.addException(p0, p1, q_ij, new_sigma, new_epsilon)

        n_total = self._conn.execute('''SELECT COUNT(*) FROM pair_12_6_es_term''').fetchone()
        n_in_exclusions= self._conn.execute('''SELECT COUNT(*)
        FROM exclusion INNER JOIN pair_12_6_es_term
        ON exclusion.p0==pair_12_6_es_term.p0 AND exclusion.p1==pair_12_6_es_term.p1''').fetchone()
        if not n_total == n_in_exclusions:
            raise NotImplementedError('All pair_12_6_es_terms must have a corresponding exclusion')

        # Desmond puts scaled 1-4 interactions in the pair_12_6_es
        # table, and then adds a corresponding exception here. We are
        # using the exception part of NonbondedForce, so we're just
        # adding the 1-4 interaction as an exception when its
        # registered, and then NOT registering it as an exception here.
        q = '''SELECT E.p0, E.p1
        FROM exclusion E LEFT OUTER JOIN pair_12_6_es_term P ON
        E.p0 = P.p0 and E.p1 = P.p1
        WHERE P.p0 is NULL'''
        # http://stackoverflow.com/questions/5464131/finding-pairs-that-do-not-exist-in-a-different-table
        for p0, p1 in self._conn.execute(q):
            if verbose:
                l = nx.algorithms.shortest_path_length(g, p0, p1)
                print 'Creating exception for a %d-%d (%s) interaction' % (p0, p1, nbnames[l])
            nb.addException(p0, p1, 0.0, 1.0, 0.0)

        return nb

    def _hasTable(self, table_name):
        '''Does our DMS file contain this table?
        '''
        return table_name in self._tables

    def _readSchemas(self):
        '''Read the schemas of each of the tables in the dms file, populating
        the `_tables` instance attribute
        '''
        tables = {}
        for table in self._conn.execute("SELECT name FROM sqlite_master WHERE type='table'"):
            names = []
            for e in self._conn.execute('PRAGMA table_info(%s)' % table):
                names.append(str(e[1]))
            tables[str(table[0])] = names
        self._tables = tables

    def _checkForUnsupportedTerms(self):
        '''Check the file for forcefield terms that are not currenty supported,
        raising a NotImplementedError
        '''
        if 'posre_harm_term' in self._tables:
            raise NotImplementedError('Position restraints are not implemented.')
        flat_bottom_potential_terms = ['stretch_fbhw_term', 'angle_fbhw_term',
                                       'improper_fbhw_term', 'posre_fbhw_term']
        if any((t in self._tables) for t in flat_bottom_potential_terms):
            raise NotImplementedError('Flat bottom potential terms '
                                      'are not implemeneted')

    def close(self):
        '''Close the SQL connection
        '''
        if self._open:
            self._conn.close()

    def __del__(self):
        self.close()
