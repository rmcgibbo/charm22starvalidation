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
        self.c = sqlite3.connect(file)
        self._open = True

        tables = {}
        for table in self.c.execute("SELECT name FROM sqlite_master WHERE type='table'"):
            names = []
            for e in self.c.execute('PRAGMA table_info(%s)' % table):
                names.append(str(e[1]))
            tables[str(table[0])] = names
        self._tables = tables

        if 'nbtype' not in tables['particle']:
            raise ValueError('No nonbonded parameters associated with this DMS file.')

        # Build the topology
        self.topology = self._getTopology()
        self.positions = self._getPositions()
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
        
    def _getTopology(self):
        '''Build the topology of the system
        '''
        top = Topology()
        boxVectors = []
        for _, x, y, z in self.c.execute('select * FROM global_cell'):
            boxVectors.append(mm.Vec3(x, y, z)*angstrom)
        unitCellDimensions = [boxVectors[0][0], boxVectors[1][1], boxVectors[2][2]]
        top.setUnitCellDimensions(unitCellDimensions)

        atoms = {}
        lastChain = None
        lastResId = None
        c = top.addChain()
        query = '''SELECT id, name, anum, resname, resid, chain
                   FROM particle ORDER BY chain, resid'''
        for (atomId, atomName, atomNumber, resName, resId, chain) in self.c.execute(query):
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

        for p0, p1 in self.c.execute('SELECT p0, p1 from bond'):
            top.addBond(atoms[p0], atoms[p1])

        return top

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
        sys = mm.System()

        # Buld the box dimensions
        sys = mm.System()
        boxSize = self.topology.getUnitCellDimensions()
        if boxSize is not None:
            sys.setDefaultPeriodicBoxVectors((boxSize[0], 0, 0), (0, boxSize[1], 0), (0, 0, boxSize[2]))
        elif nonbondedMethod in (ff.CutoffPeriodic, ff.Ewald, ff.PME):
            raise ValueError('Illegal nonbonded method for a non-periodic system')

        # Create all of the particles
        for term in self._iterRows('particle'):
            sys.addParticle(term['mass']*dalton)

        # Add all of the forces
        bonds = self._addBondsToSystem(sys)
        angles = self._addAnglesToSystem(sys)
        perioic = self._addPeriodicTorsionsToSystem(sys)
        improper = self._addImproperTorsionsToSystem(sys)
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
        for term in self._iterRows('stretch_harm_term'):
            r0, fc = self.c.execute('SELECT r0, fc FROM stretch_harm_param WHERE id=?', (term['param'],)).fetchone()
            p0, p1 = term['p0'], term['p1']

            if term['constrained']:
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

        for term in self._iterRows('angle_harm_term'):
            theta0, fc = self.c.execute('''
                SELECT theta0, fc
                FROM angle_harm_param
                WHERE id=?''', (term['param'], )).fetchone()
            p0, p1, p2 = term['p0'], term['p1'], term['p2']

            if term['constrained']:
                l1 = self._atomBonds[p1][p0]
                l2 = self._atomBonds[p1][p2]
                length = (l1*l1 + l2*l2 - 2*l1*l2*math.cos(theta0*degToRad)).sqrt()
                sys.addConstraint(p0, p2, length)
            else:
                # Desmond writes the harmonic angle force without 1/2
                # so we need to to double the force constant
                angles.addAngle(term['p0'], term['p1'], term['p2'],
                                theta0*degToRad, 2*fc*kilocalorie_per_mole/radian**2)

        return angles

    def _addPeriodicTorsionsToSystem(self, sys):
        '''Create the torsion terms
        '''
        periodic = mm.PeriodicTorsionForce()
        sys.addForce(periodic)
        for term in self._iterRows('dihedral_trig_term'):
            phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 = self.c.execute('''
                SELECT phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6
                FROM dihedral_trig_param
                WHERE id=?''', (term['param'], )).fetchone()

            for order, fc in enumerate([fc0, fc1, fc2, fc3, fc4, fc5, fc6]):
                if fc == 0:
                    continue

                periodic.addTorsion(term['p0'], term['p1'], term['p2'], term['p3'],
                                    order, phi0*degree, fc*kilocalorie_per_mole)


    def _addImproperTorsionsToSystem(self, sys):
        '''Create the improper harmonic torsion terms
        '''
        if not self._hasTable('improper_harm_term'):
            return

        harmonicTorsion = mm.CustomTorsionForce('0.5*k*(theta-theta0)^2')
        harmonicTorsion.addPerTorsionParameter('theta0')
        harmonicTorsion.addPerTorsionParameter('k')
        sys.addForce(harmonicTorsion)
        for term in self._iterRows('improper_harm_term'):
            phi0, fc = self.c.execute('''
                SELECT phi0, fc FROM improper_harm_param
                WHERE id=?''', (term['param'],)).fetchone()
            harmonicTorsion.addTorsion(term['p0'], term['p1'], term['p2'],
                                       term['p3'], [phi0*degree, fc*kilocalorie_per_mole])

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
            size2 = self.c.execute('SELECT COUNT(*) FROM %s' % name).fetchone()[0]
            fsize = math.sqrt(size2)
            if fsize != int(fsize):
                raise ValueError('Non-square CMAPs are not supported')
            size = int(fsize)

            map = [0 for i in range(size2)]
            for term in self._iterRows(name):
                i = int((size*term['phi']/360.0 + size/2) % size)
                j = int((size*term['psi']/360.0 + size/2) % size)
                map[i+size*j] = term['energy']*kilocalorie_per_mole
            index = cmap.addMap(size, map)
            cmap_indices[name] = index


        for term in self._iterRows('torsiontorsion_cmap_term'):
            cmapid = self.c.execute('''SELECT cmapid FROM torsiontorsion_cmap_param
                                       WHERE id=?''', (term['param'],)).fetchone()[0]
            cmap.addTorsion(cmap_indices[cmapid],
                            term['p0'], term['p1'], term['p2'], term['p3'],
                            term['p4'], term['p5'], term['p6'], term['p7'])

    def _addNonbondedForceToSystem(self, sys, verbose):
        '''Create the nonbonded force
        '''
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        sys.addForce(nb)

        for term in self._iterRows('particle'):
            sigma, epsilon =  self.c.execute('''
                SELECT sigma, epsilon FROM nonbonded_param
                WHERE id=?''', (term['nbtype'],)).fetchone()
            nb.addParticle(term['charge'], sigma*angstrom, epsilon*kilocalorie_per_mole)

        if verbose:
            # Bond graph (for debugging)
            g = nx.from_edgelist(self.c.execute('SELECT p0, p1 FROM stretch_harm_term').fetchall())
            nbnames = {1: '1-2', 2:'1-3', 3:'1-4'}

        # Record the terms that are set by viparr in the pair_12_6_es column,
        # so that we can do the exceptions correctly
        pair_12_6_es_terms = set()
        for term in self._iterRows('pair_12_6_es_term'):
            p0, p1 = term['p0'], term['p1']
            if verbose:
                l = nx.algorithms.shortest_path_length(g, p0, p1)
                print 'Scaling interaction for a %d-%d (%s) interaction' % (p0, p1, nbnames[l])

            a_ij, b_ij, q_ij = self.c.execute('''
                SELECT aij, bij, qij FROM pair_12_6_es_param
                WHERE id=?''', (term['param'],)).fetchone()
            a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
            b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
            q_ij = q_ij*elementary_charge**2
            new_epsilon =  b_ij**2/(4*a_ij)
            new_sigma = (a_ij / b_ij)**(1.0/6.0)

            pair_12_6_es_terms.add((p0, p1))
            nb.addException(p0, p1, q_ij, new_sigma, new_epsilon)

            if self.c.execute('SELECT COUNT(*) FROM exclusion WHERE p0=? AND p1=?', (p0, p1)).fetchone()[0] != 1:
                raise ValueError('I can only support pair_12_6_es_terms that correspond to an exclusion')

        for term in self._iterRows('exclusion'):
            p0, p1 = term['p0'], term['p1']
            if (p0, p1) in pair_12_6_es_terms:
                # Desmond puts scaled 1-4 interactions in the pair_12_6_es
                # table, and then adds a corresponding exception here. We are
                # using the exception part of NonbondedForce, so we're just
                # adding the 1-4 interaction as an exception when its
                # registered, and then NOT registering it as an exception here.
                continue

            if verbose:
                print 'Creating exception for a %d-%d (%s) interaction' % (p0, p1, nbnames[l])
            nb.addException(p0, p1, 0.0, 1.0, 0.0)

        return nb

    def _getPositions(self):
        '''Build the array of atom positions
        '''
        positions = []
        for term in self._iterRows('particle'):
            positions.append(mm.Vec3(term['x'], term['y'], term['z'])*angstrom)
        return positions

    def _iterRows(self, table_name):
        '''Iterate through rows of one of the SQL tables, yielding dictionaries
        mapping the name of all of the entries to the values.
        '''
        for row in self.c.execute('SELECT * FROM %s' % table_name):
            yield dict(zip(self._tables[table_name], row))

    def _hasTable(self, table_name):
        '''Does our DMS file contain this table?
        '''
        return table_name in self._tables

    def close(self):
        '''Close the SQL connection
        '''
        if self._open:
            self.c.close()

    def __del__(self):
        self.close()

