import sqlite3
import warnings
from pprint import pprint

from simtk.unit import *
from simtk import openmm as mm


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
        
        # Create the harmonic bonds
        bonds = mm.HarmonicBondForce()
        sys.addForce(bonds)
        for term in self._iterRows('stretch_harm_term'):
            r0, fc = self.c.execute('SELECT r0, fc FROM stretch_harm_param WHERE id=%s' % term['param']).fetchone()
            if not term['constrained']:
                bonds.addBond(term['p0'], term['p1'], r0*angstroms,
                              fc*kilocalorie_per_mole/angstrom)
            else:
                sys.addConstraint(term['p0'], term['p1'], r0*angstoms)

        # Create the harmonic angles
        angles = mm.HarmonicAngleForce()
        sys.addForce(angles)
        for term in self._iterRows('angle_harm_term'):
            theta0, fc = self.c.execute('SELECT theta0, fc FROM angle_harm_param WHERE id=%s' % term['param']).fetchone()
            if not term['constrained']:
                angles.addAngle(term['p0'], term['p1'], term['p2'],
                                theta0*degrees, fc*kilocalorie_per_mole/angstrom)
            else:
                raise NotImplementedError('Angle restraints not implemeneted yet')

        # Create the torsion terms
        periodic = mm.PeriodicTorsionForce()
        sys.addForce(periodic)
        for term in self._iterRows('dihedral_trig_term'):
            phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 = self.c.execute('''
                SELECT phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6
                FROM dihedral_trig_param
                WHERE id=%s''' % term['param']).fetchone()

            for order, fc in enumerate([fc0, fc1, fc2, fc3, fc4, fc5, fc6]):
                if fc == 0.0:
                    continue
                periodic.addTorsion(term['p0'], term['p1'], term['p2'],
                                    term['p3'], order, phi0*degree,
                                    fc*kilocalorie_per_mole)

        # Create the improper harmonic torsion terms
        harmonicTorsion = mm.CustomTorsionForce('0.5*k*(theta-theta0)^2')
        harmonicTorsion.addPerTorsionParameter('theta0')
        harmonicTorsion.addPerTorsionParameter('k')
        sys.addForce(harmonicTorsion)
        for term in self._iterRows('improper_harm_term'):
            phi0, fc = self.c.execute('SELECT phi0, fc FROM improper_harm_param WHERE id=%s' % term['param']).fetchone()
            harmonicTorsion.addTorsion(term['p0'], term['p1'], term['p2'],
                                       term['p3'], [phi0*degree, fc*kilocalorie_per_mole])


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


        # Create the nonbonded force
        nb = mm.NonbondedForce()
        sys.addForce(nb)

        for term in self._iterRows('particle'):
            sigma, epsilon =  self.c.execute('''
                SELECT sigma, epsilon
                FROM nonbonded_param
                WHERE id=%s''' % term['nbtype']).fetchone()
            nb.addParticle(term['charge'], sigma*angstroms, epsilon*kilocalories_per_mole)

        exceptions = []
        for term in self._iterRows('pair_12_6_es_term'):
            aij, bij, qij, memo = self.c.execute('''
                SELECT aij, bij, qij, type
                FROM pair_12_6_es_param
                WHERE id=%s''' % term['param']).fetchone()
            atom0params = nb.getParticleParameters(term['p0'])
            atom1params = nb.getParticleParameters(term['p1'])
            sigma_ij_orig = (atom0params[1]+atom1params[1])/2
            epsilon_ij_orig = (atom0params[2]*atom1params[2]).sqrt()
            aij_orig = 4*epsilon_ij_orig*sigma_ij_orig**12
            bij_orig = -4*epsilon_ij_orig*sigma_ij_orig**6
            qij_orig = atom0params[0]*atom1params[0]

            # This interaction is supposed to be ontop of regular
            # nonbonded interactions, but the charges are always
            # equal to the normal charge product -- so just adding
            # this term onto the existing charge interactionn gives
            # us a double-counting of the charge-charge interaction, no?
            
            q_ij_new = qij*elementary_charge**2 + qij_orig
            a_ij_new = aij*kilocalorie_per_mole*angstroms**12 + aij_orig
            b_ij_new = bij*kilocalorie_per_mole*angstroms**6 + bij_orig

            #print 'Orig aij %s   New %s' % (aij_orig, a_ij_new.value_in_unit_system(md_unit_system))
            #print 'Orig bij %s   New %s' % (bij_orig, b_ij_new.value_in_unit_system(md_unit_system))
            #print 'Orig qij %s   New %s' % (qij_orig, q_ij_new.value_in_unit_system(md_unit_system))

            if abs(a_ij_new.value_in_unit_system(md_unit_system)) < 1e-10 or \
               abs(b_ij_new.value_in_unit_system(md_unit_system)) < 1e-10:
                sigma_new = 0
                epsilon_new = 0
            else:
                epsilon_new = b_ij_new**2 / (4*a_ij_new)
                try:
                    sigma_new = (-a_ij_new / b_ij_new)**(1.0/6.0)
                except:
                    sigma_new = 0
                    epsilon_new = 0
                    print term, 'ERROR'

            exceptions.append([term['p0'], term['p1'], qij_orig, sigma_new, epsilon])

        for exception in exceptions:
            nb.addException(*exception)

        for term in self._iterRows('exclusion'):
            # Some of these excpetions conflict with the ones specified in pair_12_6_es.
            # What does that mean? Should we calculate the pair_12_6_es interaction?
            nb.addException(term['p0'], term['p1'], 0.0, 0.0, 0.0, True)

        #print self.c.execute('SELECT * from nonbonded_info').fetchall()
        return sys
        


    def _iterRows(self, table_name):
        for row in self.c.execute('SELECT * FROM %s' % table_name):
            yield dict(zip(self._tables[table_name]['names'], row))

    def close(self):
        if self._open:
            self.c.close()

    def __del__(self):
        self.close()

def deg2rad(angle):
    return math.pi * angle / 180.0

if __name__ == '__main__':
    dms = DesmondDMSFile('1LYD.dms')
    dms.createSystem()

