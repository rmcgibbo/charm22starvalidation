all: mm-gmx.png mm-des.png des-gmx.png

clean:
	rm -rf *.gro *.top \#topol.top* *.itp *.dms dessert_* eneseq forces.dtr trajectory.dtr *.png topology-match.pdb
	rm -rf *.dat rm \#conf.*

# Everything starts from the gro file, which we take as
# the "reference geometry" of the system

conf.gro topol.top: 1LYD.pdb
	pdb2gmx -f 1LYD.pdb -ff charmm22star -water none -ignh
	editconf -f conf.gro -o conf.gro -box 10 10 10
	rm posre.itp

# Force calculation
mmforces.dat: conf.gro topol.top
	bin/mm-forces conf.gro topol.top mmforces.dat
gmxforces.dat: conf.gro topol.top
	bin/gmx-forces conf.gro topol.top gmxforces.dat
desforces.dat: conf.gro bin/desmond-forces
	bin/desmond-forces conf.gro desforces.dat
mm2forces.dat: conf.gro
	bin/mm2-forces conf.gro mm2forces.dat


# PDB files with bfactors giving the force error
mm-gmx.pdb: conf.gro mmforces.dat gmxforces.dat
	editconf -f conf.gro -o mm-gmx.pdb
	bin/replace_bfactor mm-gmx.pdb mmforces.dat gmxforces.dat
mm-des.pdb: conf.gro mmforces.dat desforces.dat
	editconf -f conf.gro -o mm-des.pdb
	bin/replace_bfactor mm-des.pdb mmforces.dat desforces.dat
des-gmx.pdb: conf.gro desforces.dat gmxforces.dat
	editconf -f conf.gro -o des-gmx.pdb
	bin/replace_bfactor des-gmx.pdb desforces.dat gmxforces.dat
mm2-des.pdb: conf.gro mm2forces.dat desforces.dat
	editconf -f conf.gro -o mm2-des.pdb
	bin/replace_bfactor mm2-des.pdb mm2forces.dat desforces.dat


# Plots
mm-gmx.png: mmforces.dat gmxforces.dat
	bin/plot-max-deltaf mmforces.dat gmxforces.dat mm-gmx.png
mm-des.png: mmforces.dat desforces.dat
	bin/plot-max-deltaf mmforces.dat desforces.dat mm-des.png
des-gmx.png: desforces.dat gmxforces.dat
	bin/plot-max-deltaf desforces.dat gmxforces.dat des-gmx.png
mm2-des.png: mm2forces.dat desforces.dat
	bin/plot-max-deltaf desforces.dat mm2forces.dat mm2-des.png
