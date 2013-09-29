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

# Plots

mm-gmx.png: mmforces.dat gmxforces.dat
	bin/plot-max-deltaf mmforces.dat gmxforces.dat mm-gmx.png

mm-des.png: mmforces.dat desforces.dat
	bin/plot-max-deltaf mmforces.dat desforces.dat mm-des.png

des-gmx.png: desforces.dat gmxforces.dat
	bin/plot-max-deltaf desforces.dat gmxforces.dat des-gmx.png
