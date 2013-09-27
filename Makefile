all: force_error.png

clean:
	rm -rf *.gro *.top \#topol.top* *.itp *.dms dessert_* eneseq forces.dtr trajectory.dtr *.png topology-match.pdb

conf.gro topol.top: 1PLX.pdb
	# build the system in gromacs
	pdb2gmx -f 1PLX.pdb -ff charmm22star -water tip3p

1PLX.dms:
	# build the system in desmond. for some reason the
	# box length wasn't being set correctly, so we're doing
	# that manually with dms-set
	viparr 1PLX.pdb 1PLX.dms -f charmm22star_aminoacids
	dms-set 1PLX.dms 1PLX.dms box.d=100

trajectory.dtr forces.dtr: 1PLX.dms
	desmond --include desmondsim.cfg --cfg boot.file=1PLX.dms

force_error.png: trajectory.dtr forces.dtr topol.top
	python compare_forces.py trajectory.dtr forces.dtr conf.gro topol.top
