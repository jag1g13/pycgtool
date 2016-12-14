all: pycgtool cgprepare em npt md samples rgyr

pycgtool:
	../../pycgtool.py -g ref.gro -x ref.xtc -m atenolol.map -b atenolol.bnd

cgprepare:
	gmx solvate -cp out.gro -cs ../../data/water.gro -o solv.gro -radius 0.21
	cp template.top topol.top
	echo "W 251" >> topol.top

em:
	gmx grompp -f em.mdp -c solv.gro -p topol.top -o em.tpr
	gmx mdrun -deffnm em

npt:
	gmx grompp -f npt.mdp -c em.gro -p topol.top -o npt.tpr
	gmx mdrun -deffnm npt -v

md:
	gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
	gmx mdrun -deffnm md -v

samples:
	../../pycgtool.py -g ref.gro -x ref.xtc -m atenolol.map -b atenolol.bnd --advanced
	mv {36KB,ref}_length.dat
	mv {36KB,ref}_angle.dat
	../../pycgtool.py -g md.gro -x md.xtc -b atenolol.bnd

comparesamples:
	./average_columns.py ref_length.dat 36KB_length.dat
	./average_columns.py ref_angle.dat 36KB_angle.dat

rgyr:
	echo 2 | gmx gyrate -f ref.xtc -s ref-for-rgyr.tpr -o ref-gyr.xvg
	echo 2 | gmx gyrate -f md.xtc -s md.tpr -o cg-gyr.xvg

clean:
	git clean -i
