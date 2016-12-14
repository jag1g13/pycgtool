all: pycgtool cgprepare em npt md

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

rgyr:
	gmx gyrate -f ref.xtc -s ref.tpr -o ref-gyr.xvg
	gmx gyrate -f md.xtc -s md.tpr -o cg-gyr.xvg

clean:
	git clean -i
