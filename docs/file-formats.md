# File Formats

## Mapping Definition File
An example of mapping definition file for the monosaccharide allose taken from `test/data/sugar.map` is shown below.

Molecule names (as present in the gro coordinate file) are used as section headers inside square brackets.
Each of the following lines describes a single coarse-grained bead mapping.
The items on a line are:
- the name of the bead
- its MARTINI bead type
- optionally the bead charge
- a list of all the atom names it should contain

All items on a line are whitespace separated.
Multiple molecules may be specified in their own sections.
It is not recommended to provide a mapping for water since MARTINI water combines four molecules into a single bead which is not yet supported by PyCGTOOL.
Note that bead charges in the MARTINI framework are by convention integers and are used only for formally charged functional groups.  An example of a molecule mapping using charges can be found in `test/data/dppc.map`.

```
; comments begin with a semicolon
[ALLA]
C1 P3 C1 O1
C2 P3 C2 O2
C3 P3 C3 O3
C4 P3 C4 O4
C5 P2 C5 C6 O6
O5 P4 O5

[SOL]
W P4 OW HW1 HW2
```

The Martini 3.0 force field has now introduced the use of virtual sites for polycylic compounds or polysaccharides in order to improve numerical stability for highly constrained structures.
Below is an example of how virtual sites can be included in the mapping file for naphthalene (see `/test/data/martini3/naphthalene.map`).
Similar to the previously described mapping syntax, a virtual site is defined with a prefix of ``@v`` as follows `@v [name] [type] [charge] [constructing beads]`.
The constructing beads refer to a list of space delimited coarse grained bead names from which the position of the virtual site will be calculated.
Currently virtual sites can be constructed from either the center of geometry or mass of the constructing sites via the `--virtual_map_center` flag.

```
[ NAPH ]
R1 TC4 C8 H8 C9 H9
R2 TC4 C6 H6 C7 H7
@v R3 TC4 R1 R2 R4 R5
R4 TC4 C1 H1 C2 H2
R5 TC4 C3 H3 C4 H4
```

## Bond Definition File
An example bond definition file for the monosaccharide allose taken from `test/data/sugar.bnd` is shown below.

As in the mapping definition file, molecule names are used as section headers inside square brackets. 
The following lines define bonds lengths, angles and dihedrals between coarse-grained beads. 
Each line is a list of bead names, using the names defined in the mapping file. 
Two bead names on a line defines a bond length, three defines an angle, and four defines a dihedral.

If no angles are defined for a molecule, PyCGTOOL will construct all angles from the list of bonds. 
This may also be enabled for dihedrals via the `--generate_dihedrals` flag, but is not recommended as in most cases coarse-grained models do not require dihedrals. 
Additionally, any angles inside a triangle of bond lengths are excluded from the output as they often cause simulation stability issues when used in conjunction with LINCS.

```
; comments begin with a semicolon
[ALLA]
C1 C2
C2 C3
C3 C4
C4 C5
C5 O5
O5 C1

C1 C2 C3
C2 C3 C4
C3 C4 C5
C4 C5 O5
C5 O5 C1
O5 C1 C2

C1 C2 C3 C4
C2 C3 C4 C5
C3 C4 C5 O5
C4 C5 O5 C1
C5 O5 C1 C2
O5 C1 C2 C3
```
