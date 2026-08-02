# Full Arguments List

The `pycgtool --help` command also shows all arguments together with usage information.

## Input Files

`topology`
: Atomistic simulation topology

`trajectory`
: Atomistic simulation trajectory

`-m`, `--mapping`
: CG mapping definition file

`-b`, `--bondset`
: Bond definition file

`-i`, `--itp`
: Optional GROMACS ITP file containing atom masses and charges

`--begin`
: First trajectory frame to use (default: `0`)

`--end`
: Last trajectory frame to use (default: all remaining frames)

## Output Files

`--out-dir`
: Directory for output files (default: current directory)

`--output-name`
: Base name for output files (default: `out`)

`--output`
: Coordinate output format (default: `gro`)

`--output-xtc`
: Write a pseudo-CG trajectory

`--output-forcefield`
: Write a GROMACS forcefield directory

`--dump-measurements`
: Write sampled bond measurements

`--dump-n-values`
: Number of measurements to sample (default: `10000`)

## Mapping Options

`--map-center`
: Mapping method: `geom`, `mass`, or `first` (default: `geom`)

`--virtual-map-center`
: Virtual site mapping method: `geom` or `mass` (default: `geom`)

`--backmapper-resname`
: Residue name for which to train an experimental backmapper

## Bond Options

`--constr-threshold`
: Convert bonds above this force constant to constraints (default: `100000`)

`--temperature`
: Reference simulation temperature (default: `310`)

`--default-fc`
: Use default MARTINI force constants

`--generate-angles`
: Generate angles from bonds (enabled by default)

`--generate-dihedrals`
: Generate dihedrals from bonds

`--length-form`, `--angle-form`, `--dihedral-form`
: Functional forms for the corresponding bonded terms

## Run Options

`--profile`
: Write profiling data to `gprof.out`

`--log-level`
: Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR`, or `CRITICAL` (default: `INFO`)
