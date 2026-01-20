# Tutorial: VASP Calculations

This tutorial covers generating VASP input files and parsing output results.

## Learning Objectives

- Generate VASP input files (INCAR, KPOINTS, POSCAR, POTCAR)
- Configure Hubbard U and vdW corrections
- Parse calculation results

## Setup

```python
from ase.io import read
from nh3sofc.calculators.vasp import VASPInputGenerator

atoms = read("structure.xyz")
```

## Generating Input Files

### Basic Relaxation

```python
vasp = VASPInputGenerator(
    atoms,
    calc_type="relax",
    work_dir="./relax"
)

files = vasp.generate_all(
    encut=520,
    kspacing=0.03,
    ediff=1e-6,
    ediffg=-0.02,
)

print("Generated files:", list(files.keys()))
```

### With Hubbard U Correction

```python
vasp = VASPInputGenerator(atoms, calc_type="relax", work_dir="./relax_U")

files = vasp.generate_all(
    encut=520,
    kspacing=0.03,
    hubbard_u={"V": 3.25, "Ti": 3.0},  # Element: U value (eV)
)
```

### With vdW-D3 Correction

```python
files = vasp.generate_all(
    encut=520,
    kspacing=0.03,
    vdw="D3BJ",  # Options: "D3", "D3BJ", "dDsC"
)
```

## Calculation Types

### Static (Single-Point)

```python
vasp = VASPInputGenerator(atoms, calc_type="static", work_dir="./static")
vasp.generate_all(encut=520)
```

### Frequency Calculation

```python
from nh3sofc.calculators.vasp import FrequencyCalculation

freq = FrequencyCalculation(
    atoms,
    work_dir="./freq",
    selective_dynamics=[48, 49, 50, 51]  # Only vibrate adsorbate
)

freq.setup(nfree=2, encut=520)
```

### NEB Calculation

```python
vasp = VASPInputGenerator(
    images,  # List of NEB images
    calc_type="neb",
    work_dir="./neb"
)

vasp.generate_all(
    encut=520,
    images=5,
    lclimb=True,  # Climbing image NEB
)
```

## Parsing Results

### Parse OUTCAR

```python
from nh3sofc.calculators.vasp import VASPOutputParser

parser = VASPOutputParser("./relax")

results = parser.parse_outcar()
print(f"Energy: {results['energy']:.4f} eV")
print(f"Converged: {results['converged']}")
print(f"Ionic steps: {results['ionic_steps']}")
```

### Get Trajectory

```python
trajectory = parser.get_trajectory()
print(f"Number of steps: {len(trajectory)}")

# Save final structure
final = trajectory[-1]
final.write("relaxed.xyz")
```

### Check Convergence

```python
status = parser.check_convergence()

if status["electronic"]:
    print("Electronic converged")
if status["ionic"]:
    print("Ionic converged")
if status["complete"]:
    print("Calculation complete")
```

## PBS Job Script

```python
from nh3sofc.jobs import PBSScriptGenerator

pbs = PBSScriptGenerator(
    job_name="vasp_relax",
    nodes=1,
    ppn=24,
    walltime="24:00:00",
    queue="normal"
)

pbs.generate(work_dir="./relax")
```

## Complete Example

```python
from ase.io import read
from nh3sofc.calculators.vasp import VASPInputGenerator, VASPOutputParser
from nh3sofc.jobs import PBSScriptGenerator

# 1. Generate inputs
atoms = read("nh3_on_surface.xyz")

vasp = VASPInputGenerator(atoms, calc_type="relax", work_dir="./calc")
vasp.generate_all(
    encut=520,
    kspacing=0.03,
    hubbard_u={"V": 3.25},
    vdw="D3BJ",
    ispin=2,
)

# 2. Generate PBS script
pbs = PBSScriptGenerator(nodes=1, ppn=24, walltime="24:00:00")
pbs.generate("./calc")

# 3. After calculation completes, parse results
parser = VASPOutputParser("./calc")
results = parser.parse_outcar()
print(f"Final energy: {results['energy']:.4f} eV")
```

## Next Steps

- [Decomposition Pathway](decomposition.md) - Set up decomposition calculations
- [Thermochemistry](thermochemistry.md) - Add frequency calculations
