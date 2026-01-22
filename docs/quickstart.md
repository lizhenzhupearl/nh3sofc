# Quick Start Guide

This guide will walk you through the basic workflow for studying NH3 decomposition on a perovskite oxynitride surface.

## Overview

A typical calculation workflow consists of:

1. **Build Structure** - Create surface from bulk, add defects
2. **Place Adsorbate** - Add NH3 or decomposition intermediates
3. **Set Up Calculation** - Generate VASP inputs
4. **Run & Analyze** - Parse results, calculate properties

## Step 1: Build an Oxynitride Surface

### Load Bulk Structure

```python
from nh3sofc.structure import BulkStructure

# From CIF file
bulk = BulkStructure.from_cif("LaVO3.cif")

# Or from ASE/pymatgen
from ase.build import bulk as ase_bulk
bulk_atoms = ase_bulk('Ni', 'fcc', a=3.52)
bulk = BulkStructure(bulk_atoms)

# Check structure
print(f"Formula: {bulk.atoms.get_chemical_formula()}")
print(f"Space group: {bulk.get_spacegroup()}")
```

### Create Surface Slab

```python
from nh3sofc.structure import SurfaceBuilder

builder = SurfaceBuilder(bulk)

# Create (001) surface with 6 layers and 15 Å vacuum
surface = builder.create_surface(
    miller_index=(0, 0, 1),
    layers=6,
    vacuum=15.0,
    fix_bottom=2  # Fix bottom 2 layers
)

print(f"Slab atoms: {len(surface.atoms)}")
print(f"Cell: {surface.atoms.cell.lengths()}")
```

### Create Oxynitride Defects

```python
from nh3sofc.structure import DefectBuilder

defect_builder = DefectBuilder(surface.atoms)

# Create LaVON2-x with 10% vacancies
oxynitride = defect_builder.create_oxynitride(
    nitrogen_fraction=0.67,      # Replace 2/3 of O with N
    vacancy_concentration=0.10   # 10% O/N vacancies
)

print(f"Composition: {oxynitride.get_chemical_formula()}")
```

## Step 2: Place NH3 Adsorbate

```python
from nh3sofc.structure import AdsorbatePlacer

placer = AdsorbatePlacer(oxynitride)

# Method 1: On-top of specific atom types
configs = placer.add_on_site(
    adsorbate="NH3",
    site_type="ontop",        # Site type: "ontop", "bridge", or "hollow"
    atom_types=["La", "V"],   # Only on La or V atoms
    height=2.0
)
print(f"Generated {len(configs)} configurations")

# Method 2: Random placement with collision detection
random_configs = placer.add_with_collision(
    adsorbate="NH3",
    n_configs=10,             # Number of configurations to generate
    min_distance=2.0,
    height=2.0
)

# Method 3: All unique sites (using CatKit if available)
try:
    all_sites = placer.add_catkit("NH3", site_type="ontop")
except ImportError:
    print("CatKit not installed, using collision-aware placement")
    all_sites = random_configs if random_configs else []
```

## Step 3: Set Up VASP Calculation

### Generate Input Files

```python
from nh3sofc.calculators.vasp import VASPInputGenerator

# Select first configuration
atoms = configs[0]

# Create generator
vasp = VASPInputGenerator(
    atoms,
    calc_type="relax",  # Options: relax, static, neb, frequency, md
    work_dir="./calc/LaVON_001_NH3"
)

# Generate all input files
files = vasp.generate_all(
    encut=520,              # Plane-wave cutoff
    kspacing=0.03,          # K-point spacing
    hubbard_u={"V": 3.25},  # DFT+U for V
    vdw="D3BJ",             # vdW correction
    is_surface=True         # Surface-specific settings
)

print("Generated files:")
for name, path in files.items():
    print(f"  {name}: {path}")
```

### Generate PBS Script

```python
from nh3sofc.jobs import create_vasp_job

pbs_file = create_vasp_job(
    work_dir="./calc/LaVON_001_NH3",
    nodes=2,
    ppn=24,
    walltime="24:00:00",
    queue="normal"
)

print(f"Submit with: qsub {pbs_file}")
```

## Step 4: Parse Results

After the VASP job completes:

```python
from nh3sofc.calculators.vasp import VASPOutputParser

parser = VASPOutputParser("./calc/LaVON_001_NH3")

# Check convergence
if parser.is_converged():
    energy = parser.get_energy()
    forces = parser.get_forces()
    final_structure = parser.get_final_structure()

    print(f"Energy: {energy:.4f} eV")
    print(f"Max force: {parser.get_max_force():.4f} eV/Å")
else:
    print("Calculation did not converge!")
```

## Step 5: Calculate Adsorption Energy

```python
from nh3sofc.analysis import AdsorptionEnergyCalculator

calc = AdsorptionEnergyCalculator()

# Set reference energies (from your calculations)
calc.set_surface_energy(-250.0)  # Clean surface energy
calc.set_gas_reference("NH3", -19.54)  # Gas phase NH3

# Calculate adsorption energy
E_ads = calc.calculate(
    e_total=-270.5,  # NH3/surface energy
    adsorbate="NH3"
)

print(f"Adsorption energy: {E_ads:.3f} eV")
```

## Complete Workflow Example

Here's a complete script for studying NH3 decomposition:

```python
from nh3sofc.structure import BulkStructure, SurfaceBuilder, DefectBuilder
from nh3sofc.structure import DecompositionBuilder
from nh3sofc.workflows import DecompositionWorkflow

# 1. Build oxynitride surface
bulk = BulkStructure.from_cif("LaVO3.cif")
surface = SurfaceBuilder(bulk).create_surface((0,0,1), layers=6, vacuum=15)
oxynitride = DefectBuilder(surface.atoms).create_oxynitride(
    nitrogen_fraction=0.67,
    vacancy_concentration=0.1
)

# 2. Add NH3 and optimize (assuming you have the relaxed structure)
# nh3_on_surface = ... (load from previous calculation)

# 3. Set up decomposition workflow
workflow = DecompositionWorkflow(
    nh3_on_slab=nh3_on_surface,
    work_dir="./decomposition",
    n_configs_per_step=5,
    calculator="vasp",
    encut=520,
    hubbard_u={"V": 3.25},
)

# 4. Generate all configurations and set up calculations
workflow.setup()

# 5. After calculations complete, analyze results
results = workflow.parse_results()
energy_profile = workflow.get_energy_profile()
workflow.print_summary()
```

## Next Steps

- [Tutorial: Surface Building](tutorials/surface_building.md)
- [Tutorial: NH3 Decomposition](tutorials/decomposition.md)
- [Tutorial: NEB Calculations](tutorials/neb.md)
- [API Reference](api/index.md)
