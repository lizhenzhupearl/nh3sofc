# Tutorial: High-Throughput Screening

This tutorial covers systematic screening of compositions, defects, and adsorption configurations.

## Learning Objectives

- Set up parameter space screening
- Manage large numbers of calculations
- Analyze screening results

## Overview

High-throughput screening automates exploration of:
- Vacancy concentrations
- Nitrogen fractions
- Dopant elements
- Adsorption sites

## Step 1: Define Parameter Space

```python
from nh3sofc.workflows import ScreeningWorkflow

workflow = ScreeningWorkflow(
    base_structure=surface,
    parameter_space={
        "vacancy_concentration": [0.0, 0.05, 0.10, 0.15, 0.20],
        "nitrogen_fraction": [0.5, 0.67, 0.75],
    },
    work_dir="./oxynitride_screening",
    encut=520,
    kspacing=0.03,
)
```

This creates 5 × 3 = 15 calculations.

## Step 2: Setup Calculations

```python
calc_dirs = workflow.setup()
print(f"Created {len(calc_dirs)} calculation directories")

# Directory structure:
# oxynitride_screening/
# ├── vac0.00_N0.50/
# ├── vac0.00_N0.67/
# ├── vac0.00_N0.75/
# ├── vac0.05_N0.50/
# ...
```

## Step 3: Monitor Progress

```python
status = workflow.get_status()

for name, state in status.items():
    print(f"{name}: {state}")

# Output:
# vac0.00_N0.50: completed
# vac0.00_N0.67: running
# vac0.05_N0.50: pending
# ...
```

## Step 4: Parse Results

```python
import pandas as pd

results = workflow.parse_results()

# Returns DataFrame with columns:
# - vacancy_concentration
# - nitrogen_fraction
# - energy
# - converged
# - ...

print(results.head())
```

## Step 5: Find Optimal

```python
optimal = workflow.find_optimal(metric="energy")
print(f"Optimal parameters: {optimal}")

# Or custom criteria
best = results.loc[results["energy"].idxmin()]
print(f"Lowest energy: {best}")
```

## Specialized Screening

### Composition Screening

```python
from nh3sofc.workflows import CompositionScreening

workflow = CompositionScreening(
    base_structure=surface,
    elements_to_vary={
        "La": ["La", "Sr", "Ba", "Ca"],
        "V": ["V", "Ti", "Nb", "Ta"],
    },
    work_dir="./dopant_screening",
)

workflow.setup()
# Creates 4 × 4 = 16 calculations
```

### Adsorbate Screening

```python
from nh3sofc.workflows import AdsorbateScreening

workflow = AdsorbateScreening(
    surface=surface,
    adsorbates=["NH3", "NH2", "NH", "N", "H"],
    sites=["ontop_La", "ontop_V", "bridge", "hollow"],
    work_dir="./adsorbate_screening",
)

paths = workflow.setup()
# Creates 5 × 4 = 20 calculations
```

### Generating Adsorbate Configurations

For fine-grained control over adsorbate placement, use `AdsorbatePlacer` directly.
All placement methods use **uniform sampling on SO(3)** for unbiased molecular
orientations, ensuring proper coverage of configurational space.

```python
from nh3sofc.structure import AdsorbatePlacer, SlabStructure

# Load surface
slab = SlabStructure.from_file("surface.vasp")
placer = AdsorbatePlacer(slab)

# Method 1: Random placement with random orientations
configs = placer.add_random(
    "NH3",
    n_configs=20,
    height=2.0,
    random_seed=42,
)
# Generates 20 configurations with uniformly sampled positions and orientations

# Method 2: Grid-based placement with multiple orientations per point
configs = placer.add_grid(
    "NH3",
    grid_size=(3, 3),    # 3x3 grid of positions
    orientations=4,       # 4 random orientations per grid point
    height=2.0,
    random_seed=42,
)
# Generates 3 × 3 × 4 = 36 configurations

# Method 3: Site-specific placement (most physically meaningful)
# Step 1: Identify surface atoms (top 20% of slab by z-coordinate)
# Step 2: Filter to only La and V atoms among those surface atoms
configs = placer.add_on_site(
    "NH3",
    atom_types=["La", "V"],  # Select La/V from surface atoms only
    n_orientations=5,         # 5 orientations per site
    height=2.0,
    z_threshold=0.2,          # Top 20% of slab = surface
    random_seed=42,
)
# Places above La and V atoms that are ON THE SURFACE
# Subsurface La/V atoms are automatically excluded

# Method 4: Collision-aware random placement
configs = placer.add_with_collision(
    "NH3",
    n_configs=10,
    min_distance=2.0,  # Minimum atom-atom distance
    height=2.0,
    random_seed=42,
)
```

#### Rotation Sampling

All methods apply uniform random rotations using proper SO(3) sampling:

```python
# The rotation uses: β = arccos(1 - 2u) for uniform spherical coverage
# This avoids clustering near the "poles" that occurs with naive sampling

# For reproducibility, always set random_seed
configs_a = placer.add_random("NH3", n_configs=10, random_seed=42)
configs_b = placer.add_random("NH3", n_configs=10, random_seed=42)
# configs_a and configs_b are identical
```

#### Filtering Unique Configurations

After generating many configurations, filter duplicates:

```python
from nh3sofc.structure.adsorbates import filter_unique_configs, save_configs

# Remove duplicates based on RMSD
unique = filter_unique_configs(configs, threshold=0.5)
print(f"Filtered {len(configs)} → {len(unique)} unique configurations")

# Save for calculations
paths = save_configs(unique, output_dir="./adsorbate_configs", format="vasp")
```

## Visualizing Results

```python
import matplotlib.pyplot as plt

results = workflow.parse_results()

# Heatmap
pivot = results.pivot(
    index="vacancy_concentration",
    columns="nitrogen_fraction",
    values="energy"
)

plt.figure(figsize=(8, 6))
plt.imshow(pivot.values, cmap='RdYlGn_r')
plt.colorbar(label='Energy (eV)')
plt.xlabel('Nitrogen Fraction')
plt.ylabel('Vacancy Concentration')
plt.savefig('screening_heatmap.png', dpi=150)
```

## Complete Example

```python
from ase.io import read
from nh3sofc.workflows import ScreeningWorkflow

# 1. Load base structure (POSCAR format)
surface = read("work/surfaces/LaVO3_001/surface.vasp", format="vasp")

# 2. Define screening with meaningful work directory
workflow = ScreeningWorkflow(
    base_structure=surface,
    parameter_space={
        "vacancy_concentration": [0.0, 0.10, 0.20],
        "nitrogen_fraction": [0.5, 0.67],
    },
    work_dir="work/screening/LaVON_composition",
    encut=520,
)

# 3. Setup (generates POSCAR files with atoms sorted by element)
workflow.setup()

# 4. After completion, analyze
results = workflow.parse_results()

# 5. Find optimal
optimal = workflow.find_optimal("energy")
print(f"Best: {optimal}")

# 6. Save results
results.to_csv("work/screening/LaVON_composition/results.csv")
```

## Best Practices

1. **Start small** - Test with subset before full screening
2. **Use database** - Store all results in ASE database
3. **Check convergence** - Verify all calculations completed
4. **Automate submission** - Use job arrays for HPC

## Next Steps

- [Surface Comparison](surface_comparison.md) - Rank screened surfaces
- [MACE Force Fields](mace.md) - Use ML for faster screening
