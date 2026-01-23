# Tutorial: Surface Building

This tutorial covers building surfaces from bulk structures, handling surface polarity, and creating oxynitride defects.

## Learning Objectives

- Load bulk structures from CIF files
- Generate surface slabs with different Miller indices
- Handle polar surfaces (critical for perovskites)
- Create symmetric slabs with specific terminations
- Create oxynitride structures with different defect placement strategies
- Generate configuration pools for screening studies

## Step 1: Load Bulk Structure

```python
from nh3sofc.structure import BulkStructure

# Load from CIF
bulk = BulkStructure.from_cif("LaVO3.cif")
print(f"Formula: {bulk.atoms.get_chemical_formula()}")
print(f"Space group: {bulk.get_spacegroup()}")

# Create supercell if needed
supercell = bulk.make_supercell((2, 2, 1))
```

## Step 2: Generate Surface Slab

### Basic Surface Creation

```python
from nh3sofc.structure import SurfaceBuilder

builder = SurfaceBuilder(bulk)

# Create (001) surface
surface = builder.create_surface(
    miller_index=(0, 0, 1),
    layers=6,
    vacuum=15.0,
    fix_bottom=2
)

print(f"Surface atoms: {len(surface.atoms)}")
print(f"Cell: {surface.atoms.cell.lengths()}")
```

### Different Terminations

```python
from nh3sofc import write_poscar
import os

# Create work directory for this surface study
work_dir = "work/surfaces/LaVO3_001"
if not os.path.exists(work_dir):
    os.makedirs(work_dir)

# Get all possible terminations
terminations = builder.get_all_terminations((0, 0, 1))

for i, term in enumerate(terminations):
    info = builder.identify_termination(term)
    print(f"Termination {i}: {info['composition']}")
    # Save as POSCAR with proper format (atoms sorted by element automatically)
    write_poscar(term.atoms, f"{work_dir}/POSCAR_term_{i}")
```

### Creating Supercells

```python
# Method 1: Using supercell parameter (recommended)
surface = builder.create_surface_with_size(
    miller_index=(0, 0, 1),
    size=(2, 2, 6),  # 2x2 supercell, 6 layers
    vacuum=15.0
)
print(f"Surface atoms: {len(surface.atoms)}")
print(f"Cell: {surface.atoms.cell.lengths()}")

# Method 2: Repeat after creation (preserves constraints)
surface = builder.create_surface(
    miller_index=(0, 0, 1),
    layers=6,
    fix_bottom=2
)
surface = surface.repeat_xy(2, 2)  # Fixed atoms are preserved
print(f"Surface atoms: {len(surface.atoms)}")
print(f"Cell: {surface.atoms.cell.lengths()}")
```

## Step 3: Handling Surface Polarity

Many oxide surfaces are polar, meaning they have alternating charged layers that create an electric dipole. This is critical for accurate DFT calculations.

### Check Surface Polarity

```python
# Analyze polarity
polarity = surface.check_polarity()

print(f"Is polar: {polarity['is_polar']}")
print(f"Dipole (z-component): {polarity['dipole_z']:.2f} e·Å")
print(f"Layer charges: {polarity['layer_charges']}")
print(f"Recommendation: {polarity['recommendation']}")
```

### Analyze Layers

```python
# Identify atomic layers (uses automatic tolerance based on covalent radii)
layers = surface.identify_layers()  # tolerance="auto" by default

for i, layer in enumerate(layers):
    print(f"Layer {i}: z={layer['z']:.2f} Å, composition={layer['composition']}")

# Check what tolerance was calculated
tol = surface.estimate_layer_tolerance()
print(f"Auto-calculated tolerance: {tol:.2f} Å")

# For complex materials like perovskites, auto tolerance correctly groups
# atoms that belong to the same layer (e.g., V and O in VO2 layers)
# You can override with explicit tolerance if needed:
layers_custom = surface.identify_layers(tolerance=0.5)

# Get layer spacing
spacings = surface.get_layer_spacing()
print(f"Interlayer distances: {spacings}")
```

### Create Symmetric Slabs

For polar surfaces, symmetric slabs (same termination on both sides) help cancel the dipole. NH3SOFC creates truly symmetric slabs by trimming to matching top/bottom layer compositions:

```python
# Create symmetric slab (automatically trims to matching terminations)
supercell = bulk.make_supercell((2, 2, 1))
builder = SurfaceBuilder(supercell)
symmetric_surface = builder.create_symmetric_slab(
    #termination={"V":1, "O":2},
    miller_index=(0, 0, 1),
    layers=4,
    vacuum=15.0,
    fix_bottom=2
)

# Verify top and bottom layers match
layers = symmetric_surface.identify_layers()
print(f"Top layer: {layers[-1]['composition']}")
print(f"Bottom layer: {layers[0]['composition']}")

# Verify reduced dipole
polarity = symmetric_surface.check_polarity()
print(f"Dipole after symmetrization: {polarity['dipole_z']:.2f} e·Å")
print(f"Total number of atoms: {symmetric_surface.atoms.get_number_of_atoms()}")
```

#### Requesting Specific Terminations

You can request a specific termination for your symmetric slab. The termination is specified as an **element ratio**, so `{"La": 1, "O": 1}` and `{"La": 8, "O": 8}` give identical results:

```python
# Create LaO-terminated symmetric slab
lao_symmetric = builder.create_symmetric_slab(
    miller_index=(0, 0, 1),
    layers=7,
    vacuum=15.0,
    termination={"La": 1, "O": 1},  # LaO composition ratio (1:1)
    min_layers=5,
    fix_bottom=2
)

# Create VO2-terminated symmetric slab
vo2_symmetric = builder.create_symmetric_slab(
    miller_index=(0, 0, 1),
    layers=7,
    vacuum=15.0,
    termination={"V": 1, "O": 2},  # VO2 composition ratio (1:2)
    min_layers=5,
    fix_bottom=2
)
```

When no termination is specified, the builder automatically tries all unique layer compositions and picks the one that creates the best symmetric slab (most layers). This usually works well, but for precise control, specifying the termination explicitly is recommended.

#### Manual Trimming for Fine Control

For more control, you can create an oversized slab and trim it yourself:

```python
# Create oversized slab
oversized = builder.create_surface(
    miller_index=(0, 0, 1),
    layers=6,
    vacuum=15.0
)

# Trim to symmetric LaO termination
symmetric = oversized.trim_to_symmetric_termination(
    termination={"La": 1, "O": 1},
    min_layers=5
)

# Verify the result
layers = symmetric.identify_layers()
print(f"Total number of atoms: {symmetric.atoms.get_number_of_atoms()}")
print(f"Total layers: {len(layers)}")
print(f"Top: {layers[-1]['composition']}")
print(f"Bottom: {layers[0]['composition']}")
```

## Step 4: Stoichiometry Validation

Check if your surface has the expected composition:

```python
# Get stoichiometry
stoich = surface.get_stoichiometry()
print(f"Normalized stoichiometry: {stoich}")

# Check against expected
result = surface.check_stoichiometry(
    expected={"La": 1, "V": 1, "O": 3},
    tolerance=0.1
)
print(f"Is stoichiometric: {result['is_stoichiometric']}")
if result['warnings']:
    print(f"Warnings: {result['warnings']}")

# Layer-by-layer stoichiometry
layers = surface.get_layer_stoichiometry()
for layer in layers:
    print(f"z={layer['z']:.1f}: {layer['composition']}")
```

## Step 5: Create Oxynitride Structure

```python
from nh3sofc.structure import DefectBuilder

defect = DefectBuilder(symmetric_surface)

# Create oxynitride: replace 2/3 of O with N, add 10% vacancies
oxynitride = defect.create_oxynitride(
    nitrogen_fraction=0.67,      # 2/3 O → N
    vacancy_concentration=0.10,  # 10% vacancies (x in LaVON2-x)
)

print(f"Formula: {oxynitride.get_chemical_formula()}")
```

### Defect Placement Strategies

Control where defects are placed using `placement` and preference parameters.
The placement uses **probability-weighted selection** - atoms in preferred regions
have higher probability of being selected, but selection is still stochastic,
so each configuration is unique.

```python
# Random placement (default) - uniform probability
oxynitride_random = defect.create_oxynitride(
    vacancy_element="N",
    nitrogen_fraction=0.67,
    vacancy_concentration=0.10,
    placement="random"
)

# Surface-populated: N preferentially at surface
# surface_n_preference controls how strongly N prefers surface:
#   0.5 = random distribution
#   0.7 = ~70% of surface anions will be N (default)
#   0.9 = ~90% of surface anions will be N
oxynitride_surface = defect.create_oxynitride(
    vacancy_element="N",
    nitrogen_fraction=0.67,
    vacancy_concentration=0.10,
    placement="surface",
    surface_n_preference=0.8,   # 80% of surface anions are N
    vacancy_preference=0.6,     # Vacancies slightly prefer surface
)

# Near-N: vacancies preferentially near existing N atoms
# vacancy_preference controls how strongly vacancies prefer being near N
oxynitride_near_n = defect.create_oxynitride(
    vacancy_element="N",
    nitrogen_fraction=0.67,
    vacancy_concentration=0.10,
    placement="near_N",
    vacancy_preference=0.7,     # Vacancies moderately prefer being near N
)
```

#### Preference Parameters

| Parameter | Range | Description |
|-----------|-------|-------------|
| `surface_n_preference` | 0.5-1.0 | For "surface" placement: fraction of surface anions that are N |
| `vacancy_preference` | 0.5-1.0 | Preference strength for vacancy placement (surface or near-N) |

- **0.5** = Random (no preference)
- **0.7** = Moderate preference (default)
- **0.9-1.0** = Strong preference (most defects in target region)

### Generate Configuration Pool

For screening studies, generate multiple configurations with all strategies.
Each configuration is unique due to probability-weighted selection:

```python
from nh3sofc import save_configurations

# Generate pool of candidate structures with preference settings
pool = defect.create_oxynitride_pool(
    nitrogen_fraction=0.67,
    vacancy_concentration=0.10,
    n_configs_per_strategy=5,       # 5 configs × 3 strategies = 15 total
    surface_n_preference=0.8,       # For "surface": 80% surface anions are N
    vacancy_preference=0.6,         # Moderate vacancy preference
)

print(f"Generated {len(pool)} configurations")

# Extract atoms and create descriptive names
configs = [config['atoms'] for config in pool]
names = [f"{config['placement']}_{config['config_id']}" for config in pool]

# Option 1: Auto-generate work directory (creates folder like "work/slab_La8V8N16O8_40atoms/")
result = save_configurations(configs, names=names, prefix="oxynitride")
print(f"Saved to: {result['work_dir']}")

# Option 2: Specify work directory explicitly
result = save_configurations(configs, "work/oxynitrides/LaVON2_001", names=names)

for config, p in zip(pool, result["configs"]):
    print(f"Saved {p['poscar'].name}: {config['atoms'].get_chemical_formula()}")

# Each config dict contains:
# - atoms: the Atoms object
# - placement: "random", "surface", or "near_N"
# - surface_n_preference, vacancy_preference: the preference values used
# - nitrogen_fraction, vacancy_concentration
# - config_id: unique identifier
```

### Analyze Defect Distribution

Use statistical analysis to verify that defects are distributed as expected:

```python
from nh3sofc.structure import (
    analyze_defect_distribution,
    analyze_oxynitride_pool,
    print_defect_analysis,
)

# Analyze a single configuration
stats = analyze_defect_distribution(
    oxynitride_surface,
    z_threshold=0.3,  # Top 30% of slab height is "surface"
)

print(f"N atoms in surface region: {stats['n_surface']} / {stats['n_total']}")
print(f"Surface N/(N+O) ratio: {stats['surface_n_ratio']:.1%}")
print(f"Bulk N/(N+O) ratio: {stats['bulk_n_ratio']:.1%}")

# With vacancy analysis (requires original structure for comparison)
stats = analyze_defect_distribution(
    oxynitride_surface,
    reference_atoms=symmetric_surface.atoms,  # Original before defects
    z_threshold=0.3,
    near_n_cutoff=3.0,  # Count vacancies within 3 Å of N
)

print(f"Vacancies in surface: {stats['vacancy_surface']} / {stats['vacancy_total']}")
print(f"Vacancies near N: {stats['vacancy_near_n']} ({stats['vacancy_near_n_fraction']:.1%})")

# Analyze entire pool and compare strategies
pool_stats = analyze_oxynitride_pool(pool, z_threshold=0.3)

# Print formatted summary
print_defect_analysis(pool_stats, title="Oxynitride Pool Analysis")
```

Example output:
```
============================================================
Oxynitride Pool Analysis
============================================================

Total configurations: 15

RANDOM Strategy (5 configs):
  Surface N ratio: 48.2% ± 5.3%
  Bulk N ratio: 51.1% ± 4.8%
  N in surface: 32.5% ± 3.2%

SURFACE Strategy (5 configs):
  Surface N ratio: 72.4% ± 4.1%
  Bulk N ratio: 38.6% ± 3.9%
  N in surface: 48.3% ± 2.8%

NEAR_N Strategy (5 configs):
  Surface N ratio: 50.1% ± 6.2%
  Bulk N ratio: 49.8% ± 5.1%
  Vacancies near N: 68.3% ± 8.5%
============================================================
```

## Step 6: Visualize and Save

```python
from nh3sofc import save_structure, write_poscar
from ase.visualize import view

# Option 1: Auto-generate meaningful work directory (recommended)
# Creates folder like "work/slab_La8V8O24_40atoms/"
paths = save_structure(surface.atoms, name="surface", formats=["poscar", "cif"])
print(f"Saved to: {paths['work_dir']}")

# Option 2: Specify work directory explicitly
work_dir = "work/surfaces/LaVO3_001"
save_structure(bulk.atoms, work_dir, "bulk", formats=["poscar", "cif"])
save_structure(surface.atoms, work_dir, "surface", formats=["poscar", "cif"])
save_structure(oxynitride, work_dir, "oxynitride", formats=["poscar", "cif"])

# Option 3: Save individual POSCAR files directly
write_poscar(bulk.atoms, f"{work_dir}/POSCAR_bulk")
write_poscar(surface.atoms, f"{work_dir}/POSCAR_surface")
write_poscar(oxynitride, f"{work_dir}/POSCAR_oxynitride")

# Visualize (if GUI available)
# view(oxynitride)
```

## Best Practices

### Surface Construction

1. **Check polarity** - Always check if your surface is polar using `check_polarity()`
2. **Use symmetric slabs** - For polar surfaces, use `symmetric=True` or `create_symmetric_slab()`
3. **Sufficient layers** - Use at least 5-7 layers for accurate surface properties
4. **Vacuum spacing** - 15-20 Å is typically sufficient to avoid image interactions
5. **Fix bottom atoms** - Fix 2-3 bottom layers to simulate bulk behavior
6. **Layer identification** - Use `identify_layers()` with default auto-tolerance; it correctly groups atoms in complex materials like perovskites where atoms in the same layer may have slightly different z-positions

### Polarity Handling

| Surface Type | Polarity | Recommendation |
|-------------|----------|----------------|
| Perovskite (001) | Polar | Use symmetric slab or dipole correction |
| Perovskite (110) | Polar | Use symmetric slab |

### VASP Settings for Polar Surfaces

If you cannot use a symmetric slab, add dipole corrections to VASP:

```python
from nh3sofc.calculators.vasp import VASPInputGenerator

vasp = VASPInputGenerator(surface.atoms, calc_type="relax")
vasp.set_surface_settings()  # Adds IDIPOL=3, LDIPOL=.TRUE.
```

## Complete Example

```python
from nh3sofc.structure import (
    BulkStructure,
    SurfaceBuilder,
    DefectBuilder
)
from nh3sofc import save_structure, save_configurations

# 1. Load bulk perovskite
bulk = BulkStructure.from_cif("LaVO3.cif")

# 2. Create supercell for sufficient surface area
supercell = bulk.make_supercell((2, 2, 1))

# 3. Create symmetric surface with LaO termination
builder = SurfaceBuilder(supercell)

surface = builder.create_symmetric_slab(
    miller_index=(0, 0, 1),
    termination={"La": 1, "O": 1},  # LaO termination
    layers=7,
    vacuum=15.0,
    fix_bottom=2,
)

# 4. Check polarity
polarity = surface.check_polarity()
print(f"Dipole: {polarity['dipole_z']:.2f} e·Å")

# Save the clean surface (auto-generates work_dir like "work/slab_La8V8O24_40atoms/")
paths = save_structure(surface.atoms, name="clean_surface", prefix="LaO_term")
work_dir = paths["work_dir"]
print(f"Clean surface saved to: {work_dir}")

# 5. Generate oxynitride configuration pool
defect = DefectBuilder(surface)
pool = defect.create_oxynitride_pool(
    nitrogen_fraction=0.67,
    vacancy_concentration=0.10,
    n_configs_per_strategy=5,
)

print(f"Generated {len(pool)} configurations")

# 6. Save all configurations (in same work_dir as clean surface)
configs = [config['atoms'] for config in pool]
names = [f"{config['placement']}_{config['config_id']}" for config in pool]

result = save_configurations(configs, work_dir, names=names)

for config, p in zip(pool, result["configs"]):
    print(f"Saved {p['poscar'].name}: {config['atoms'].get_chemical_formula()}")

print(f"\nAll configurations saved to: {work_dir}")
```

## Next Steps

- [Adsorbate Placement](adsorbate_placement.md) - Place NH3 on the surface
- [VASP Calculations](vasp_calculations.md) - Optimize the structure
