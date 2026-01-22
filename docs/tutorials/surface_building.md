# Tutorial: Surface Building

This tutorial covers building surfaces from bulk structures, handling surface polarity, and creating oxynitride defects.

## Learning Objectives

- Load bulk structures from CIF files
- Generate surface slabs with different Miller indices
- Handle polar surfaces (critical for perovskites)
- Use specialized builders for different crystal structures
- Create oxynitride structures with controlled defect concentrations

## Step 1: Load Bulk Structure

```python
from nh3sofc.structure import BulkStructure

# Load from CIF
bulk = BulkStructure.from_cif("LaVO3.cif")
print(f"Formula: {bulk.atoms.get_chemical_formula()}")
print(f"Space group: {bulk.get_spacegroup()}")

# Create supercell if needed
supercell = bulk.make_supercell((2, 2, 2))
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
# Get all possible terminations
terminations = builder.get_all_terminations((0, 0, 1))

for i, term in enumerate(terminations):
    info = builder.identify_termination(term)
    print(f"Termination {i}: {info['composition']}")
    term.write(f"surface_term_{i}.xyz")
```

### Creating Supercells

```python
# Method 1: Using supercell parameter (recommended)
surface = builder.create_surface_with_size(
    miller_index=(0, 0, 1),
    size=(2, 2, 6),  # 2x2 supercell, 6 layers
    vacuum=15.0
)

# Method 2: Repeat after creation (preserves constraints)
surface = builder.create_surface(
    miller_index=(0, 0, 1),
    layers=6,
    fix_bottom=2
)
surface = surface.repeat_xy(2, 2)  # Fixed atoms are preserved
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
symmetric_surface = builder.create_symmetric_slab(
    miller_index=(0, 0, 1),
    layers=7,
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
    layers=12,
    vacuum=15.0
)

# Trim to symmetric LaO termination
symmetric = oversized.trim_to_symmetric_termination(
    termination={"La": 1, "O": 1},
    min_layers=5
)

# Verify the result
layers = symmetric.identify_layers()
print(f"Total layers: {len(layers)}")
print(f"Top: {layers[-1]['composition']}")
print(f"Bottom: {layers[0]['composition']}")
```

## Step 4: Specialized Surface Builders

NH3SOFC provides specialized builders for common crystal structure types:

### Perovskite Surfaces (ABO3)

For LaVO3, SrTiO3, BaTiO3, etc.:

```python
from nh3sofc.structure import PerovskiteSurfaceBuilder

# Create builder with automatic site detection
builder = PerovskiteSurfaceBuilder(bulk)  # Auto-detects La as A-site, V as B-site

# Or specify explicitly
builder = PerovskiteSurfaceBuilder(bulk, A_site="La", B_site="V", anion="O")

# See available terminations
print(builder.get_termination_options((0, 0, 1)))  # ['LaO', 'VO2']

# Create surface with specific termination
surface = builder.create_surface(
    miller_index=(0, 0, 1),
    termination="LaO",      # Named termination (also accepts "AO" or "BO2")
    layers=7,
    symmetric=True,         # Symmetric slab (recommended for polar surfaces)
    fix_bottom=2,
    supercell=(2, 2)        # 2x2 supercell
)

# Analyze the surface
analysis = builder.analyze_surface(surface)
print(f"Termination: {analysis['termination_type']}")
print(f"Surface composition: {analysis['surface_composition']}")
```

When `symmetric=True`, the builder creates a truly symmetric slab where both top and bottom layers have the same composition. For example, with `termination="LaO"`:

```python
# Verify both surfaces are LaO-terminated
layers = surface.identify_layers()
print(f"Top layer: {layers[-1]['composition']}")     # {'La': N, 'O': N}
print(f"Bottom layer: {layers[0]['composition']}")   # {'La': N, 'O': N}

# VO2 termination works the same way
vo2_surface = builder.create_surface(
    miller_index=(0, 0, 1),
    termination="VO2",
    layers=7,
    symmetric=True
)
layers = vo2_surface.identify_layers()
print(f"Top: {layers[-1]['composition']}")    # {'V': N, 'O': 2N}
print(f"Bottom: {layers[0]['composition']}")  # {'V': N, 'O': 2N}
```

### Rocksalt Surfaces (MX)

For NiO, MgO, CoO, etc.:

```python
from nh3sofc.structure import RocksaltSurfaceBuilder
from ase.build import bulk

# Create NiO bulk
nio = bulk('NiO', 'rocksalt', a=4.17)

builder = RocksaltSurfaceBuilder(nio, cation="Ni", anion="O")

# (001) and (110) are non-polar
surface_001 = builder.create_surface(
    miller_index=(0, 0, 1),
    layers=6
)

# (111) is polar - use symmetric slab
surface_111 = builder.create_surface(
    miller_index=(1, 1, 1),
    layers=7,
    symmetric=True  # Important for polar (111)
)
```

### Fluorite Surfaces (MX2)

For CeO2, ZrO2, etc.:

```python
from nh3sofc.structure import FluoriteSurfaceBuilder
from ase.build import bulk

# Create CeO2 bulk
ceo2 = bulk('CeO2', 'fluorite', a=5.41)

builder = FluoriteSurfaceBuilder(ceo2, cation="Ce", anion="O")

# (111) is most stable and non-polar
surface_111 = builder.create_surface(
    miller_index=(1, 1, 1),
    layers=6
)

# (110) and (100) are polar
surface_110 = builder.create_surface(
    miller_index=(1, 1, 0),
    symmetric=True  # Needed for polar surfaces
)
```

## Step 5: Stoichiometry Validation

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

## Step 6: Create Oxynitride Structure

```python
from nh3sofc.structure import DefectBuilder

defect = DefectBuilder(surface)

# Create oxynitride: replace 2/3 of O with N, add 10% vacancies
oxynitride = defect.create_oxynitride(
    nitrogen_fraction=0.67,      # 2/3 O → N
    vacancy_concentration=0.10,  # 10% vacancies (x in LaVON2-x)
    prefer_surface=True          # Preferentially modify surface atoms
)

print(f"Formula: {oxynitride.atoms.get_chemical_formula()}")
```

### Create Specific Vacancies

```python
# Create oxygen vacancy at specific site
vacancy_structure = defect.create_vacancy(
    element="O",
    n_vacancies=2
)
```

## Step 7: Visualize and Save

```python
from ase.io import write
from ase.visualize import view

# Save structures
write("bulk.xyz", bulk.atoms)
write("surface.xyz", surface.atoms)
write("oxynitride.xyz", oxynitride.atoms)

# Visualize (if GUI available)
# view(oxynitride.atoms)
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
| Rocksalt (001) | Non-polar | Standard slab OK |
| Rocksalt (111) | Polar | Use symmetric slab |
| Fluorite (111) | Non-polar | Standard slab OK |
| Fluorite (110) | Polar | Use symmetric slab |

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
    PerovskiteSurfaceBuilder,
    DefectBuilder
)
from ase.io import write

# 1. Load bulk perovskite
bulk = BulkStructure.from_cif("LaVO3.cif")

# 2. Create surface with specialized builder
builder = PerovskiteSurfaceBuilder(bulk, A_site="La", B_site="V")

surface = builder.create_surface(
    miller_index=(0, 0, 1),
    termination="LaO",
    layers=7,
    symmetric=True,
    fix_bottom=2,
    supercell=(2, 2)
)

# 3. Check polarity
polarity = surface.check_polarity()
print(f"Dipole: {polarity['dipole_z']:.2f} e·Å")

# 4. Create oxynitride
defect = DefectBuilder(surface)
oxynitride = defect.create_oxynitride(
    nitrogen_fraction=0.67,
    vacancy_concentration=0.10,
)

# 5. Validate stoichiometry
stoich = oxynitride.check_stoichiometry()
print(f"Final composition: {oxynitride.formula}")

# 6. Save
write("LaVON2_001_symmetric.xyz", oxynitride.atoms)
```

## Next Steps

- [Adsorbate Placement](adsorbate_placement.md) - Place NH3 on the surface
- [VASP Calculations](vasp_calculations.md) - Optimize the structure
