# Tutorial: Surface Building

This tutorial covers building surfaces from bulk structures and creating oxynitride defects.

## Learning Objectives

- Load bulk structures from CIF files
- Generate surface slabs with different Miller indices
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
    print(f"Termination {i}: {term.get_chemical_formula()}")
    term.write(f"surface_term_{i}.xyz")
```

## Step 3: Create Oxynitride Structure

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

## Step 4: Visualize and Save

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

1. **Check symmetry** - Verify surface symmetry after creation
2. **Sufficient layers** - Use at least 4-6 layers for accurate surface properties
3. **Vacuum spacing** - 15 Å is typically sufficient to avoid image interactions
4. **Fix bottom atoms** - Fix 2-3 bottom layers to simulate bulk behavior

## Complete Example

```python
from nh3sofc.structure import BulkStructure, SurfaceBuilder, DefectBuilder
from ase.io import write

# 1. Load bulk
bulk = BulkStructure.from_cif("LaVO3.cif")

# 2. Create surface
builder = SurfaceBuilder(bulk)
surface = builder.create_surface(
    miller_index=(0, 0, 1),
    layers=6,
    vacuum=15.0,
)

# 3. Create oxynitride
defect = DefectBuilder(surface)
oxynitride = defect.create_oxynitride(
    nitrogen_fraction=0.67,
    vacancy_concentration=0.10,
)

# 4. Save
write("LaVON2_001.xyz", oxynitride.atoms)
```

## Next Steps

- [Adsorbate Placement](adsorbate_placement.md) - Place NH3 on the surface
- [VASP Calculations](vasp_calculations.md) - Optimize the structure
