# Tutorial: Surface Stability Analysis

This tutorial covers calculating and comparing surface energies to determine
which crystal facets are thermodynamically most stable.

## Learning Objectives

- Understand surface energy and its calculation
- Build surfaces for different Miller indices
- Calculate surface energies from DFT results
- Rank surfaces by thermodynamic stability
- Handle polar and non-stoichiometric surfaces

## Theory: Surface Energy

Surface energy (γ) measures the energy cost of creating a surface:

```
γ = (E_slab - n × E_bulk) / (2 × A)
```

Where:
- `E_slab` = Total energy of the slab (from DFT)
- `E_bulk` = Bulk energy per atom (from DFT)
- `n` = Number of atoms in the slab
- `A` = Surface area
- Factor of 2 accounts for two surfaces (top and bottom)

**Units:** Typically reported in J/m² (1 eV/Å² = 16.02 J/m²)

**Lower surface energy = more stable surface**

## Step 1: Build the Bulk Structure

```python
from ase.spacegroup import crystal
from nh3sofc.structure import SurfaceBuilder

# Example: CeO2 fluorite structure (Fm-3m, space group 225)
ceo2_bulk = crystal(
    symbols=['Ce', 'O'],
    basis=[(0, 0, 0), (0.25, 0.25, 0.25)],
    spacegroup=225,
    cellpar=[5.411, 5.411, 5.411, 90, 90, 90]
)

print(f"Bulk: {ceo2_bulk.get_chemical_formula()}")
print(f"Atoms per unit cell: {len(ceo2_bulk)}")
```

## Step 2: Build Surfaces for Different Miller Indices

```python
import numpy as np
from nh3sofc.structure import SurfaceBuilder

builder = SurfaceBuilder(ceo2_bulk)

# Common low-index surfaces
miller_indices = [(1,1,1), (1,1,0), (1,0,0)]

surfaces = {}
for miller in miller_indices:
    miller_str = ''.join(map(str, miller))

    slab = builder.create_surface(
        miller_index=miller,
        layers=6,          # Enough layers for convergence
        vacuum=15.0,       # Vacuum thickness in Å
    )

    # Calculate surface area
    cell = slab.atoms.cell
    area = np.linalg.norm(np.cross(cell[0], cell[1]))

    surfaces[miller_str] = {
        'slab': slab,
        'n_atoms': len(slab.atoms),
        'area': area,
    }

    # Check stoichiometry
    n_ce = sum(1 for a in slab.atoms if a.symbol == 'Ce')
    n_o = sum(1 for a in slab.atoms if a.symbol == 'O')

    print(f"({miller_str}): {len(slab.atoms)} atoms, "
          f"area = {area:.2f} Å², Ce{n_ce}O{n_o}")
```

## Step 3: Run DFT Calculations

Save structures and run VASP calculations:

```python
from nh3sofc.calculators import VASPInputGenerator
from nh3sofc import write_poscar

# Save bulk for reference calculation
write_poscar(ceo2_bulk, "work/CeO2_bulk/POSCAR")

# Generate VASP inputs for each surface
for miller_str, data in surfaces.items():
    work_dir = f"work/CeO2_{miller_str}"

    vasp = VASPInputGenerator(
        data['slab'].atoms,
        calc_type='relax',
        work_dir=work_dir,
        encut=520,
        kspacing=0.03,
    )
    vasp.generate_all()

    print(f"Generated inputs: {work_dir}")
```

## Step 4: Calculate Surface Energies

After DFT calculations complete:

```python
from nh3sofc.analysis import SurfaceEnergyCalculator

# Get bulk energy per atom from your DFT calculation
# E_bulk = total energy / n_atoms from bulk OUTCAR
e_bulk_per_atom = -9.876  # eV/atom (replace with your value)

# Initialize calculator
calc = SurfaceEnergyCalculator(e_bulk_per_atom=e_bulk_per_atom)

# Calculate surface energy for each facet
results = {}

for miller_str, data in surfaces.items():
    # Get slab energy from VASP OUTCAR
    # e_slab = parse_outcar(f"work/CeO2_{miller_str}/OUTCAR")
    e_slab = -712.345  # Example - replace with actual value

    # Calculate in eV/Å²
    gamma_ev = calc.calculate(
        e_slab=e_slab,
        n_atoms=data['n_atoms'],
        area=data['area'],
        symmetric=True,  # Top and bottom surfaces equivalent
    )

    # Convert to J/m²
    gamma_j = calc.calculate_j_m2(
        e_slab=e_slab,
        n_atoms=data['n_atoms'],
        area=data['area'],
        symmetric=True,
    )

    results[miller_str] = {
        'gamma_ev_A2': gamma_ev,
        'gamma_J_m2': gamma_j,
        'area': data['area'],
    }
```

## Step 5: Rank Surfaces by Stability

```python
# Sort by surface energy (ascending = most stable first)
ranked = sorted(results.items(), key=lambda x: x[1]['gamma_J_m2'])

print("\n" + "="*60)
print("Surface Stability Ranking")
print("="*60)
print(f"\nBulk reference: {e_bulk_per_atom:.4f} eV/atom")
print("\n{:<10} {:>12} {:>12} {:>15}".format(
    "Surface", "γ (eV/Å²)", "γ (J/m²)", "Stability"
))
print("-"*60)

for i, (miller, data) in enumerate(ranked):
    stability = "Most stable" if i == 0 else ""
    print(f"({miller})      {data['gamma_ev_A2']:>12.4f} "
          f"{data['gamma_J_m2']:>12.2f} {stability:>15}")

print("\n" + "="*60)
stability_order = " > ".join([f"({m})" for m, _ in ranked])
print(f"Stability order: {stability_order}")
print("="*60)
```

## Step 6: Visualize Results

Use the built-in plotting function:

```python
from nh3sofc.analysis import plot_surface_stability

# Collect surface energies
energies = {miller: data['gamma_J_m2'] for miller, data in results.items()}

# Plot (sorted by stability automatically)
fig, ax = plot_surface_stability(
    energies,
    title="CeO2 Surface Stability",
    filename="surface_stability.png"
)

# Output: Bar chart with green (stable) to red (unstable) colors
```

The plot shows:
- Bars sorted by stability (most stable on left)
- Color gradient: green (stable) → red (unstable)
- Energy values labeled on each bar

## Handling Special Cases

### Polar Surfaces

Some surfaces (e.g., CeO2 (100)) are polar and may require:
- Reconstruction
- Non-stoichiometric termination
- Chemical potential corrections

```python
# For non-stoichiometric slabs, use chemical potentials
def surface_energy_nonstoich(e_slab, n_ce, n_o, area,
                              e_bulk_ceo2, mu_o):
    """
    γ = (E_slab - n_Ce×μ_Ce - n_O×μ_O) / (2A)

    With constraint: μ_Ce + 2μ_O = E_bulk(CeO2)
    So: μ_Ce = E_bulk - 2μ_O
    """
    mu_ce = e_bulk_ceo2 - 2 * mu_o
    e_ref = n_ce * mu_ce + n_o * mu_o
    gamma = (e_slab - e_ref) / (2 * area)
    return gamma
```

### Asymmetric Slabs

If top and bottom surfaces are different:

```python
gamma = calc.calculate(
    e_slab=e_slab,
    n_atoms=n_atoms,
    area=area,
    symmetric=False,  # Only one surface counted
)
```

### Convergence Testing

Surface energy should converge with slab thickness:

```python
# Test convergence with number of layers
layers_to_test = [4, 6, 8, 10]
convergence = {}

for n_layers in layers_to_test:
    slab = builder.create_surface(
        miller_index=(1,1,1),
        layers=n_layers,
        vacuum=15.0,
    )
    # Run DFT and calculate γ
    # convergence[n_layers] = gamma

# Plot convergence
plt.plot(layers_to_test, [convergence[n] for n in layers_to_test], 'o-')
plt.xlabel('Number of layers')
plt.ylabel('Surface energy (J/m²)')
plt.title('Convergence test')
```

## Complete Example: CeO2

```python
from ase.spacegroup import crystal
from ase.io import read
import numpy as np
from nh3sofc.structure import SurfaceBuilder
from nh3sofc.analysis import SurfaceEnergyCalculator
from nh3sofc.calculators import VASPInputGenerator

# 1. Create bulk CeO2
ceo2_bulk = crystal(
    symbols=['Ce', 'O'],
    basis=[(0, 0, 0), (0.25, 0.25, 0.25)],
    spacegroup=225,
    cellpar=[5.411, 5.411, 5.411, 90, 90, 90]
)

# 2. Build surfaces
builder = SurfaceBuilder(ceo2_bulk)
miller_indices = [(1,1,1), (1,1,0), (1,0,0)]

for miller in miller_indices:
    miller_str = ''.join(map(str, miller))
    slab = builder.create_surface(miller, layers=6, vacuum=15.0)

    # Save for DFT
    vasp = VASPInputGenerator(
        slab.atoms,
        calc_type='relax',
        work_dir=f"work/CeO2_{miller_str}",
    )
    vasp.generate_all(encut=520)

# 3. After DFT: Parse energies and calculate γ
# e_bulk_per_atom = parse_outcar("work/CeO2_bulk/OUTCAR") / len(ceo2_bulk)
e_bulk_per_atom = -9.876  # Replace with actual

calc = SurfaceEnergyCalculator(e_bulk_per_atom)

results = {}
for miller in miller_indices:
    miller_str = ''.join(map(str, miller))
    slab = read(f"work/CeO2_{miller_str}/CONTCAR")
    # e_slab = parse_outcar(f"work/CeO2_{miller_str}/OUTCAR")
    e_slab = -700.0  # Replace with actual

    area = np.linalg.norm(np.cross(slab.cell[0], slab.cell[1]))
    gamma = calc.calculate_j_m2(e_slab, len(slab), area)
    results[miller_str] = gamma

# 4. Rank
ranked = sorted(results.items(), key=lambda x: x[1])
print("Stability:", " > ".join([f"({m})" for m, _ in ranked]))
```

## Reference: Typical Surface Energies

| Material | Surface | γ (J/m²) | Notes |
|----------|---------|----------|-------|
| CeO2 | (111) | 0.7-1.0 | Most stable, O-terminated |
| CeO2 | (110) | 1.0-1.5 | Intermediate |
| CeO2 | (100) | 2.0-2.5 | Polar, may reconstruct |
| TiO2 rutile | (110) | 0.5-0.8 | Most stable |
| TiO2 rutile | (100) | 0.9-1.2 | |
| Pt | (111) | 1.5-2.0 | FCC most stable |
| Pt | (100) | 2.0-2.5 | |

## Next Steps

- [Surface Comparison](surface_comparison.md) - Compare catalytic activity
- [Adsorbate Placement](adsorbate_placement.md) - Add molecules to stable surfaces
- [Defect Building](defects.md) - Create oxygen vacancies
