# Doped Ceria Tutorial (GDC, SDC, PDC)

This tutorial covers creating acceptor-doped ceria structures for solid oxide fuel cell applications.

## Overview

Doped ceria is crucial for ionic conductivity in SOFCs:

- **GDC**: Gadolinium-doped ceria (Ce₁₋ₓGdₓO₂₋δ)
- **SDC**: Samarium-doped ceria (Ce₁₋ₓSmₓO₂₋δ)
- **PDC**: Praseodymium-doped ceria (Ce₁₋ₓPrₓO₂₋δ)

### Chemistry

Trivalent dopants (M³⁺) substitute Ce⁴⁺ sites as acceptor dopants. Charge compensation requires oxygen vacancy formation:

```
2 Ce⁴⁺ → 2 M³⁺ + V_O^{2+}
```

This gives a 2:1 dopant-to-vacancy ratio for charge neutrality.

!!! note "Odd Number of Dopants"
    When the number of dopants is odd (e.g., 9 Gd³⁺), the ideal vacancy count would be
    non-integer (4.5). The code uses floor division, creating 4 vacancies in this case.

    The uncompensated negative charge physically corresponds to Ce⁴⁺ → Ce³⁺ reduction
    (small polaron formation), which is realistic for reducible CeO₂.

    **Options:**

    - Adjust `dopant_fraction` to get an even number of dopants for exact 2:1 stoichiometry
    - Accept the implicit Ce³⁺ formation (physically reasonable)

## Basic Usage

### Creating Gd-Doped Ceria (GDC)

```python
from ase.build import bulk
from nh3sofc.structure import DopantBuilder

# Create CeO2 supercell
ceo2 = bulk('CeO2', 'fluorite', a=5.41)
ceo2_slab = ceo2 * (2, 2, 2)  # 8 Ce, 16 O

# Initialize builder
builder = DopantBuilder(ceo2_slab)

# Create 10% Gd-doped ceria
gdc = builder.create_doped_structure(
    dopant="Gd",
    dopant_fraction=0.10,  # 10% of Ce → Gd
    random_seed=42,
)

# Check composition
symbols = gdc.get_chemical_symbols()
n_gd = symbols.count("Gd")
n_ce = symbols.count("Ce")
n_o = symbols.count("O")

print(f"Gd atoms: {n_gd}")        # ~1 (10% of 8)
print(f"Ce atoms: {n_ce}")        # ~7
print(f"O atoms: {n_o}")          # 16 - vacancies
print(f"Charge balanced: 2 Gd → 1 vacancy")
```

### Using DopedCeriaStructure

The `DopedCeriaStructure` class tracks dopant metadata:

```python
from nh3sofc.structure import DopedCeriaStructure

# Create from pristine ceria
doped = DopedCeriaStructure.from_ceria(
    ceo2_slab,
    dopant="Sm",
    dopant_fraction=0.15,
    random_seed=42,
)

print(doped)  # DopedCeriaStructure(SDC, Sm_frac=15.0%, n_vac=1)
print(doped.get_dopant_name())  # "SDC"
print(doped.get_stoichiometry())  # {"Ce": 0.875, "Sm": 0.125, "O": 1.9375}
```

## Placement Strategies

Control where dopants and vacancies are placed:

### Random Placement (Default)

Uniform random distribution - represents high-temperature equilibrium:

```python
gdc_random = builder.create_doped_structure(
    dopant="Gd",
    dopant_fraction=0.20,
    dopant_placement="random",
    vacancy_placement="random",
)
```

### Surface-Preferring Placement

For surface-segregated dopants or vacancies:

```python
gdc_surface = builder.create_doped_structure(
    dopant="Gd",
    dopant_fraction=0.20,
    dopant_placement="surface",
    dopant_preference=0.8,  # Strong surface bias
    vacancy_placement="surface",
    vacancy_preference=0.7,
)
```

### Associated Dopant-Vacancy Pairs

Vacancies near dopants (low-temperature regime):

```python
gdc_associated = builder.create_doped_structure(
    dopant="Gd",
    dopant_fraction=0.20,
    vacancy_placement="near_dopant",
    vacancy_preference=0.9,  # Strong association
)
```

### Bulk-Preferring Placement

For bulk-segregated defects:

```python
gdc_bulk = builder.create_doped_structure(
    dopant="Gd",
    dopant_fraction=0.20,
    dopant_placement="bulk",
    vacancy_placement="bulk",
)
```

## Pr Mixed Valence

Praseodymium is unique - it can be Pr³⁺ or Pr⁴⁺:

```python
# 20% Pr, but only half is Pr³⁺ (contributes to vacancies)
pdc = builder.create_doped_structure(
    dopant="Pr",
    dopant_fraction=0.20,       # 20% total Pr
    pr_trivalent_fraction=0.5,  # 50% is Pr³⁺
    random_seed=42,
)

# Effective: 10% Pr³⁺ → 5% vacancies (half of normal)
```

## Generating Configuration Pools

For screening studies, generate multiple configurations:

```python
# Pool with different vacancy placement strategies
pool = builder.create_doped_pool(
    dopant="Gd",
    dopant_fraction=0.15,
    n_configs_per_strategy=5,
    strategies=["random", "surface", "near_dopant"],
    random_seed=42,
)

print(f"Generated {len(pool)} configurations")  # 15 configs

# Each config has metadata
for config in pool[:3]:
    print(f"Config {config['config_id']}: {config['vacancy_placement']}")
```

## Concentration Series

Generate structures with varying dopant levels:

```python
from nh3sofc.structure import generate_dopant_series

series = generate_dopant_series(
    ceo2_slab,
    dopant="Gd",
    dopant_fractions=[0.05, 0.10, 0.15, 0.20, 0.25],
    n_configs=3,
    random_seed=42,
)

print(f"Generated {len(series)} structures")  # 15 (5 fractions × 3 configs)

# Group by fraction
from collections import defaultdict
by_fraction = defaultdict(list)
for s in series:
    by_fraction[s["dopant_fraction"]].append(s["atoms"])
```

## Analysis

### Single Structure Analysis

```python
from nh3sofc.structure import analyze_dopant_distribution, print_dopant_analysis

# Analyze doped structure
stats = analyze_dopant_distribution(
    gdc,
    dopant="Gd",
    reference_atoms=ceo2_slab,  # Original for vacancy counting
    z_threshold=0.3,            # Top 30% = surface
    near_dopant_cutoff=3.5,     # Å
)

# Print formatted report
print_dopant_analysis(stats, title="GDC Analysis")
```

Output:
```
============================================================
GDC Analysis
============================================================

Dopant: Gd
  Total Gd atoms: 2
  Gd in surface (top 30%): 1 (50.0%)
  Gd in bulk: 1
  Dopant fraction: 25.0%

Host Cations:
  Ce atoms remaining: 6
  O atoms: 15

Vacancy Distribution:
  Total vacancies: 1
  Vacancies in surface: 1 (100.0%)
  Vacancies in bulk: 0
  Vacancies near Gd (< 3.5 Å): 1 (100.0%)

Charge Balance Check:
  Expected vacancies (n_dopant/2): 1
  Actual vacancies: 1
  Status: Charge balanced
============================================================
```

### Key Metrics

| Metric | Description |
|--------|-------------|
| `dopant_surface_fraction` | Fraction of dopants in surface region |
| `vacancy_surface_fraction` | Fraction of vacancies in surface region |
| `vacancy_near_dopant_fraction` | Fraction of vacancies within cutoff of dopants |
| `dopant_fraction` | Actual dopant/(dopant+host) ratio |

## Complete Workflow Example

```python
from ase.build import bulk
from ase.io import write
from nh3sofc.structure import (
    DopantBuilder,
    DopedCeriaStructure,
    analyze_dopant_distribution,
    print_dopant_analysis,
)
from nh3sofc.calculators.vasp import VASPInputGenerator
import os

# 1. Create CeO2 surface
ceo2 = bulk('CeO2', 'fluorite', a=5.41) * (3, 3, 3)

# 2. Create 10% GDC with surface-preferring vacancies
builder = DopantBuilder(ceo2)
gdc = builder.create_doped_structure(
    dopant="Gd",
    dopant_fraction=0.10,
    vacancy_placement="surface",
    vacancy_preference=0.75,
    random_seed=42,
)

# 3. Analyze distribution
stats = analyze_dopant_distribution(gdc, "Gd", reference_atoms=ceo2)
print_dopant_analysis(stats)

# 4. Save structure
write("GDC_10pct.vasp", gdc, format="vasp")

# 5. Generate VASP inputs
os.makedirs("calc_gdc", exist_ok=True)
vasp = VASPInputGenerator(
    gdc,
    calc_type="relax",
    work_dir="calc_gdc",
)
vasp.generate_all(
    encut=520,
    hubbard_u={"Ce": 5.0},  # DFT+U for Ce 4f electrons
    vdw="D3BJ",
)
```

## Supported Dopants

| Dopant | Name | Ionic Radius (Å) | Common Use |
|--------|------|------------------|------------|
| Gd | Gadolinium | 1.053 | GDC - highest conductivity |
| Sm | Samarium | 1.079 | SDC - high conductivity |
| Pr | Praseodymium | 1.126 | PDC - mixed ionic-electronic |
| Y | Yttrium | 1.019 | YDC - stable |
| La | Lanthanum | 1.160 | LDC |
| Nd | Neodymium | 1.109 | NDC |

## Tips

1. **Dopant fraction**: Typical SOFC concentrations are 10-20% (x = 0.10-0.20)

2. **Charge balance**: Always verify with `analyze_dopant_distribution()`

3. **Surface vs bulk**: Use `z_threshold` to define what counts as "surface"

4. **Multiple configs**: Generate pools for statistical sampling of configurations

5. **Pr doping**: Remember to set `pr_trivalent_fraction` for mixed-valence studies

## Next Steps

- [VASP Calculations](vasp_calculations.md) - Run DFT on doped structures
- [High-Throughput Screening](screening.md) - Screen multiple configurations
- [Exsolution Simulation](exsolution.md) - Metal particles on doped ceria
