# Doped Fluorite Materials Tutorial (GDC, SDC, YSZ, ScSZ)

This tutorial covers creating acceptor-doped fluorite oxide structures for solid oxide fuel cell applications.

## Overview

Doped fluorite oxides are crucial for ionic conductivity in SOFCs:

### Ceria-Based (CeO2)

- **GDC**: Gadolinium-doped ceria (Ce₁₋ₓGdₓO₂₋δ)
- **SDC**: Samarium-doped ceria (Ce₁₋ₓSmₓO₂₋δ)
- **PDC**: Praseodymium-doped ceria (Ce₁₋ₓPrₓO₂₋δ)
- **YDC**: Yttrium-doped ceria (Ce₁₋ₓYₓO₂₋δ)

### Zirconia-Based (ZrO2)

- **YSZ**: Yttria-stabilized zirconia (Zr₁₋ₓYₓO₂₋δ)
- **ScSZ**: Scandia-stabilized zirconia (Zr₁₋ₓScₓO₂₋δ)
- **CSZ**: Calcia-stabilized zirconia (Zr₁₋ₓCaₓO₂₋δ)
- **MSZ**: Magnesia-stabilized zirconia (Zr₁₋ₓMgₓO₂₋δ)

### Chemistry

Acceptor dopants substitute M⁴⁺ sites (Ce⁴⁺, Zr⁴⁺). Charge compensation requires oxygen vacancy formation:

**Trivalent dopants (Gd³⁺, Y³⁺, Sc³⁺):**
```
2 M⁴⁺ → 2 D³⁺ + V_O^{2+}  (2:1 dopant:vacancy ratio)
```

**Divalent dopants (Ca²⁺, Mg²⁺):**
```
M⁴⁺ → D²⁺ + V_O^{2+}  (1:1 dopant:vacancy ratio)
```

!!! note "Odd Number of Trivalent Dopants"
    When the number of trivalent dopants is odd (e.g., 9 Gd³⁺), the ideal vacancy count would be
    non-integer (4.5). The code uses floor division, creating 4 vacancies in this case.

    The uncompensated negative charge physically corresponds to M⁴⁺ → M³⁺ reduction
    (small polaron formation), which is realistic for reducible oxides like CeO₂.

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

### Creating Yttria-Stabilized Zirconia (YSZ)

```python
from ase.build import bulk
from nh3sofc.structure import DopantBuilder

# Create ZrO2 supercell (cubic fluorite phase)
zro2 = bulk('ZrO2', 'fluorite', a=5.07)
zro2_slab = zro2 * (2, 2, 2)  # 8 Zr, 16 O

# Initialize builder
builder = DopantBuilder(zro2_slab)

# Create 8% Y-doped zirconia (8YSZ)
ysz = builder.create_doped_structure(
    dopant="Y",
    dopant_fraction=0.08,
    host_cation="Zr",  # Specify the host cation
    random_seed=42,
)

# Check composition
symbols = ysz.get_chemical_symbols()
print(f"Y atoms: {symbols.count('Y')}")
print(f"Zr atoms: {symbols.count('Zr')}")
```

### Creating Scandia-Stabilized Zirconia (ScSZ)

```python
# ScSZ has the highest ionic conductivity among zirconia-based electrolytes
scsz = builder.create_doped_structure(
    dopant="Sc",
    dopant_fraction=0.10,  # 10 mol% Sc2O3
    host_cation="Zr",
    random_seed=42,
)
```

### Creating Calcia-Stabilized Zirconia (CSZ)

```python
# Divalent dopants have 1:1 dopant:vacancy ratio
csz = builder.create_doped_structure(
    dopant="Ca",
    dopant_fraction=0.15,
    host_cation="Zr",
    random_seed=42,
)

# 1 Ca²⁺ → 1 vacancy (not 2:1 like trivalent dopants)
```

### Using DopedFluoriteStructure

The `DopedFluoriteStructure` class tracks dopant metadata:

```python
from nh3sofc.structure import DopedFluoriteStructure

# Ceria-based (GDC)
gdc = DopedFluoriteStructure.from_ceria(
    ceo2_slab,
    dopant="Gd",
    dopant_fraction=0.10,
    random_seed=42,
)
print(gdc)  # DopedFluoriteStructure(GDC, Gd_frac=10.0%, n_vac=0)
print(gdc.get_material_name())  # "GDC"

# Zirconia-based (YSZ)
ysz = DopedFluoriteStructure.from_zirconia(
    zro2_slab,
    dopant="Y",
    dopant_fraction=0.08,
    random_seed=42,
)
print(ysz.get_material_name())  # "YSZ"
print(ysz.host_cation)  # "Zr"

# General factory (works with any host)
doped = DopedFluoriteStructure.from_parent(
    parent_structure,
    dopant="Sc",
    dopant_fraction=0.10,
    host_cation="Zr",
    random_seed=42,
)
```

!!! tip "Backwards Compatibility"
    `DopedCeriaStructure` is still available as a deprecated alias.
    Existing code will continue to work but will show deprecation warnings.

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

## Mixed Valence Dopants (Pr, Tb)

Praseodymium and Terbium can exist as M³⁺ or M⁴⁺:

```python
# 20% Pr, but only half is Pr³⁺ (contributes to vacancies)
pdc = builder.create_doped_structure(
    dopant="Pr",
    dopant_fraction=0.20,       # 20% total Pr
    pr_trivalent_fraction=0.5,  # 50% is Pr³⁺
    random_seed=42,
)
# Effective: 10% Pr³⁺ → 5% vacancies (half of normal)

# Same applies to Tb
tdc = builder.create_doped_structure(
    dopant="Tb",
    dopant_fraction=0.15,
    pr_trivalent_fraction=0.7,  # 70% is Tb³⁺
    random_seed=42,
)
```

!!! note "Parameter Name"
    The parameter is named `pr_trivalent_fraction` for historical reasons but
    applies to both Pr and Tb dopants. The alias `trivalent_fraction` is also
    available in `DopedFluoriteStructure.from_ceria()`.

## Generating Configuration Pools

For screening studies, generate multiple configurations:

```python
# Pool with different vacancy placement strategies
pool = builder.create_doped_pool(
    dopant="Gd",
    dopant_fraction=0.15,
    n_configs=5,
    strategies=["random", "surface", "near_dopant"],
    random_seed=42,
)

print(f"Generated {len(pool)} configurations")  # 15 configs

# Each config has metadata
for config in pool[:3]:
    print(f"Config {config['config_id']}: {config['vacancy_placement']}")
```

### Zirconia Pools

```python
# YSZ configuration pool
ysz_pool = builder.create_doped_pool(
    dopant="Y",
    dopant_fraction=0.08,
    n_configs=3,
    strategies=["random", "surface"],
    host_cation="Zr",
    random_seed=42,
)
```

## Concentration Series

Generate structures with varying dopant levels by passing a list of fractions:

```python
# Use create_doped_pool with a list of fractions
pool = builder.create_doped_pool(
    dopant="Gd",
    dopant_fraction=[0.05, 0.10, 0.15, 0.20, 0.25],
    n_configs=3,
    random_seed=42,
)

print(f"Generated {len(pool)} structures")  # 15 (5 fractions × 3 configs)

# Group by fraction
from collections import defaultdict
by_fraction = defaultdict(list)
for s in pool:
    by_fraction[s["dopant_fraction"]].append(s["atoms"])
```

!!! tip "Default Strategy for Concentration Series"
    When `dopant_fraction` is a list, the default strategy is `["random"]` since
    concentration series typically don't need strategy variation. Override with
    `strategies=["random", "surface"]` if needed.

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
    DopedFluoriteStructure,
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

### YSZ Workflow Example

```python
# 1. Create ZrO2 structure
zro2 = bulk('ZrO2', 'fluorite', a=5.07) * (3, 3, 3)

# 2. Create 8YSZ
ysz = DopedFluoriteStructure.from_zirconia(
    zro2,
    dopant="Y",
    dopant_fraction=0.08,
    vacancy_placement="random",
    random_seed=42,
)

print(f"Created {ysz.get_material_name()}")
print(f"Host: {ysz.host_cation}, Dopant: {ysz.dopant}")
print(f"Vacancies: {ysz.n_vacancies}")

# Save structure
write("YSZ_8mol.vasp", ysz.atoms, format="vasp")
```

## Supported Dopants

### Trivalent Dopants (2:1 dopant:vacancy ratio)

| Dopant | Name | Ionic Radius (Å) | Common Use |
|--------|------|------------------|------------|
| Gd | Gadolinium | 1.053 | GDC - highest conductivity for ceria |
| Sm | Samarium | 1.079 | SDC - high conductivity |
| Pr | Praseodymium | 1.126 | PDC - mixed ionic-electronic |
| Y | Yttrium | 1.019 | YDC, YSZ - most common for zirconia |
| La | Lanthanum | 1.160 | LDC |
| Nd | Neodymium | 1.109 | NDC |
| Tb | Terbium | 1.040 | TDC - mixed valence |
| Sc | Scandium | 0.870 | ScSZ - highest conductivity for zirconia |

### Divalent Dopants (1:1 dopant:vacancy ratio)

| Dopant | Name | Ionic Radius (Å) | Common Use |
|--------|------|------------------|------------|
| Ca | Calcium | 1.120 | CSZ - cost-effective |
| Mg | Magnesium | 0.890 | MSZ |

## Supported Host Materials

| Host | Oxide | Lattice (Å) | Notes |
|------|-------|-------------|-------|
| Ce | CeO2 | 5.41 | Reducible, high electronic conductivity at low pO2 |
| Zr | ZrO2 | 5.07 | Stable, purely ionic conductor |
| Hf | HfO2 | 5.11 | Similar to zirconia |
| Th | ThO2 | 5.60 | Nuclear applications |

## Tips

1. **Dopant fraction**: Typical SOFC concentrations are 8-20% (x = 0.08-0.20)

2. **Charge balance**: Always verify with `analyze_dopant_distribution()`

3. **Surface vs bulk**: Use `z_threshold` to define what counts as "surface"

4. **Multiple configs**: Generate pools for statistical sampling of configurations

5. **Mixed valence**: Remember to set `pr_trivalent_fraction` for Pr/Tb studies

6. **YSZ vs ScSZ**: ScSZ has higher conductivity but Sc is more expensive

7. **Divalent dopants**: Remember they create more vacancies per dopant (1:1 vs 2:1)

## Next Steps

- [VASP Calculations](vasp_calculations.md) - Run DFT on doped structures
- [High-Throughput Screening](screening.md) - Screen multiple configurations
- [Exsolution Simulation](exsolution.md) - Metal particles on doped fluorites
