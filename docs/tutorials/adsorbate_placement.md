# Tutorial: Adsorbate Placement

This tutorial covers 6 different methods for placing adsorbate molecules on surfaces.

## Learning Objectives

- Understand the 6 adsorbate placement methods
- Generate multiple adsorption configurations
- Filter and select configurations for calculations

## The 6 Placement Methods

| Method | Use Case |
|--------|----------|
| `add_simple()` | Manual placement at known position |
| `add_random()` | Exploratory sampling |
| `add_grid()` | Systematic coverage |
| `add_on_site()` | Specific atomic sites |
| `add_with_collision()` | Safe random placement |
| `add_catkit()` | Automated site detection |

## Setup

```python
from ase.io import read
from nh3sofc.structure import AdsorbatePlacer

# Load your surface
surface = read("relaxed_surface.xyz")

# Create placer
placer = AdsorbatePlacer(surface)
```

## Method 1: Simple (Manual)

Place at specific (x, y) coordinates.

```python
# Place NH3 at position (2.5, 2.5) Å, 2.0 Å above surface
result = placer.add_simple(
    adsorbate="NH3",
    position=(2.5, 2.5),
    height=2.0,
    rotation=(0, 0, 0)  # Euler angles (alpha, beta, gamma) in radians
)
```

## Method 2: Random

Generate random configurations with rotation.

```python
# Generate 10 random configurations
configs = placer.add_random(
    adsorbate="NH3",
    n_configs=10,
    height=2.0,
    random_seed=42  # Reproducibility
)

print(f"Generated {len(configs)} configurations")
```

## Method 3: Grid

Systematic grid-based placement.

```python
# 3x3 grid of positions
configs = placer.add_grid(
    adsorbate="NH3",
    grid_size=(3, 3),
    orientations=1,  # One orientation per grid point
    height=2.0
)
# Returns 9 configurations
```

## Method 4: On-Site (Atomic Sites)

Place on specific atomic sites (ontop, bridge, hollow).

```python
# On top of La atoms only
configs = placer.add_on_site(
    adsorbate="NH3",
    site_type="ontop",
    atom_types=["La"],
    height=2.0
)

# On top of all metal atoms
configs = placer.add_on_site(
    adsorbate="NH3",
    site_type="ontop",
    atom_types=["La", "V"],
    height=2.0
)

# Bridge sites between O atoms
configs = placer.add_on_site(
    adsorbate="H",
    site_type="bridge",
    atom_types=["O"],
    height=1.0
)
```

## Method 5: Collision-Aware

Random placement with minimum distance checking.

```python
# Ensure no atoms closer than 2.0 Å
configs = placer.add_with_collision(
    adsorbate="NH3",
    n_configs=10,
    min_distance=2.0,
    height=2.0
)
```

## Method 6: CatKit Integration

Use CatKit for automated site detection.

```python
# Requires catkit installation
configs = placer.add_catkit(
    adsorbate="NH3",
    site_type="ontop"  # or "bridge", "hollow"
)
```

## Supported Adsorbates

| Name | Formula | Default Height |
|------|---------|----------------|
| `"NH3"` | NH3 | 2.0 Å |
| `"NH2"` | NH2 | 1.8 Å |
| `"NH"` | NH | 1.5 Å |
| `"N"` | N | 1.2 Å |
| `"H"` | H | 1.0 Å |
| `"H2"` | H2 | 2.5 Å |
| `"H2O"` | H2O | 2.0 Å |
| `"O"` | O | 1.2 Å |
| `"OH"` | OH | 1.5 Å |

## Filtering Configurations

```python
from nh3sofc.structure import filter_by_rmsd

# Remove similar configurations
unique_configs = filter_by_rmsd(
    configs,
    threshold=0.5  # Å
)

print(f"Unique configurations: {len(unique_configs)}")
```

## Complete Example

```python
from ase.io import read, write
from nh3sofc.structure import AdsorbatePlacer

# Load surface
surface = read("surface.xyz")

# Create placer
placer = AdsorbatePlacer(surface)

# Method 1: On top of La atoms
la_configs = placer.add_on_site("NH3", site_type="ontop", atom_types=["La"])

# Method 2: On top of V atoms
v_configs = placer.add_on_site("NH3", site_type="ontop", atom_types=["V"])

# Combine all configs
all_configs = la_configs + v_configs

# Save for calculations
for i, config in enumerate(all_configs):
    write(f"nh3_config_{i:03d}.xyz", config)

print(f"Generated {len(all_configs)} configurations")
```

## Next Steps

- [VASP Calculations](vasp_calculations.md) - Optimize the configurations
- [Decomposition Pathway](decomposition.md) - Study NH3 decomposition
