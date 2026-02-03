# Tutorial: Adsorbate Placement

This tutorial covers 6 different methods for placing adsorbate molecules on surfaces.

## Learning Objectives

- Understand the 6 adsorbate placement methods
- Generate multiple adsorption configurations with proper rotation sampling
- Filter and select configurations for calculations

## The 6 Placement Methods

| Method | Use Case |
|--------|----------|
| `add_simple()` | Manual placement at known position |
| `add_random()` | Exploratory sampling |
| `add_grid()` | Systematic coverage |
| `add_on_site()` | Specific atomic sites (most physically meaningful) |
| `add_with_collision()` | Safe random placement |
| `add_catkit()` | Automated site detection via CatKit |

## Setup

```python
from ase.io import read
from nh3sofc.structure import AdsorbatePlacer

# Load your surface (POSCAR or CIF format recommended)
surface = read("work/surfaces/LaVO3_001/surface.vasp", format="vasp")

# Create placer
placer = AdsorbatePlacer(surface)
```

## Rotation Sampling

All placement methods (except `add_simple()`) apply **uniform random rotations on SO(3)**
for unbiased molecular orientations. This uses proper spherical sampling:

```python
# Internally, rotations use: β = arccos(1 - 2u) for uniform coverage
# This avoids clustering near the "poles" that occurs with naive sampling
```

For reproducibility, always set `random_seed`:

```python
configs_a = placer.add_random("NH3", n_configs=10, random_seed=42)
configs_b = placer.add_random("NH3", n_configs=10, random_seed=42)
# configs_a and configs_b are identical
```

## Method 1: Simple (Manual)

Place at specific (x, y) coordinates with optional rotation.

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

Generate random configurations with uniform rotation sampling.

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

Systematic grid-based placement with multiple orientations per point.

```python
# 3x3 grid with 4 orientations per point
configs = placer.add_grid(
    adsorbate="NH3",
    grid_size=(3, 3),
    orientations=4,   # 4 random orientations per grid point
    height=2.0,
    random_seed=42
)
# Returns 3 × 3 × 4 = 36 configurations
```

## Method 4: On-Site (Atomic Sites)

Place on specific atomic sites. This is the **most physically meaningful** method.

**Important:** Only atoms in the **top bilayer** (within `layer_tolerance` of the
topmost atom) are considered as surface sites. Bulk atoms are automatically excluded.

```python
# On top of La and V atoms in the surface bilayer
configs = placer.add_on_site(
    adsorbate="NH3",
    site_type="ontop",
    atom_types=["La", "V"],
    n_orientations=5,     # 5 orientations per site
    height=2.0,
    layer_tolerance=2.0,  # Top 2 Å = surface bilayer (default)
    random_seed=42
)

# For only the topmost layer (e.g., VO2 termination only)
configs = placer.add_on_site(
    adsorbate="NH3",
    site_type="ontop",
    atom_types=["V"],
    n_orientations=3,
    layer_tolerance=1.0,  # Only topmost ~1 Å
    random_seed=42
)

# Bridge sites (requires CatKit or uses geometric approximation)
configs = placer.add_on_site(
    adsorbate="H",
    site_type="bridge",
    atom_types=["O"],
    height=1.0
)
```

### Surface Detection

The `layer_tolerance` parameter controls which atoms are considered "on the surface":

| `layer_tolerance` | Effect |
|-------------------|--------|
| 1.0 Å | Topmost atomic layer only |
| 2.0 Å (default) | Top bilayer (e.g., VO2 + LaO in perovskites) |
| 3.0 Å | Top 2-3 layers |

## Method 5: Collision-Aware

Random placement with minimum distance checking.

```python
# Ensure no atoms closer than 2.0 Å
configs = placer.add_with_collision(
    adsorbate="NH3",
    n_configs=10,
    min_distance=2.0,
    height=2.0,
    random_seed=42
)
```

## Method 6: CatKit Integration

Use CatKit for automated site detection with rotation support.

```python
# Requires catkit installation: pip install catkit
configs = placer.add_catkit(
    adsorbate="NH3",
    site_type="ontop",      # or "bridge", "hollow", "4fold"
    n_orientations=3,       # 3 orientations per site
    random_seed=42
)
```

## Supported Adsorbates

| Name | Formula | Default Height |
|------|---------|----------------|
| `"NH3"` | NH₃ | 2.0 Å |
| `"NH2"` | NH₂ | 1.8 Å |
| `"NH"` | NH | 1.5 Å |
| `"N"` | N | 1.2 Å |
| `"H"` | H | 1.0 Å |
| `"H2"` | H₂ | 2.5 Å |
| `"H2O"` | H₂O | 2.0 Å |
| `"O"` | O | 1.2 Å |
| `"OH"` | OH | 1.5 Å |

## Filtering Configurations

Remove duplicate configurations based on RMSD:

```python
from nh3sofc.structure.adsorbates import filter_unique_configs, save_configs

# Remove similar configurations (RMSD < 0.5 Å)
unique_configs = filter_unique_configs(configs, threshold=0.5)
print(f"Filtered {len(configs)} → {len(unique_configs)} unique")

# Save to files
paths = save_configs(unique_configs, output_dir="./configs", format="vasp")
```

## Getting Site Information

Inspect available adsorption sites before placement:

```python
info = placer.get_site_info(layer_tolerance=2.0)

print(f"Surface atoms: {info['n_surface_atoms']}")
print(f"By element: {info['element_counts']}")
# Output: {'La': 1, 'V': 1, 'O': 3}
```

## Complete Example

```python
from ase.io import read
from nh3sofc.structure import AdsorbatePlacer
from nh3sofc.structure.adsorbates import filter_unique_configs, save_configs

# Load surface
surface = read("work/surfaces/LaVO3_001/surface.vasp", format="vasp")
placer = AdsorbatePlacer(surface)

# Check available sites
info = placer.get_site_info()
print(f"Surface sites: {info['element_counts']}")

# Generate configurations on La and V sites
configs = placer.add_on_site(
    "NH3",
    atom_types=["La", "V"],
    n_orientations=5,
    random_seed=42
)
print(f"Generated: {len(configs)} configurations")

# Filter duplicates
unique = filter_unique_configs(configs, threshold=0.5)
print(f"Unique: {len(unique)} configurations")

# Save for VASP calculations
paths = save_configs(unique, "work/adsorbates/NH3_on_LaVO3", format="vasp")
```

## Next Steps

- [VASP Calculations](vasp_calculations.md) - Optimize the configurations
- [Decomposition Pathway](decomposition.md) - Study NH3 decomposition
