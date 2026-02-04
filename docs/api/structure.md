# Structure Module

The `nh3sofc.structure` module provides tools for building and manipulating atomic structures.

## BulkStructure

Load and manipulate bulk crystal structures.

```python
from nh3sofc.structure import BulkStructure
```

### Constructor

```python
BulkStructure(atoms: Atoms)
```

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `atoms` | `Atoms` | ASE Atoms object |

### Class Methods

#### from_cif

```python
@classmethod
BulkStructure.from_cif(
    filepath: str,
    primitive: bool = False
) -> BulkStructure
```

Load structure from CIF file.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `filepath` | `str` | - | Path to CIF file |
| `primitive` | `bool` | `False` | Convert to primitive cell |

**Example:**

```python
bulk = BulkStructure.from_cif("LaVO3.cif")
print(f"Formula: {bulk.atoms.get_chemical_formula()}")
```

### Methods

#### make_supercell

```python
make_supercell(size: tuple) -> BulkStructure
```

Create supercell.

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `size` | `tuple` | Supercell size (nx, ny, nz) |

**Example:**

```python
supercell = bulk.make_supercell((2, 2, 2))
```

#### get_spacegroup

```python
get_spacegroup(symprec: float = 1e-5) -> str
```

Get space group symbol.

---

## SurfaceBuilder

Create surface slabs from bulk structures.

```python
from nh3sofc.structure import SurfaceBuilder
```

### Constructor

```python
SurfaceBuilder(bulk_structure: Union[BulkStructure, Atoms])
```

### Methods

#### create_surface

```python
create_surface(
    miller_index: tuple,
    layers: int = 6,
    vacuum: float = 15.0,
    fix_bottom: int = 0
) -> SlabStructure
```

Create surface slab.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `miller_index` | `tuple` | - | Miller indices (h, k, l) |
| `layers` | `int` | `6` | Number of atomic layers |
| `vacuum` | `float` | `15.0` | Vacuum thickness (Å) |
| `fix_bottom` | `int` | `0` | Layers to fix at bottom |

**Example:**

```python
builder = SurfaceBuilder(bulk)
surface = builder.create_surface(
    miller_index=(0, 0, 1),
    layers=6,
    vacuum=15.0,
    fix_bottom=2
)
```

#### get_all_terminations

```python
get_all_terminations(miller_index: tuple) -> List[Atoms]
```

Get all possible terminations for a surface.

#### create_symmetric_slab

```python
create_symmetric_slab(
    miller_index: tuple,
    layers: int = 7,
    vacuum: float = 15.0,
    fix_bottom: int = 2
) -> SlabStructure
```

Create symmetric slab with same termination on both surfaces (helps cancel dipole).

---

## SlabStructure

Surface slab structure with polarity and layer analysis methods.

```python
from nh3sofc.structure import SlabStructure
```

### Polarity Methods

#### check_polarity

```python
check_polarity(
    formal_charges: dict = None,
    tolerance: float = 0.5
) -> dict
```

Check if surface is polar.

**Returns:**

```python
{
    "is_polar": bool,
    "dipole_moment": float,  # e·Å
    "dipole_z": float,       # z-component
    "top_layer_charge": float,
    "bottom_layer_charge": float,
    "layer_charges": list,
    "recommendation": str
}
```

#### calculate_dipole_moment

```python
calculate_dipole_moment(formal_charges: dict = None) -> tuple
```

Calculate electric dipole moment. Returns `(magnitude, direction_vector)`.

### Layer Analysis Methods

#### identify_layers

```python
identify_layers(tolerance: float = 0.5) -> List[dict]
```

Identify atomic layers. Returns list of layer info dicts.

#### get_layer_spacing

```python
get_layer_spacing() -> List[float]
```

Get interlayer distances in Angstrom.

### Stoichiometry Methods

#### check_stoichiometry

```python
check_stoichiometry(
    expected: dict = None,
    tolerance: float = 0.1
) -> dict
```

Check if stoichiometry matches expected.

---

## DefectBuilder

Create point defects and oxynitride structures.

```python
from nh3sofc.structure import DefectBuilder
```

### Constructor

```python
DefectBuilder(structure: Union[SlabStructure, Atoms])
```

### Methods

#### create_oxynitride

```python
create_oxynitride(
    nitrogen_fraction: float = 0.67,
    vacancy_concentration: float = 0.0,
    vacancy_element: str = "N",
    placement: str = "random",
    random_seed: int = None
) -> Atoms
```

Create oxynitride by substituting O with N and creating vacancies.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `nitrogen_fraction` | `float` | `0.67` | Fraction of O to replace with N |
| `vacancy_concentration` | `float` | `0.0` | Fraction of anion sites to leave vacant |
| `vacancy_element` | `str` | `"N"` | Element to create vacancies in ("N" or "O") |
| `placement` | `str` | `"random"` | Defect placement strategy (see below) |
| `random_seed` | `int` | `None` | Random seed for reproducibility |

**Placement strategies:**

| Strategy | Description |
|----------|-------------|
| `"random"` | Random placement throughout structure |
| `"surface"` | Preferentially place defects near surface (high z) |
| `"near_N"` | Place vacancies near existing N atoms |

**Example:**

```python
defect = DefectBuilder(surface)

# Random placement (default)
oxynitride = defect.create_oxynitride(
    nitrogen_fraction=0.67,
    vacancy_concentration=0.10
)

# Surface-populated defects
oxynitride_surface = defect.create_oxynitride(
    nitrogen_fraction=0.67,
    vacancy_concentration=0.10,
    placement="surface"
)
```

#### create_oxynitride_pool

```python
create_oxynitride_pool(
    nitrogen_fraction: float = 0.67,
    vacancy_concentration: float = 0.0,
    vacancy_element: str = "N",
    n_configs_per_strategy: int = 3,
    strategies: List[str] = None,
    random_seed: int = None
) -> List[Dict]
```

Generate a pool of oxynitride configurations with different placement strategies.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `nitrogen_fraction` | `float` | `0.67` | Fraction of O to replace with N |
| `vacancy_concentration` | `float` | `0.0` | Fraction of anion sites to leave vacant |
| `vacancy_element` | `str` | `"N"` | Element to create vacancies in |
| `n_configs_per_strategy` | `int` | `3` | Configurations per strategy |
| `strategies` | `List[str]` | `["random", "surface", "near_N"]` | Strategies to use |
| `random_seed` | `int` | `None` | Base random seed |

**Returns:**

```python
[{
    "atoms": Atoms,
    "placement": str,
    "nitrogen_fraction": float,
    "vacancy_concentration": float,
    "config_id": int,
}, ...]
```

**Example:**

```python
defect = DefectBuilder(surface)
pool = defect.create_oxynitride_pool(
    nitrogen_fraction=0.67,
    vacancy_concentration=0.1,
    n_configs_per_strategy=5
)
print(f"Generated {len(pool)} configurations")  # 15 configs (5 x 3 strategies)
```

#### create_vacancy

```python
create_vacancy(
    element: str,
    concentration: float,
    indices: List[int] = None,
    random_seed: int = None
) -> Atoms
```

Create vacancies by removing atoms.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `element` | `str` | - | Element to remove (e.g., "O") |
| `concentration` | `float` | - | Fraction of atoms to remove (0.0 to 1.0) |
| `indices` | `List[int]` | `None` | Specific indices to remove (if None, random selection) |
| `random_seed` | `int` | `None` | Random seed for reproducibility |

---

## AdsorbatePlacer

Place adsorbate molecules on surfaces using 6 methods.

```python
from nh3sofc.structure import AdsorbatePlacer
```

### Constructor

```python
AdsorbatePlacer(
    slab: Atoms,
    adsorbate_info: dict = None
)
```

### Methods

#### add_simple

```python
add_simple(
    adsorbate: str,
    position: tuple,
    height: float = 2.0,
    orientation: tuple = None
) -> Atoms
```

Manual placement at specific (x, y) position.

**Example:**

```python
placer = AdsorbatePlacer(slab)
result = placer.add_simple("NH3", position=(2.5, 2.5), height=2.0)
```

#### add_random

```python
add_random(
    adsorbate: str,
    n_configs: int = 10,
    height: float = 2.0,
    seed: int = None
) -> List[Atoms]
```

Random placement with rotation.

#### add_grid

```python
add_grid(
    adsorbate: str,
    nx: int = 3,
    ny: int = 3,
    height: float = 2.0
) -> List[Atoms]
```

Grid-based systematic placement.

#### add_on_site

```python
add_on_site(
    adsorbate: str,
    site_type: str = "ontop",
    atom_types: List[str] = None,
    height: float = 2.0
) -> List[Atoms]
```

Place on specific atomic sites.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `adsorbate` | `str` | - | Molecule name ("NH3", "H2O", etc.) |
| `site_type` | `str` | `"ontop"` | Site type: "ontop", "bridge", "hollow" |
| `atom_types` | `List[str]` | `None` | Filter by atom types |
| `height` | `float` | `2.0` | Height above surface (Å) |

**Example:**

```python
# Place NH3 on top of La and V atoms only
configs = placer.add_on_site(
    "NH3",
    site_type="ontop",
    atom_types=["La", "V"]
)
```

#### add_with_collision

```python
add_with_collision(
    adsorbate: str,
    n_configs: int = 10,
    min_distance: float = 2.0,
    height: float = 2.0
) -> List[Atoms]
```

Random placement with collision detection.

#### add_catkit

```python
add_catkit(
    adsorbate: str,
    site_type: str = "ontop"
) -> List[Atoms]
```

Use CatKit for automated site detection.

---

## DecompositionBuilder

Generate NH3 decomposition intermediate configurations.

```python
from nh3sofc.structure import DecompositionBuilder
```

### Constructor

```python
DecompositionBuilder(
    nh3_on_slab: Atoms,
    random_seed: int = None
)
```

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `nh3_on_slab` | `Atoms` | Optimized NH3 on surface structure |
| `random_seed` | `int` | Random seed for reproducibility |

### Methods

#### create_NH2_H_configs

```python
create_NH2_H_configs(
    n_configs: int = 10,
    h_height: float = 1.0,
    min_h_distance: float = 1.5,
    max_h_distance: float = 4.0
) -> List[Atoms]
```

Generate NH2* + H* configurations.

**Example:**

```python
decomp = DecompositionBuilder(nh3_on_slab)
nh2_h = decomp.create_NH2_H_configs(n_configs=10)
```

#### create_NH_2H_configs

```python
create_NH_2H_configs(
    n_configs: int = 10,
    h_h_min: float = 1.5
) -> List[Atoms]
```

Generate NH* + 2H* configurations.

#### create_N_3H_configs

```python
create_N_3H_configs(
    n_configs: int = 5,
    cluster_h: bool = False
) -> List[Atoms]
```

Generate N* + 3H* configurations.

### Helper Functions

#### generate_decomposition_pathway

```python
generate_decomposition_pathway(
    nh3_on_slab: Atoms,
    n_configs_per_step: int = 5,
    random_seed: int = None
) -> Dict[str, List[Atoms]]
```

Generate all decomposition intermediates.

**Returns:**

```python
{
    "NH3": [atoms],
    "NH2_H": [atoms, atoms, ...],
    "NH_2H": [atoms, atoms, ...],
    "N_3H": [atoms, atoms, ...]
}
```

---

## Supported Adsorbates

| Name | Formula | Atoms | Default Height |
|------|---------|-------|---------------|
| `"NH3"` | NH3 | 4 | 2.0 Å |
| `"NH2"` | NH2 | 3 | 1.8 Å |
| `"NH"` | NH | 2 | 1.5 Å |
| `"N"` | N | 1 | 1.2 Å |
| `"H"` | H | 1 | 1.0 Å |
| `"H2"` | H2 | 2 | 2.5 Å |
| `"H2O"` | H2O | 3 | 2.0 Å |
| `"O"` | O | 1 | 1.2 Å |
| `"OH"` | OH | 2 | 1.5 Å |

---

## ExsolutionBuilder

Create structures for exsolution processes where transition metals migrate from perovskite bulk to surface.

```python
from nh3sofc.structure import ExsolutionBuilder
```

### Constructor

```python
ExsolutionBuilder(
    structure: Union[SlabStructure, Atoms],
    a_site_elements: List[str] = None,
    b_site_elements: List[str] = None
)
```

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `structure` | `SlabStructure` or `Atoms` | - | Perovskite surface slab |
| `a_site_elements` | `List[str]` | `["La", "Sr", "Ba", "Ca"]` | A-site elements |
| `b_site_elements` | `List[str]` | `["Ti", "V", "Mn", "Fe", "Co", "Ni"]` | B-site elements |

### Methods

#### identify_perovskite_sites

```python
identify_perovskite_sites() -> Dict[str, List[int]]
```

Identify A-site, B-site, and O-site atom indices.

**Returns:**

```python
{
    "A_site": [0, 1, 2, ...],
    "B_site": [4, 5, ...],
    "O_site": [8, 9, 10, ...]
}
```

#### create_defective_perovskite

```python
create_defective_perovskite(
    a_site_vacancy_fraction: float = 0.0,
    b_site_vacancy_fraction: float = 0.0,
    oxygen_vacancy_fraction: float = 0.0,
    b_site_element: str = None,
    random_seed: int = None
) -> Atoms
```

Create defective perovskite with cation and anion vacancies.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `a_site_vacancy_fraction` | `float` | `0.0` | Fraction of A-site to remove |
| `b_site_vacancy_fraction` | `float` | `0.0` | Fraction of B-site to remove |
| `oxygen_vacancy_fraction` | `float` | `0.0` | Fraction of O to remove |
| `b_site_element` | `str` | `None` | Specific B-site element for vacancies |
| `random_seed` | `int` | `None` | Random seed for reproducibility |

**Example:**

```python
builder = ExsolutionBuilder(surface)
defective = builder.create_defective_perovskite(
    a_site_vacancy_fraction=0.05,
    b_site_vacancy_fraction=0.1,
    oxygen_vacancy_fraction=0.08,
    b_site_element="Ni"
)
```

#### create_surface_segregation

```python
create_surface_segregation(
    metal: str,
    n_atoms: int = 1,
    random_seed: int = None
) -> Atoms
```

Move B-site metal atoms from bulk to surface layer.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `metal` | `str` | - | Element to segregate ("Ni", "Co", "Fe") |
| `n_atoms` | `int` | `1` | Number of atoms to move |
| `random_seed` | `int` | `None` | Random seed |

#### create_nanoparticle

```python
create_nanoparticle(
    metal: str,
    n_atoms: int,
    shape: str = "hemispherical",
    position: str = "hollow",
    interface_distance: float = 2.0,
    socketed: bool = True,
    random_seed: int = None
) -> Atoms
```

Place metallic nanoparticle on surface.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `metal` | `str` | - | Metal element ("Ni", "Co", "Fe") |
| `n_atoms` | `int` | - | Cluster size (1, 4, 7, 13, 19) |
| `shape` | `str` | `"hemispherical"` | Shape: "hemispherical" or "icosahedral" |
| `position` | `str` | `"hollow"` | Position: "hollow", "ontop", "bridge", "random" |
| `interface_distance` | `float` | `2.0` | Metal-surface distance (Å) |
| `socketed` | `bool` | `True` | Remove atoms under particle (real exsolution) |
| `random_seed` | `int` | `None` | Random seed |

**Example:**

```python
builder = ExsolutionBuilder(surface)
with_particle = builder.create_nanoparticle(
    metal="Ni",
    n_atoms=13,
    shape="hemispherical",
    socketed=True
)
```

#### create_exsolution_pathway

```python
create_exsolution_pathway(
    metal: str,
    particle_size: int,
    vacancy_fraction: float = 0.1,
    random_seed: int = None
) -> List[Dict[str, Any]]
```

Generate complete exsolution pathway structures.

**Returns:**

```python
[
    {"atoms": Atoms, "stage": "pristine", "description": "..."},
    {"atoms": Atoms, "stage": "defective", "description": "..."},
    {"atoms": Atoms, "stage": "segregated", "description": "..."},
    {"atoms": Atoms, "stage": "exsolved", "description": "..."}
]
```

#### get_adsorption_sites

```python
get_adsorption_sites(
    structure: Atoms = None,
    particle_indices: List[int] = None
) -> Dict[str, List[Dict]]
```

Identify adsorption sites on exsolved particle system.

**Returns:**

```python
{
    "metal_top": [...],        # Top of metal particle
    "interface_edge": [...],   # Metal-oxide boundary
    "vacancy_site": [...],     # Near O vacancies
    "oxide_surface": [...]     # Remaining perovskite
}
```

---

## ExsolutionStructure

Structure class with exsolution metadata tracking.

```python
from nh3sofc.structure import ExsolutionStructure
```

### Constructor

```python
ExsolutionStructure(
    atoms: Atoms,
    exsolution_metal: str = None,
    nanoparticle_size: int = None,
    a_site_vacancy_concentration: float = 0.0,
    b_site_vacancy_concentration: float = 0.0,
    oxygen_vacancy_concentration: float = 0.0,
    exsolution_state: str = "pristine"
)
```

### Attributes

| Name | Type | Description |
|------|------|-------------|
| `exsolution_metal` | `str` | Metal being exsolved (Ni, Co, Fe) |
| `nanoparticle_size` | `int` | Number of atoms in particle |
| `exsolution_state` | `str` | State: pristine, defective, segregated, exsolved |
| `oxygen_vacancy_concentration` | `float` | O vacancy fraction |

---

## Helper Functions

### create_metallic_cluster

```python
create_metallic_cluster(
    metal: str,
    n_atoms: int,
    shape: str = "hemispherical"
) -> Atoms
```

Create isolated metallic cluster.

**Example:**

```python
from nh3sofc.structure import create_metallic_cluster

ni13 = create_metallic_cluster("Ni", 13, "hemispherical")
```

### generate_exsolution_series

```python
generate_exsolution_series(
    base_structure: Atoms,
    metal: str,
    vacancy_range: List[float],
    particle_sizes: List[int],
    n_configs: int = 3,
    random_seed: int = None
) -> List[Dict]
```

Generate series for screening.

**Example:**

```python
from nh3sofc.structure import generate_exsolution_series

series = generate_exsolution_series(
    surface,
    metal="Ni",
    vacancy_range=[0.0, 0.05, 0.1],
    particle_sizes=[1, 13],
    n_configs=3
)
```

---

## DopantBuilder

Create acceptor-doped ceria structures (GDC, SDC, PDC) with charge-compensating oxygen vacancies.

```python
from nh3sofc.structure import DopantBuilder
```

### Constructor

```python
DopantBuilder(structure: Union[SlabStructure, Atoms])
```

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `structure` | `SlabStructure` or `Atoms` | CeO2 structure (bulk or slab) to dope |

### Methods

#### create_doped_structure

```python
create_doped_structure(
    dopant: str,
    dopant_fraction: float,
    host_cation: str = "Ce",
    auto_compensate: bool = True,
    vacancy_placement: str = "random",
    dopant_placement: str = "random",
    dopant_preference: float = 0.7,
    vacancy_preference: float = 0.7,
    pr_trivalent_fraction: float = 1.0,
    z_threshold: float = 0.3,
    random_seed: int = None
) -> Atoms
```

Create doped structure with charge-compensating vacancies.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `dopant` | `str` | - | Dopant element ("Sm", "Gd", "Pr", "Y", "La", "Nd") |
| `dopant_fraction` | `float` | - | Fraction of host cations to replace (0.0 to 1.0) |
| `host_cation` | `str` | `"Ce"` | Cation to replace |
| `auto_compensate` | `bool` | `True` | Create compensating vacancies automatically |
| `vacancy_placement` | `str` | `"random"` | Vacancy placement strategy (see below) |
| `dopant_placement` | `str` | `"random"` | Dopant placement strategy |
| `dopant_preference` | `float` | `0.7` | Bias strength for dopant placement (0.5=uniform, 1.0=strong) |
| `vacancy_preference` | `float` | `0.7` | Bias strength for vacancy placement |
| `pr_trivalent_fraction` | `float` | `1.0` | For Pr: fraction that is Pr³⁺ (rest is Pr⁴⁺) |
| `z_threshold` | `float` | `0.3` | Fraction of slab considered "surface" |
| `random_seed` | `int` | `None` | Random seed for reproducibility |

**Placement strategies:**

| Strategy | Description |
|----------|-------------|
| `"random"` | Uniform random distribution (ignores preference) |
| `"surface"` | Preferentially place at high z (surface) |
| `"bulk"` | Preferentially place at low z (bulk) |
| `"near_dopant"` | Place vacancies near dopant atoms (vacancy only) |

**Example:**

```python
builder = DopantBuilder(ceo2_slab)

# 10% Gd-doped ceria (GDC) with random placement
gdc = builder.create_doped_structure(
    dopant="Gd",
    dopant_fraction=0.10,
)

# GDC with vacancies preferring dopant sites
gdc_associated = builder.create_doped_structure(
    dopant="Gd",
    dopant_fraction=0.10,
    vacancy_placement="near_dopant",
    vacancy_preference=0.8,
)

# Pr-doped ceria with mixed valence
pdc = builder.create_doped_structure(
    dopant="Pr",
    dopant_fraction=0.20,
    pr_trivalent_fraction=0.5,  # 50% Pr³⁺, 50% Pr⁴⁺
)
```

#### create_doped_pool

```python
create_doped_pool(
    dopant: str,
    dopant_fraction: float,
    n_configs_per_strategy: int = 3,
    strategies: List[str] = None,
    host_cation: str = "Ce",
    dopant_placement: str = "random",
    dopant_preference: float = 0.7,
    vacancy_preference: float = 0.7,
    pr_trivalent_fraction: float = 1.0,
    z_threshold: float = 0.3,
    random_seed: int = None
) -> List[Dict]
```

Generate multiple configurations with different placement strategies.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `dopant` | `str` | - | Dopant element |
| `dopant_fraction` | `float` | - | Fraction of host cations to replace |
| `n_configs_per_strategy` | `int` | `3` | Configurations per strategy |
| `strategies` | `List[str]` | `["random", "surface", "near_dopant"]` | Strategies to use |

**Returns:**

```python
[{
    "atoms": Atoms,
    "dopant": str,
    "dopant_fraction": float,
    "vacancy_placement": str,
    "config_id": int,
}, ...]
```

**Example:**

```python
pool = builder.create_doped_pool(
    dopant="Sm",
    dopant_fraction=0.15,
    n_configs_per_strategy=5,
    strategies=["random", "surface", "near_dopant"],
)
print(f"Generated {len(pool)} configurations")  # 15 configs
```

---

## DopedCeriaStructure

Structure class for doped ceria materials with metadata tracking.

```python
from nh3sofc.structure import DopedCeriaStructure
```

### Constructor

```python
DopedCeriaStructure(
    atoms: Atoms,
    dopant: str,
    dopant_fraction: float,
    n_vacancies: int,
    parent_formula: str = None
)
```

### Class Methods

#### from_ceria

```python
@classmethod
DopedCeriaStructure.from_ceria(
    ceria: Union[Atoms, BaseStructure],
    dopant: str,
    dopant_fraction: float,
    host_cation: str = "Ce",
    vacancy_placement: str = "random",
    dopant_placement: str = "random",
    random_seed: int = None
) -> DopedCeriaStructure
```

Create doped ceria from pristine CeO2 structure.

**Example:**

```python
doped = DopedCeriaStructure.from_ceria(
    ceo2_slab,
    dopant="Gd",
    dopant_fraction=0.10,
)
print(doped.get_dopant_name())  # "GDC"
```

### Methods

#### get_stoichiometry

```python
get_stoichiometry() -> Dict[str, float]
```

Calculate stoichiometry relative to total cation sites.

**Example:**

```python
stoich = doped.get_stoichiometry()
# {"Ce": 0.9, "Gd": 0.1, "O": 1.95}
```

#### get_dopant_name

```python
get_dopant_name() -> str
```

Get material abbreviation (GDC, SDC, PDC, etc.).

### Attributes

| Name | Type | Description |
|------|------|-------------|
| `dopant` | `str` | Dopant element symbol |
| `dopant_fraction` | `float` | Fraction of host cations replaced |
| `n_vacancies` | `int` | Number of oxygen vacancies created |
| `parent_formula` | `str` | Original CeO2 formula |

---

## Dopant Analysis Functions

### analyze_dopant_distribution

```python
analyze_dopant_distribution(
    atoms: Atoms,
    dopant: str,
    reference_atoms: Atoms = None,
    z_threshold: float = 0.3,
    near_dopant_cutoff: float = 3.5
) -> Dict[str, Any]
```

Analyze dopant and vacancy distribution in a structure.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `atoms` | `Atoms` | - | Doped structure to analyze |
| `dopant` | `str` | - | Dopant element symbol |
| `reference_atoms` | `Atoms` | `None` | Original structure for vacancy counting |
| `z_threshold` | `float` | `0.3` | Surface region fraction |
| `near_dopant_cutoff` | `float` | `3.5` | Distance cutoff (Å) for near-dopant analysis |

**Returns:**

```python
{
    "dopant": str,
    "dopant_total": int,
    "dopant_surface": int,
    "dopant_surface_fraction": float,
    "vacancy_total": int,
    "vacancy_near_dopant": int,
    "vacancy_near_dopant_fraction": float,
    ...
}
```

**Example:**

```python
from nh3sofc.structure import analyze_dopant_distribution

stats = analyze_dopant_distribution(gdc, "Gd", reference_atoms=ceo2)
print(f"Gd in surface: {stats['dopant_surface_fraction']:.1%}")
print(f"Vacancies near Gd: {stats['vacancy_near_dopant_fraction']:.1%}")
```

### generate_dopant_series

```python
generate_dopant_series(
    structure: Union[Atoms, BaseStructure],
    dopant: str,
    dopant_fractions: List[float],
    n_configs: int = 3,
    host_cation: str = "Ce",
    random_seed: int = None
) -> List[Dict]
```

Generate structures with varying dopant concentrations.

**Example:**

```python
from nh3sofc.structure import generate_dopant_series

series = generate_dopant_series(
    ceo2_slab,
    dopant="Gd",
    dopant_fractions=[0.05, 0.10, 0.15, 0.20],
    n_configs=5,
)
```

### print_dopant_analysis

```python
print_dopant_analysis(stats: Dict[str, Any], title: str = "Dopant Distribution Analysis")
```

Print formatted summary of dopant distribution analysis.

---

## Dopant Constants

Available in `nh3sofc.core.constants`:

```python
from nh3sofc.core.constants import ACCEPTOR_DOPANTS, HOST_CATIONS, CHARGE_COMPENSATION

# Supported acceptor dopants for CeO2
ACCEPTOR_DOPANTS = {
    "Sm": {"charge": +3, "ionic_radius": 1.079, "name": "Samarium"},
    "Gd": {"charge": +3, "ionic_radius": 1.053, "name": "Gadolinium"},
    "Pr": {"charge": +3, "ionic_radius": 1.126, "name": "Praseodymium"},
    "Y":  {"charge": +3, "ionic_radius": 1.019, "name": "Yttrium"},
    "La": {"charge": +3, "ionic_radius": 1.160, "name": "Lanthanum"},
    "Nd": {"charge": +3, "ionic_radius": 1.109, "name": "Neodymium"},
}

# Host cations
HOST_CATIONS = {
    "Ce": {"charge": +4, "ionic_radius": 0.970},
    "Zr": {"charge": +4, "ionic_radius": 0.840},
}

# Charge compensation: 2 M³⁺ + 1 V_O^{2+} → neutral
CHARGE_COMPENSATION = {
    3: (2, 1),  # 2 trivalent dopants : 1 vacancy
    2: (1, 1),  # 1 divalent dopant : 1 vacancy
}
```
