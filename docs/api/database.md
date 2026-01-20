# Database Module

The `nh3sofc.database` module provides tools for organizing calculations and storing results.

## NamingConvention

Standardized naming for files and directories.

```python
from nh3sofc.database import NamingConvention
```

### Constructor

```python
NamingConvention(
    material: str,
    miller: tuple = None,
    termination: str = None,
    calc_type: str = None,
    suffix: str = None
)
```

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `material` | `str` | Material name (e.g., "LaVO3", "LaVON2") |
| `miller` | `tuple` | Miller indices (h, k, l) |
| `termination` | `str` | Surface termination (e.g., "LaO", "VO2") |
| `calc_type` | `str` | Calculation type |
| `suffix` | `str` | Additional suffix |

### Methods

#### get_directory_name

```python
get_directory_name() -> str
```

Generate standardized directory name.

**Format:** `{material}_{miller}_{termination}_{calc_type}_{suffix}`

#### get_filename

```python
get_filename(extension: str = "") -> str
```

Generate standardized filename.

#### from_atoms

```python
@classmethod
from_atoms(atoms: Atoms, calc_type: str = None) -> NamingConvention
```

Create from ASE Atoms object.

**Example:**

```python
from nh3sofc.database import NamingConvention

# Manual creation
naming = NamingConvention(
    material="LaVON2",
    miller=(0, 0, 1),
    termination="LaO",
    calc_type="relax",
    suffix="vac0.10"
)

dir_name = naming.get_directory_name()
# "LaVON2_001_LaO_relax_vac0.10"

filename = naming.get_filename(".traj")
# "LaVON2_001_LaO_relax_vac0.10.traj"

# From atoms
from ase.io import read
atoms = read("structure.xyz")
naming = NamingConvention.from_atoms(atoms, calc_type="static")
```

### Naming Conventions

| Component | Format | Examples |
|-----------|--------|----------|
| Material | Chemical formula | LaVO3, LaVON2, SrTiO3 |
| Miller | 3-digit | 001, 110, 111 |
| Termination | Element+O | LaO, VO2, TiO2 |
| Calc Type | lowercase | relax, static, neb, freq |
| Suffix | descriptive | vac0.10, NH3, config001 |

---

## NH3SOFCDatabase

ASE database wrapper with NH3-SOFC-specific functionality.

```python
from nh3sofc.database import NH3SOFCDatabase
```

### Constructor

```python
NH3SOFCDatabase(
    db_path: str = "nh3sofc.db",
    create: bool = True
)
```

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `db_path` | `str` | `"nh3sofc.db"` | Path to SQLite database |
| `create` | `bool` | `True` | Create if doesn't exist |

### Methods

#### add_structure

```python
add_structure(
    atoms: Atoms,
    material: str = None,
    miller: tuple = None,
    termination: str = None,
    calc_type: str = None,
    energy: float = None,
    adsorbate: str = None,
    **kwargs
) -> int
```

Add structure to database.

**Returns:** Row ID of added structure.

**Example:**

```python
db = NH3SOFCDatabase("my_calculations.db")

row_id = db.add_structure(
    atoms=optimized_atoms,
    material="LaVON2",
    miller=(0, 0, 1),
    termination="LaO",
    calc_type="relax",
    energy=-245.678,
    adsorbate="NH3",
    vacancy_concentration=0.10,
    nitrogen_fraction=0.67,
)
print(f"Added structure with ID: {row_id}")
```

#### get_structure

```python
get_structure(row_id: int) -> Atoms
```

Retrieve structure by ID.

#### query

```python
query(
    material: str = None,
    miller: tuple = None,
    calc_type: str = None,
    adsorbate: str = None,
    **kwargs
) -> List[Atoms]
```

Query structures matching criteria.

**Example:**

```python
# Find all NH3 adsorption calculations on LaVON2
structures = db.query(
    material="LaVON2",
    adsorbate="NH3",
    calc_type="relax"
)

print(f"Found {len(structures)} matching structures")
```

#### get_energies

```python
get_energies(
    material: str = None,
    miller: tuple = None,
    **kwargs
) -> Dict[str, float]
```

Get energies for matching structures.

#### get_lowest_energy

```python
get_lowest_energy(
    material: str = None,
    adsorbate: str = None,
    **kwargs
) -> Tuple[Atoms, float]
```

Get structure with lowest energy.

**Example:**

```python
# Find most stable NH3 adsorption configuration
atoms, energy = db.get_lowest_energy(
    material="LaVON2",
    adsorbate="NH3"
)
print(f"Lowest energy: {energy:.4f} eV")
```

#### add_trajectory

```python
add_trajectory(
    trajectory: List[Atoms],
    base_name: str,
    **common_kwargs
) -> List[int]
```

Add multiple structures from trajectory.

#### export_to_xyz

```python
export_to_xyz(
    filename: str,
    query_kwargs: dict = None
)
```

Export matching structures to XYZ file.

#### get_statistics

```python
get_statistics() -> dict
```

Get database statistics.

**Returns:**

```python
{
    "total_structures": int,
    "by_material": Dict[str, int],
    "by_calc_type": Dict[str, int],
    "by_adsorbate": Dict[str, int],
}
```

**Example:**

```python
stats = db.get_statistics()
print(f"Total structures: {stats['total_structures']}")
print("By material:")
for mat, count in stats['by_material'].items():
    print(f"  {mat}: {count}")
```

---

## Database Schema

The database stores structures with the following key-value pairs:

| Key | Type | Description |
|-----|------|-------------|
| `material` | `str` | Material formula |
| `miller` | `str` | Miller indices (e.g., "m001") |
| `termination` | `str` | Surface termination |
| `calc_type` | `str` | Calculation type |
| `adsorbate` | `str` | Adsorbate molecule |
| `vacancy_concentration` | `float` | Oxygen vacancy fraction |
| `nitrogen_fraction` | `float` | N/(N+O) ratio |
| `converged` | `bool` | Calculation converged |

Additional custom keys can be stored as needed.

---

## Usage Examples

### Building a Calculation Database

```python
from nh3sofc.database import NH3SOFCDatabase, NamingConvention
from nh3sofc.calculators.vasp import VASPOutputParser
from glob import glob

db = NH3SOFCDatabase("project.db")

# Add results from completed calculations
for calc_dir in glob("./calculations/*/"):
    parser = VASPOutputParser(calc_dir)

    if parser.check_convergence()["converged"]:
        results = parser.parse_outcar()
        atoms = parser.get_final_structure()

        # Extract metadata from directory name
        naming = NamingConvention.from_directory(calc_dir)

        db.add_structure(
            atoms=atoms,
            material=naming.material,
            miller=naming.miller,
            calc_type=naming.calc_type,
            energy=results["energy"],
        )

print(f"Database now contains {db.get_statistics()['total_structures']} structures")
```

### Comparing Surfaces

```python
from nh3sofc.database import NH3SOFCDatabase
from nh3sofc.analysis import SurfaceComparator

db = NH3SOFCDatabase("project.db")

# Get decomposition energies for different terminations
surfaces = {}

for term in ["LaO", "VO2"]:
    energies = db.get_energies(
        material="LaVON2",
        termination=term,
    )

    surfaces[term] = {
        step: energies.get(f"{term}_{step}", 0.0)
        for step in ["NH3", "NH2_H", "NH_2H", "N_3H"]
    }

comparator = SurfaceComparator(surfaces)
ranking = comparator.rank_by_energy_span()
```

### Exporting for Machine Learning

```python
db = NH3SOFCDatabase("project.db")

# Export all relaxed structures for ML training
db.export_to_xyz(
    "training_data.xyz",
    query_kwargs={"calc_type": "relax", "converged": True}
)
```

---

## Command-Line Interface

```bash
# List database contents
nh3sofc-db list project.db

# Query specific structures
nh3sofc-db query project.db --material LaVON2 --adsorbate NH3

# Export to XYZ
nh3sofc-db export project.db output.xyz --calc-type relax

# Show statistics
nh3sofc-db stats project.db
```
