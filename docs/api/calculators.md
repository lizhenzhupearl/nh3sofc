# Calculators Module

The `nh3sofc.calculators` module provides interfaces for VASP DFT and MACE ML force field calculations.

## VASP Calculators

### VASPInputGenerator

Generate VASP input files (INCAR, KPOINTS, POSCAR, POTCAR).

```python
from nh3sofc.calculators.vasp import VASPInputGenerator
```

#### Constructor

```python
VASPInputGenerator(
    atoms: Atoms,
    calc_type: str = "relax",
    work_dir: str = ".",
    potcar_dir: str = None
)
```

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `atoms` | `Atoms` | - | ASE Atoms object |
| `calc_type` | `str` | `"relax"` | Calculation type |
| `work_dir` | `str` | `"."` | Working directory |
| `potcar_dir` | `str` | `None` | POTCAR directory (uses VASP_PP_PATH env) |

**Calculation Types:**

| Type | Description | Key Settings |
|------|-------------|--------------|
| `"relax"` | Geometry optimization | IBRION=2, NSW=300 |
| `"static"` | Single-point energy | NSW=0 |
| `"frequency"` | Vibrational analysis | IBRION=5, NFREE=2 |
| `"neb"` | NEB transition state | IMAGES, LCLIMB |
| `"md"` | Molecular dynamics | IBRION=0, SMASS |

#### Methods

##### generate_all

```python
generate_all(
    encut: float = 520,
    ediff: float = 1e-6,
    ediffg: float = -0.02,
    kspacing: float = 0.03,
    hubbard_u: dict = None,
    vdw: str = None,
    **kwargs
) -> dict
```

Generate all VASP input files.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `encut` | `float` | `520` | Plane-wave cutoff (eV) |
| `ediff` | `float` | `1e-6` | Electronic convergence |
| `ediffg` | `float` | `-0.02` | Ionic convergence (eV/Å) |
| `kspacing` | `float` | `0.03` | K-point spacing (1/Å) |
| `hubbard_u` | `dict` | `None` | Hubbard U values, e.g., `{"V": 3.25}` |
| `vdw` | `str` | `None` | vdW correction: "D3", "D3BJ", "dDsC" |

**Returns:** Dictionary with paths to generated files.

**Example:**

```python
from ase.io import read
from nh3sofc.calculators.vasp import VASPInputGenerator

atoms = read("structure.xyz")

vasp = VASPInputGenerator(atoms, calc_type="relax", work_dir="./calc")
files = vasp.generate_all(
    encut=520,
    kspacing=0.03,
    hubbard_u={"V": 3.25, "Ti": 3.0},
    vdw="D3BJ",
    ispin=2,
    lorbit=11,
)

print(f"Generated: {list(files.keys())}")
```

##### generate_incar

```python
generate_incar(**kwargs) -> str
```

Generate INCAR file content.

##### generate_kpoints

```python
generate_kpoints(kspacing: float = 0.03) -> str
```

Generate KPOINTS file content.

##### generate_potcar

```python
generate_potcar() -> str
```

Generate POTCAR by concatenating pseudopotentials.

---

### VASPOutputParser

Parse VASP output files.

```python
from nh3sofc.calculators.vasp import VASPOutputParser
```

#### Constructor

```python
VASPOutputParser(work_dir: str = ".")
```

#### Methods

##### parse_outcar

```python
parse_outcar() -> dict
```

Parse OUTCAR for energy, forces, stress.

**Returns:**

```python
{
    "energy": float,           # Total energy (eV)
    "forces": np.ndarray,      # Forces (eV/Å)
    "stress": np.ndarray,      # Stress tensor (kBar)
    "converged": bool,         # Electronic convergence
    "ionic_steps": int,        # Number of ionic steps
    "magmom": float,           # Total magnetic moment
}
```

##### parse_vasprun

```python
parse_vasprun() -> dict
```

Parse vasprun.xml for detailed results.

##### get_trajectory

```python
get_trajectory() -> List[Atoms]
```

Extract trajectory from OUTCAR.

##### check_convergence

```python
check_convergence() -> dict
```

Check electronic and ionic convergence.

**Example:**

```python
parser = VASPOutputParser("./calc")

results = parser.parse_outcar()
print(f"Energy: {results['energy']:.4f} eV")
print(f"Converged: {results['converged']}")

# Get relaxation trajectory
traj = parser.get_trajectory()
print(f"Ionic steps: {len(traj)}")
```

---

### FrequencyCalculation

Set up and parse vibrational frequency calculations.

```python
from nh3sofc.calculators.vasp import FrequencyCalculation
```

#### Constructor

```python
FrequencyCalculation(
    atoms: Atoms,
    work_dir: str = ".",
    selective_dynamics: List[int] = None
)
```

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `atoms` | `Atoms` | Optimized structure |
| `work_dir` | `str` | Working directory |
| `selective_dynamics` | `List[int]` | Indices of atoms to vibrate |

#### Methods

##### setup

```python
setup(
    nfree: int = 2,
    potim: float = 0.015,
    **vasp_kwargs
) -> dict
```

Set up frequency calculation.

##### parse_frequencies

```python
parse_frequencies() -> dict
```

Parse frequencies from OUTCAR.

**Returns:**

```python
{
    "frequencies": List[float],      # All frequencies (cm^-1)
    "real_frequencies": List[float], # Real frequencies
    "imaginary": List[float],        # Imaginary frequencies
    "zpe": float,                    # Zero-point energy (eV)
}
```

**Example:**

```python
# After geometry optimization
freq = FrequencyCalculation(
    optimized_atoms,
    work_dir="./freq",
    selective_dynamics=adsorbate_indices  # Only adsorbate atoms
)

freq.setup(nfree=2, encut=520)

# After VASP calculation completes...
results = freq.parse_frequencies()
print(f"ZPE: {results['zpe']:.4f} eV")
print(f"Frequencies: {results['real_frequencies']}")

if results['imaginary']:
    print(f"WARNING: Imaginary frequencies found: {results['imaginary']}")
```

---

## MACE Calculators

### MACECalculatorWrapper

Wrapper for MACE ML force field calculator.

```python
from nh3sofc.calculators.mace import MACECalculatorWrapper
```

#### Constructor

```python
MACECalculatorWrapper(
    model_path: str = None,
    foundation_model: str = "medium",
    device: str = "auto",
    default_dtype: str = "float64"
)
```

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `model_path` | `str` | `None` | Path to custom MACE model |
| `foundation_model` | `str` | `"medium"` | Foundation model: "small", "medium", "large" |
| `device` | `str` | `"auto"` | Device: "cpu", "cuda", "auto" |
| `default_dtype` | `str` | `"float64"` | Numerical precision |

#### Methods

##### get_calculator

```python
get_calculator() -> Calculator
```

Get ASE calculator instance.

##### calculate

```python
calculate(atoms: Atoms) -> dict
```

Calculate energy, forces, stress.

**Example:**

```python
from nh3sofc.calculators.mace import MACECalculatorWrapper

# Using foundation model
mace = MACECalculatorWrapper(foundation_model="medium")

# Or custom trained model
mace = MACECalculatorWrapper(model_path="./my_model.model")

# Attach to atoms
calc = mace.get_calculator()
atoms.calc = calc

energy = atoms.get_potential_energy()
forces = atoms.get_forces()
```

---

### MACEEnsemble

Ensemble of MACE models for uncertainty estimation.

```python
from nh3sofc.calculators.mace import MACEEnsemble
```

#### Constructor

```python
MACEEnsemble(
    model_paths: List[str],
    device: str = "auto"
)
```

#### Methods

##### calculate_with_uncertainty

```python
calculate_with_uncertainty(atoms: Atoms) -> dict
```

Calculate with uncertainty estimates.

**Returns:**

```python
{
    "energy": float,
    "energy_std": float,
    "forces": np.ndarray,
    "forces_std": np.ndarray,
    "max_force_std": float,
}
```

**Example:**

```python
ensemble = MACEEnsemble([
    "./model_1.model",
    "./model_2.model",
    "./model_3.model",
])

results = ensemble.calculate_with_uncertainty(atoms)
print(f"Energy: {results['energy']:.4f} ± {results['energy_std']:.4f} eV")

if results['max_force_std'] > 0.1:
    print("High uncertainty - consider adding to training data")
```

---

### TrainingDataExtractor

Extract training data from VASP calculations.

```python
from nh3sofc.calculators.mace import TrainingDataExtractor
```

#### Constructor

```python
TrainingDataExtractor(vasp_dirs: List[str] = None)
```

#### Methods

##### extract_from_directory

```python
extract_from_directory(
    vasp_dir: str,
    include_trajectory: bool = True
) -> List[Atoms]
```

Extract structures with energies and forces.

##### extract_all

```python
extract_all() -> List[Atoms]
```

Extract from all directories.

##### write_xyz

```python
write_xyz(
    filename: str,
    atoms_list: List[Atoms] = None
)
```

Write training data to extended XYZ format.

##### filter_by_energy

```python
filter_by_energy(
    atoms_list: List[Atoms],
    max_energy_per_atom: float = 0.0
) -> List[Atoms]
```

Filter high-energy configurations.

**Example:**

```python
extractor = TrainingDataExtractor()

# Extract from multiple calculations
for calc_dir in glob("./calculations/*/"):
    extractor.extract_from_directory(calc_dir)

# Get all training data
training_data = extractor.extract_all()
print(f"Extracted {len(training_data)} configurations")

# Filter and save
filtered = extractor.filter_by_energy(training_data, max_energy_per_atom=-2.0)
extractor.write_xyz("training_data.xyz", filtered)
```

---

### MACETrainingConfig

Generate MACE training configuration.

```python
from nh3sofc.calculators.mace import MACETrainingConfig
```

#### Constructor

```python
MACETrainingConfig(
    train_file: str,
    valid_file: str = None,
    model_name: str = "MACE_model"
)
```

#### Methods

##### generate_config

```python
generate_config(
    r_max: float = 5.0,
    num_radial_basis: int = 8,
    num_cutoff_basis: int = 5,
    max_L: int = 2,
    num_channels: int = 128,
    max_epochs: int = 1000,
    batch_size: int = 10,
    **kwargs
) -> dict
```

Generate training configuration dictionary.

##### write_config

```python
write_config(filename: str)
```

Write configuration to YAML file.

##### generate_training_script

```python
generate_training_script(filename: str = "train.sh")
```

Generate shell script for training.

**Example:**

```python
config = MACETrainingConfig(
    train_file="train.xyz",
    valid_file="valid.xyz",
    model_name="LaVON_MACE"
)

config.generate_config(
    r_max=6.0,
    max_epochs=500,
    batch_size=5,
)

config.write_config("mace_config.yaml")
config.generate_training_script("train_mace.sh")
```

---

## Preset Parameters

### VASP Presets

```python
from nh3sofc.calculators.vasp import get_preset

# Available presets
preset = get_preset("relax")      # Geometry optimization
preset = get_preset("static")     # Single-point
preset = get_preset("frequency")  # Vibrational analysis
preset = get_preset("neb")        # NEB calculation
preset = get_preset("md")         # Molecular dynamics
```

### Default Hubbard U Values

```python
from nh3sofc.core.constants import HUBBARD_U

# Default values (eV)
HUBBARD_U = {
    "V": 3.25,
    "Ti": 3.0,
    "Mn": 3.9,
    "Fe": 4.0,
    "Co": 3.3,
    "Ni": 6.4,
}
```
