# Workflows Module

The `nh3sofc.workflows` module provides high-level workflow classes for common computational tasks.

## RelaxationWorkflow

Geometry optimization workflow.

```python
from nh3sofc.workflows import RelaxationWorkflow
```

### Constructor

```python
RelaxationWorkflow(
    atoms: Atoms,
    work_dir: str = ".",
    calculator: str = "vasp",
    **calc_kwargs
)
```

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `atoms` | `Atoms` | - | Initial structure |
| `work_dir` | `str` | `"."` | Working directory |
| `calculator` | `str` | `"vasp"` | Calculator: "vasp" or "mace" |
| `**calc_kwargs` | - | - | Calculator parameters |

### Methods

#### setup

```python
setup() -> dict
```

Set up calculation files.

#### run

```python
run(submit: bool = False) -> Atoms
```

Run optimization (MACE) or generate files (VASP).

#### parse_results

```python
parse_results() -> dict
```

Parse optimization results.

**Example:**

```python
from ase.io import read
from nh3sofc.workflows import RelaxationWorkflow

atoms = read("initial.xyz")

# VASP relaxation
workflow = RelaxationWorkflow(
    atoms,
    work_dir="./relax",
    calculator="vasp",
    encut=520,
    kspacing=0.03,
    hubbard_u={"V": 3.25},
)

workflow.setup()

# After VASP completes...
results = workflow.parse_results()
relaxed = results["atoms"]
print(f"Final energy: {results['energy']:.4f} eV")
```

---

## DecompositionWorkflow

NH3 decomposition pathway workflow.

```python
from nh3sofc.workflows import DecompositionWorkflow
```

### Constructor

```python
DecompositionWorkflow(
    nh3_on_slab: Atoms,
    work_dir: str = "./decomposition",
    n_configs_per_step: int = 5,
    calculator: str = "vasp",
    **calc_kwargs
)
```

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `nh3_on_slab` | `Atoms` | - | Optimized NH3 on surface |
| `work_dir` | `str` | `"./decomposition"` | Working directory |
| `n_configs_per_step` | `int` | `5` | Configurations per step |
| `calculator` | `str` | `"vasp"` | Calculator type |

### Methods

#### setup

```python
setup() -> Dict[str, List[str]]
```

Generate all intermediate configurations and calculation files.

**Returns:** Dictionary mapping step names to list of calculation directories.

#### get_status

```python
get_status() -> Dict[str, dict]
```

Check calculation status for all steps.

#### parse_results

```python
parse_results() -> Dict[str, List[dict]]
```

Parse all calculation results.

#### get_lowest_energies

```python
get_lowest_energies() -> Dict[str, float]
```

Get lowest energy for each step.

#### get_energy_profile

```python
get_energy_profile(reference: str = "NH3") -> Dict[str, float]
```

Get relative energy profile.

#### get_reaction_energies

```python
get_reaction_energies() -> Dict[str, float]
```

Get step-by-step reaction energies.

#### print_summary

```python
print_summary()
```

Print formatted results summary.

**Example:**

```python
from ase.io import read
from nh3sofc.workflows import DecompositionWorkflow

nh3_on_slab = read("relaxed_nh3_on_surface.traj")

workflow = DecompositionWorkflow(
    nh3_on_slab=nh3_on_slab,
    work_dir="./decomposition_study",
    n_configs_per_step=5,
    encut=520,
    hubbard_u={"V": 3.25},
    vdw="D3BJ",
)

# Setup all calculations
paths = workflow.setup()
for step, dirs in paths.items():
    print(f"{step}: {len(dirs)} configurations")

# After calculations complete...
profile = workflow.get_energy_profile()
workflow.print_summary()
```

### Directory Structure

```
decomposition_study/
├── initial_configs/
│   ├── NH3.xyz
│   ├── NH2_H_000.xyz
│   └── ...
├── NH3/
│   └── config_000/
├── NH2_H/
│   ├── config_000/
│   ├── config_001/
│   └── ...
├── NH_2H/
│   └── ...
└── N_3H/
    └── ...
```

---

## NEBWorkflow

Nudged Elastic Band workflow for transition states.

```python
from nh3sofc.workflows import NEBWorkflow
```

### Constructor

```python
NEBWorkflow(
    initial: Atoms,
    final: Atoms,
    work_dir: str = "./neb",
    n_images: int = 5,
    calculator: str = "vasp",
    **calc_kwargs
)
```

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `initial` | `Atoms` | - | Initial state structure |
| `final` | `Atoms` | - | Final state structure |
| `work_dir` | `str` | `"./neb"` | Working directory |
| `n_images` | `int` | `5` | Number of NEB images |
| `calculator` | `str` | `"vasp"` | Calculator type |

### Methods

#### setup

```python
setup(
    interpolation: str = "idpp",
    climb: bool = True
) -> dict
```

Set up NEB calculation.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `interpolation` | `str` | `"idpp"` | Interpolation: "linear" or "idpp" |
| `climb` | `bool` | `True` | Use climbing image NEB |

#### run_mace

```python
run_mace(
    fmax: float = 0.05,
    steps: int = 500
) -> List[Atoms]
```

Run NEB with MACE calculator.

#### parse_results

```python
parse_results() -> dict
```

Parse NEB results.

**Returns:**

```python
{
    "images": List[Atoms],
    "energies": List[float],
    "barrier_forward": float,
    "barrier_reverse": float,
    "ts_image": int,
    "ts_structure": Atoms,
}
```

#### plot_energy_profile

```python
plot_energy_profile(filename: str = None)
```

Plot NEB energy profile.

**Example:**

```python
from ase.io import read
from nh3sofc.workflows import NEBWorkflow

initial = read("nh3_on_surface.traj")
final = read("nh2_h_on_surface.traj")

neb = NEBWorkflow(
    initial=initial,
    final=final,
    work_dir="./neb_nh3_dissociation",
    n_images=7,
    encut=520,
)

# For VASP
neb.setup(interpolation="idpp", climb=True)

# Or run with MACE
neb_images = neb.run_mace(fmax=0.03)

results = neb.parse_results()
print(f"Forward barrier: {results['barrier_forward']:.3f} eV")
print(f"Reverse barrier: {results['barrier_reverse']:.3f} eV")

neb.plot_energy_profile("neb_profile.png")
```

---

## FrequencyWorkflow

Vibrational frequency and thermochemistry workflow.

```python
from nh3sofc.workflows import FrequencyWorkflow
```

### Constructor

```python
FrequencyWorkflow(
    atoms: Atoms,
    work_dir: str = "./frequency",
    adsorbate_indices: List[int] = None,
    calculator: str = "vasp",
    **calc_kwargs
)
```

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `atoms` | `Atoms` | Optimized structure |
| `work_dir` | `str` | Working directory |
| `adsorbate_indices` | `List[int]` | Indices of adsorbate atoms |
| `calculator` | `str` | Calculator type |

### Methods

#### setup

```python
setup(nfree: int = 2) -> dict
```

Set up frequency calculation.

#### parse_results

```python
parse_results() -> dict
```

Parse frequencies and calculate thermochemistry.

**Returns:**

```python
{
    "frequencies": List[float],
    "real_frequencies": List[float],
    "imaginary": List[float],
    "zpe": float,
    "thermo": HarmonicThermo,
}
```

#### get_thermochemistry

```python
get_thermochemistry(
    temperature: float = 673.0,
    pressure: float = 1.0
) -> dict
```

Calculate thermodynamic properties.

**Returns:**

```python
{
    "zpe": float,
    "U_vib": float,
    "S_vib": float,
    "G": float,
    "H": float,
}
```

**Example:**

```python
from nh3sofc.workflows import FrequencyWorkflow

workflow = FrequencyWorkflow(
    optimized_atoms,
    work_dir="./freq",
    adsorbate_indices=[48, 49, 50, 51],  # NH3 atom indices
    encut=520,
)

workflow.setup()

# After VASP completes...
results = workflow.parse_results()
print(f"Frequencies (cm^-1): {results['real_frequencies']}")
print(f"ZPE: {results['zpe']:.4f} eV")

# Thermochemistry at 673 K
thermo = workflow.get_thermochemistry(temperature=673)
print(f"G(673K): {thermo['G']:.4f} eV")
```

---

## ScreeningWorkflow

High-throughput parameter screening.

```python
from nh3sofc.workflows import ScreeningWorkflow
```

### Constructor

```python
ScreeningWorkflow(
    base_structure: Atoms,
    parameter_space: Dict[str, List],
    work_dir: str = "./screening",
    calculator: str = "vasp",
    **calc_kwargs
)
```

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `base_structure` | `Atoms` | Base structure to modify |
| `parameter_space` | `dict` | Parameters to scan |
| `work_dir` | `str` | Working directory |

### Methods

#### setup

```python
setup() -> Dict[str, str]
```

Generate all calculation directories.

#### get_status

```python
get_status() -> Dict[str, str]
```

Check status of all calculations.

#### parse_results

```python
parse_results() -> pd.DataFrame
```

Parse all results into DataFrame.

#### find_optimal

```python
find_optimal(metric: str = "energy") -> dict
```

Find optimal parameters.

**Example:**

```python
from nh3sofc.workflows import ScreeningWorkflow

workflow = ScreeningWorkflow(
    base_structure=surface,
    parameter_space={
        "vacancy_concentration": [0.0, 0.05, 0.10, 0.15, 0.20],
        "nitrogen_fraction": [0.5, 0.67, 0.75],
    },
    work_dir="./oxynitride_screening",
    encut=520,
)

# Setup all calculations
calc_dirs = workflow.setup()
print(f"Created {len(calc_dirs)} calculations")

# After completion...
results = workflow.parse_results()
optimal = workflow.find_optimal(metric="adsorption_energy")
print(f"Optimal parameters: {optimal}")
```

---

## CompositionScreening

Specialized workflow for composition screening.

```python
from nh3sofc.workflows import CompositionScreening
```

### Constructor

```python
CompositionScreening(
    base_structure: Atoms,
    elements_to_vary: Dict[str, List[str]],
    work_dir: str = "./composition_screening",
    **calc_kwargs
)
```

**Example:**

```python
workflow = CompositionScreening(
    base_structure=surface,
    elements_to_vary={
        "La": ["La", "Sr", "Ba"],
        "V": ["V", "Ti", "Nb"],
    },
    work_dir="./dopant_screening",
)

workflow.setup()
```

---

## AdsorbateScreening

Screen multiple adsorbates and sites.

```python
from nh3sofc.workflows import AdsorbateScreening
```

### Constructor

```python
AdsorbateScreening(
    surface: Atoms,
    adsorbates: List[str],
    sites: List[str] = None,
    work_dir: str = "./adsorbate_screening",
    **calc_kwargs
)
```

**Example:**

```python
workflow = AdsorbateScreening(
    surface=surface,
    adsorbates=["NH3", "NH2", "NH", "N", "H"],
    sites=["ontop_La", "ontop_V", "bridge", "hollow"],
    work_dir="./adsorbate_study",
)

paths = workflow.setup()
# Creates calculations for all adsorbate-site combinations
```

---

## Workflow Utilities

### submit_all

```python
from nh3sofc.workflows import submit_all

submit_all(
    work_dir: str,
    job_system: str = "pbs",
    dry_run: bool = False
)
```

Submit all jobs in a workflow directory.

### check_all_complete

```python
from nh3sofc.workflows import check_all_complete

complete = check_all_complete(work_dir)
if complete:
    print("All calculations finished")
```

### collect_results

```python
from nh3sofc.workflows import collect_results

results = collect_results(
    work_dir,
    output_format="dataframe"  # or "dict" or "json"
)
```
