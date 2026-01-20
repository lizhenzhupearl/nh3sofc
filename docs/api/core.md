# Core Module

The `nh3sofc.core` module provides fundamental constants, parameters, and base classes.

## Constants

Physical and chemical constants used throughout the package.

```python
from nh3sofc.core.constants import *
```

### Physical Constants

| Constant | Value | Unit | Description |
|----------|-------|------|-------------|
| `KB_EV` | 8.617333262e-5 | eV/K | Boltzmann constant |
| `H_EV_S` | 4.135667696e-15 | eV·s | Planck constant |
| `C_CMS` | 2.99792458e10 | cm/s | Speed of light |
| `EV_TO_J` | 1.602176634e-19 | J/eV | Energy conversion |
| `EV_TO_KJ_MOL` | 96.485 | kJ/mol/eV | Energy conversion |
| `AMU_TO_KG` | 1.66054e-27 | kg/amu | Mass conversion |

**Example:**

```python
from nh3sofc.core.constants import KB_EV, H_EV_S

# Calculate rate constant at 673 K
T = 673  # K
barrier = 1.0  # eV
prefactor = KB_EV * T / H_EV_S
k = prefactor * np.exp(-barrier / (KB_EV * T))
```

### Default VASP Parameters

```python
from nh3sofc.core.constants import DEFAULT_VASP_PARAMS

# Default settings
DEFAULT_VASP_PARAMS = {
    "encut": 520,
    "ediff": 1e-6,
    "ediffg": -0.02,
    "ismear": 0,
    "sigma": 0.05,
    "kspacing": 0.03,
    "ispin": 2,
    "lorbit": 11,
    "ncore": 4,
    "lreal": "Auto",
    "algo": "Normal",
    "prec": "Accurate",
}
```

### Hubbard U Values

```python
from nh3sofc.core.constants import HUBBARD_U

# Default Hubbard U values (eV)
HUBBARD_U = {
    "V": 3.25,
    "Ti": 3.0,
    "Mn": 3.9,
    "Fe": 4.0,
    "Co": 3.3,
    "Ni": 6.4,
    "Cr": 3.5,
    "Cu": 4.0,
}
```

### Gas Phase Reference Energies

```python
from nh3sofc.core.constants import GAS_PHASE_ENERGIES

# Reference energies (eV) - update with your calculated values
GAS_PHASE_ENERGIES = {
    "NH3": -19.54,
    "N2": -16.64,
    "H2": -6.77,
    "H2O": -14.22,
    "O2": -9.86,
}
```

### Element Properties

```python
from nh3sofc.core.constants import ELEMENT_MASSES, ATOMIC_RADII

# Masses in amu
ELEMENT_MASSES = {
    "H": 1.008,
    "N": 14.007,
    "O": 15.999,
    "La": 138.905,
    "V": 50.942,
    # ... more elements
}

# Covalent radii in Å
ATOMIC_RADII = {
    "H": 0.31,
    "N": 0.71,
    "O": 0.66,
    # ... more elements
}
```

---

## Base Classes

Abstract base classes for structures and workflows.

### BaseStructure

Base class for structure objects.

```python
from nh3sofc.core.base import BaseStructure
```

```python
class BaseStructure:
    """Base class for structure manipulation."""

    def __init__(self, atoms: Atoms):
        self.atoms = atoms

    def copy(self) -> 'BaseStructure':
        """Return a deep copy."""
        ...

    def write(self, filename: str, format: str = None):
        """Write structure to file."""
        ...

    def get_composition(self) -> Dict[str, int]:
        """Get element counts."""
        ...

    def get_center_of_mass(self) -> np.ndarray:
        """Get center of mass."""
        ...
```

### BaseWorkflow

Base class for workflow objects.

```python
from nh3sofc.core.base import BaseWorkflow
```

```python
class BaseWorkflow:
    """Base class for calculation workflows."""

    def __init__(self, work_dir: str = "."):
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def setup(self) -> dict:
        """Set up calculation. Override in subclass."""
        raise NotImplementedError

    def run(self) -> Any:
        """Run calculation. Override in subclass."""
        raise NotImplementedError

    def parse_results(self) -> dict:
        """Parse results. Override in subclass."""
        raise NotImplementedError

    def get_status(self) -> str:
        """Check calculation status."""
        ...
```

### BaseAnalyzer

Base class for analysis objects.

```python
from nh3sofc.core.base import BaseAnalyzer
```

```python
class BaseAnalyzer:
    """Base class for result analysis."""

    def __init__(self):
        self.results = {}

    def analyze(self, data: Any) -> dict:
        """Perform analysis. Override in subclass."""
        raise NotImplementedError

    def plot(self, filename: str = None):
        """Generate plot. Override in subclass."""
        raise NotImplementedError

    def to_dataframe(self) -> pd.DataFrame:
        """Convert results to DataFrame."""
        ...
```

---

## Utility Functions

### Energy Conversions

```python
from nh3sofc.core.utils import ev_to_kj_mol, kj_mol_to_ev, ev_to_kcal_mol

energy_ev = 1.0
energy_kj = ev_to_kj_mol(energy_ev)   # 96.485 kJ/mol
energy_kcal = ev_to_kcal_mol(energy_ev)  # 23.06 kcal/mol
```

### Temperature/Frequency Conversions

```python
from nh3sofc.core.utils import cm_to_ev, ev_to_cm, kelvin_to_ev

freq_cm = 3400  # cm^-1
freq_ev = cm_to_ev(freq_cm)  # ~0.42 eV

T = 673  # K
kT = kelvin_to_ev(T)  # ~0.058 eV
```

### Structure Utilities

```python
from nh3sofc.core.utils import (
    get_surface_atoms,
    get_subsurface_atoms,
    get_adsorbate_indices,
    calculate_rmsd,
)

# Find surface atoms (top layer)
surface_indices = get_surface_atoms(slab, tolerance=0.5)

# Get adsorbate atom indices
ads_indices = get_adsorbate_indices(slab_with_adsorbate, slab)

# Calculate RMSD between structures
rmsd = calculate_rmsd(atoms1, atoms2)
```

---

## Configuration

### Package Configuration

```python
from nh3sofc.core.config import Config

# Get configuration
config = Config()

# Access settings
vasp_pp_path = config.get("vasp_pp_path")
default_calculator = config.get("default_calculator", "vasp")

# Set configuration
config.set("vasp_pp_path", "/path/to/potentials")
config.save()
```

### Configuration File

Create `~/.nh3sofc/config.yaml`:

```yaml
# VASP settings
vasp_pp_path: /path/to/vasp/potentials
vasp_command: mpirun vasp_std

# MACE settings
mace_model_path: /path/to/models
default_foundation_model: medium

# Job submission
job_system: pbs
default_nodes: 1
default_ppn: 24
default_walltime: "24:00:00"

# Calculation defaults
default_encut: 520
default_kspacing: 0.03
use_vdw: true
vdw_type: D3BJ
```

### Environment Variables

| Variable | Description |
|----------|-------------|
| `VASP_PP_PATH` | Path to VASP pseudopotentials |
| `NH3SOFC_CONFIG` | Path to config file |
| `NH3SOFC_DB` | Default database path |

---

## Logging

```python
from nh3sofc.core.logging import get_logger

logger = get_logger(__name__)

logger.info("Starting calculation")
logger.debug("Debug information")
logger.warning("Warning message")
logger.error("Error occurred")
```

Configure logging level:

```python
from nh3sofc.core.logging import set_log_level

set_log_level("DEBUG")  # or "INFO", "WARNING", "ERROR"
```
