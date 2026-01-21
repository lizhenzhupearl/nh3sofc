# Installation

## Requirements

### Python Version

- Python 3.9 or higher

### Core Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `ase` | ≥3.22.0 | Atomic Simulation Environment |
| `numpy` | ≥1.21.0 | Numerical computing |
| `scipy` | ≥1.7.0 | Scientific computing |

### Optional Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `mace-torch` | ≥0.3.0 | MACE ML force field |
| `spglib` | ≥2.0.0 | Space group analysis |
| `matplotlib` | ≥3.5.0 | Plotting |
| `pymatgen` | ≥2022.0.0 | Materials analysis (optional) |

## Installation Methods

### Method 1: Development Install (Recommended)

For active development or to get the latest features:

```bash
# Clone the repository
git clone https://github.com/lizhenzhupearl/nh3sofc.git
cd nh3sofc

# Create a virtual environment (optional but recommended)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install in development mode
pip install -e .
```

### Method 2: Direct Install from GitHub

```bash
pip install git+https://github.com/lizhenzhupearl/nh3sofc.git
```

### Method 3: Install from PyPI (Future)

```bash
pip install nh3sofc  # Not yet available
```

## Installing Optional Dependencies

### For MACE ML Force Fields

```bash
# Install PyTorch first (CPU version)
pip install torch --index-url https://download.pytorch.org/whl/cpu

# Or GPU version (CUDA 11.8)
pip install torch --index-url https://download.pytorch.org/whl/cu118

# Then install MACE
pip install mace-torch
```

### For Advanced Analysis

```bash
pip install matplotlib spglib pymatgen
```

### Complete Installation

```bash
pip install -e ".[all]"
```

## VASP Configuration

The package requires VASP pseudopotential files. Set the environment variable:

```bash
# Add to your ~/.bashrc or ~/.zshrc
export VASP_PP_PATH="/path/to/vasp/potentials"
```

The directory structure should be:

```
$VASP_PP_PATH/
├── potpaw_PBE/
│   ├── La/POTCAR
│   ├── V/POTCAR
│   ├── O/POTCAR
│   ├── N/POTCAR
│   ├── H/POTCAR
│   └── ...
└── potpaw_LDA/  # Optional
    └── ...
```

## PBS/HPC Configuration

For cluster job submission, configure your PBS settings:

```python
from nh3sofc.jobs import VASPJobScript

# Default settings can be customized
job = VASPJobScript(
    work_dir="./calc",
    nodes=2,
    ppn=24,
    walltime="24:00:00",
    queue="normal",
)
```

## Verifying Installation

Run the following to verify your installation:

```python
# Test basic imports
import nh3sofc
print(f"nh3sofc version: {nh3sofc.__version__ if hasattr(nh3sofc, '__version__') else 'dev'}")

# Test structure module
from nh3sofc.structure import BulkStructure, SurfaceBuilder, AdsorbatePlacer
print("Structure module: OK")

# Test calculators
from nh3sofc.calculators.vasp import VASPInputGenerator, VASPOutputParser
print("VASP calculators: OK")

# Test workflows
from nh3sofc.workflows import RelaxationWorkflow, DecompositionWorkflow
print("Workflows: OK")

# Test analysis
from nh3sofc.analysis import AdsorptionEnergyCalculator, SurfaceComparator
print("Analysis: OK")

# Test MACE (optional)
try:
    from nh3sofc.calculators.mace import MACECalculatorWrapper
    wrapper = MACECalculatorWrapper()
    print(f"MACE available: {wrapper._mace_available}")
except ImportError:
    print("MACE: Not installed (optional)")

print("\nInstallation successful!")
```

## Troubleshooting

### Common Issues

#### 1. Import Error: No module named 'ase'

```bash
pip install ase
```

#### 2. VASP POTCAR not found

Ensure `VASP_PP_PATH` is set correctly:

```bash
echo $VASP_PP_PATH
ls $VASP_PP_PATH/potpaw_PBE/
```

#### 3. MACE installation fails

Install PyTorch separately first:

```bash
pip install torch
pip install mace-torch
```

#### 4. spglib import error

```bash
pip install spglib
```

### Getting Help

- Open an issue on [GitHub](https://github.com/lizhenzhupearl/nh3sofc/issues)
- Check the [FAQ](faq.md)

## Updating

To update to the latest version:

```bash
cd nh3sofc
git pull
pip install -e .
```

## Uninstalling

```bash
pip uninstall nh3sofc
```
