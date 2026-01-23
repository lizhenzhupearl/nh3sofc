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
| `torch` | ≥2.0.0 | PyTorch for ML models |
| `mace-torch` | ≥0.3.0 | MACE ML force field |
| `spglib` | ≥2.0.0 | Space group analysis |
| `matplotlib` | ≥3.5.0 | Plotting and visualization |

### Optional Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `pymatgen` | ≥2022.0.0 | Advanced materials analysis |

## Installation Methods

### Method 1: Development Install (Recommended)

For active development or to get the latest features:

```bash
# Clone the repository
git clone https://github.com/lizhenzhupearl/nh3sofc.git
cd nh3sofc

# Create a virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install PyTorch first (choose one based on your system)
# CPU version:
pip install torch --index-url https://download.pytorch.org/whl/cpu

# Or GPU version (CUDA 11.8):
# pip install torch --index-url https://download.pytorch.org/whl/cu118
pip install mace-torch
pip install spglib
# Install nh3sofc with all dependencies
pip install -e .
```

### Method 2: Direct Install from GitHub

```bash
# Install PyTorch first (see above for CPU/GPU options)
pip install torch --index-url https://download.pytorch.org/whl/cpu

# Install nh3sofc
pip install git+https://github.com/lizhenzhupearl/nh3sofc.git
```

### Method 3: Install from PyPI (Future)

```bash
pip install nh3sofc  # Not yet available
```

### Installing Optional Dependencies

For advanced materials analysis with pymatgen:

```bash
pip install pymatgen
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
# Test core dependencies
import numpy as np
print(f"NumPy version: {np.__version__}")

import scipy
print(f"SciPy version: {scipy.__version__}")

import ase
print(f"ASE version: {ase.__version__}")

import torch
print(f"PyTorch version: {torch.__version__}")
print(f"CUDA available: {torch.cuda.is_available()}")

import mace
print(f"MACE-torch: OK")

import spglib
print(f"spglib version: {spglib.__version__}")

import matplotlib
print(f"Matplotlib version: {matplotlib.__version__}")

# Test nh3sofc package
import nh3sofc
print(f"nh3sofc version: {nh3sofc.__version__ if hasattr(nh3sofc, '__version__') else 'dev'}")

# Test structure module
from nh3sofc.structure import BulkStructure, SurfaceBuilder, AdsorbatePlacer
print("Structure module: OK")

# Test calculators
from nh3sofc.calculators.vasp import VASPInputGenerator, VASPOutputParser
print("VASP calculators: OK")

from nh3sofc.calculators.mace import MACECalculatorWrapper
print("MACE calculator: OK")

# Test workflows
from nh3sofc.workflows import RelaxationWorkflow, DecompositionWorkflow
print("Workflows: OK")

# Test analysis
from nh3sofc.analysis import AdsorptionEnergyCalculator, SurfaceComparator
print("Analysis: OK")

print("\n✓ Installation successful!")
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

#### 3. MACE or PyTorch installation fails

Ensure PyTorch is installed first with the correct version for your system:

```bash
# CPU version
pip install torch --index-url https://download.pytorch.org/whl/cpu

# Or GPU version (CUDA 11.8)
pip install torch --index-url https://download.pytorch.org/whl/cu118

# Then install MACE
pip install mace-torch
```

#### 4. spglib import error

```bash
pip install spglib
```

### Getting Help

- Open an issue on [GitHub](https://github.com/lizhenzhupearl/nh3sofc/issues)
- Check the [Tutorials](tutorials/index.md) for common workflows

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
