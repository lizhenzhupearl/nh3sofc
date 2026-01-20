# API Reference

Complete API documentation for the nh3sofc package.

## Package Structure

```
nh3sofc/
├── core/           # Base classes and constants
├── structure/      # Structure building
├── calculators/    # VASP and MACE interfaces
├── jobs/           # HPC job management
├── workflows/      # Automated workflows
├── analysis/       # Post-processing tools
└── database/       # Data management
```

## Modules

### Core

| Module | Description |
|--------|-------------|
| [`nh3sofc.core.constants`](core.md#constants) | Physical constants, default parameters |
| [`nh3sofc.core.base`](core.md#base-classes) | Base classes for structures and workflows |

### Structure Building

| Class | Description |
|-------|-------------|
| [`BulkStructure`](structure.md#bulkstructure) | Load and manipulate bulk crystals |
| [`SurfaceBuilder`](structure.md#surfacebuilder) | Create surface slabs |
| [`DefectBuilder`](structure.md#defectbuilder) | Create oxynitride defects |
| [`AdsorbatePlacer`](structure.md#adsorbateplacer) | Place adsorbate molecules |
| [`DecompositionBuilder`](structure.md#decompositionbuilder) | Generate decomposition intermediates |

### Calculators

| Class | Description |
|-------|-------------|
| [`VASPInputGenerator`](calculators.md#vaspinputgenerator) | Generate VASP input files |
| [`VASPOutputParser`](calculators.md#vaspoutputparser) | Parse VASP results |
| [`FrequencyCalculation`](calculators.md#frequencycalculation) | Vibrational analysis |
| [`MACECalculatorWrapper`](calculators.md#macecalculatorwrapper) | MACE ML force field |

### Workflows

| Class | Description |
|-------|-------------|
| [`RelaxationWorkflow`](workflows.md#relaxationworkflow) | Geometry optimization |
| [`DecompositionWorkflow`](workflows.md#decompositionworkflow) | NH3 decomposition pathway |
| [`NEBWorkflow`](workflows.md#nebworkflow) | Transition state search |
| [`FrequencyWorkflow`](workflows.md#frequencyworkflow) | Thermochemistry |
| [`ScreeningWorkflow`](workflows.md#screeningworkflow) | High-throughput screening |

### Analysis

| Class | Description |
|-------|-------------|
| [`AdsorptionEnergyCalculator`](analysis.md#adsorptionenergycalculator) | Adsorption energies |
| [`HarmonicThermo`](analysis.md#harmonicthermo) | Adsorbate thermodynamics |
| [`RDSAnalyzer`](analysis.md#rdsanalyzer) | Rate-determining step |
| [`SurfaceComparator`](analysis.md#surfacecomparator) | Compare catalyst surfaces |
| [`MicroKineticModel`](analysis.md#microkineticmodel) | Microkinetic modeling |

### Database

| Class | Description |
|-------|-------------|
| [`NamingConvention`](database.md#namingconvention) | Standardized naming |
| [`NH3SOFCDatabase`](database.md#nh3sofcdatabase) | ASE database wrapper |

## Quick Reference

### Common Imports

```python
# Structure
from nh3sofc.structure import (
    BulkStructure,
    SurfaceBuilder,
    DefectBuilder,
    AdsorbatePlacer,
    DecompositionBuilder,
)

# Calculators
from nh3sofc.calculators.vasp import (
    VASPInputGenerator,
    VASPOutputParser,
    FrequencyCalculation,
)
from nh3sofc.calculators.mace import (
    MACECalculatorWrapper,
    TrainingDataExtractor,
)

# Workflows
from nh3sofc.workflows import (
    RelaxationWorkflow,
    DecompositionWorkflow,
    NEBWorkflow,
    FrequencyWorkflow,
    ScreeningWorkflow,
)

# Analysis
from nh3sofc.analysis import (
    AdsorptionEnergyCalculator,
    HarmonicThermo,
    IdealGasThermo,
    RDSAnalyzer,
    SurfaceComparator,
    MicroKineticModel,
)

# Database
from nh3sofc.database import (
    NamingConvention,
    NH3SOFCDatabase,
)
```

### Physical Constants

```python
from nh3sofc.core.constants import (
    KB_EV,          # 8.617e-5 eV/K - Boltzmann constant
    H_EV_S,         # 4.136e-15 eV·s - Planck constant
    C_CMS,          # 2.998e10 cm/s - Speed of light
    EV_TO_J,        # 1.602e-19 J/eV - Energy conversion
    HUBBARD_U,      # Default Hubbard U values
    GAS_PHASE_ENERGIES,  # Reference gas energies
)
```

### Default VASP Parameters

```python
from nh3sofc.core.constants import DEFAULT_VASP_PARAMS

# Defaults:
# encut=520, ediff=1e-6, ediffg=-0.02
# ismear=0, sigma=0.05, kspacing=0.03
# ispin=2, lorbit=11, ncore=4
```
