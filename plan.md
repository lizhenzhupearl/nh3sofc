# NH3-SOFC Simulation Package Implementation Plan

## Overview
Python package (`nh3sofc`) for automating DFT (VASP) and ML force field (MACE) calculations for NH3 solid oxide fuel cell research, with focus on perovskite oxynitrides (LaVO3 → LaVON₂₋ₓ).

## Complete Workflow

### STEP 1: Build Slabs
1.1 Build surface from bulk CIF
1.2 Place molecules (6 methods: manual, random, grid, site-specific, collision-aware, CatKit)
1.3 Beyond 0K: Gibbs free energy considerations

### STEP 2: Perform Calculations
- Co-adsorption approach: NH3* → NH2*+H* → NH*+2H* → N*+3H*
- Multiple configurations per step (different H positions)
- H2 formation and desorption

### STEP 3: Understand Results
- Thermodynamic criteria (reaction energies, H2 desorption)
- Kinetic criteria (RDS, activation barriers)
- Surface comparison framework (volcano plots, scoring)

### STEP 4: Rate-Limiting Step Analysis
- Thermodynamic approximation (highest ΔE)
- Energy span model
- NEB transition states
- Microkinetic modeling

## Package Structure

```
nh3sofc/
├── __init__.py
├── core/
│   ├── constants.py          # Physical constants, default parameters
│   └── base.py               # Base classes
├── structure/
│   ├── bulk.py               # Bulk structure loading/manipulation
│   ├── surface.py            # Surface generation (Miller indices, terminations)
│   ├── defects.py            # O→N substitution, vacancy creation with concentration x
│   ├── adsorbates.py         # 6 placement methods
│   └── decomposition.py      # NH3 decomposition configuration generator
├── calculators/
│   ├── vasp/
│   │   ├── inputs.py         # INCAR, KPOINTS, POTCAR, POSCAR generation
│   │   ├── outputs.py        # Parse OUTCAR, vasprun.xml, CONTCAR
│   │   ├── presets.py        # Preset settings for relax/NEB/MD/frequency
│   │   └── frequency.py      # Vibrational analysis setup (IBRION=5,6,7)
│   └── mace/
│       ├── interface.py      # MACE calculator wrapper
│       ├── training.py       # Training data generation
│       └── active_learning.py # VASP→MACE→uncertainty→VASP loop
├── jobs/
│   ├── pbs.py                # PBS script generation
│   └── manager.py            # Job submission, monitoring, restart
├── workflows/
│   ├── relaxation.py         # Geometry optimization workflow
│   ├── adsorption.py         # Adsorption energy calculations
│   ├── decomposition.py      # NH3 decomposition pathway (co-adsorption)
│   ├── neb.py                # NEB reaction path workflow
│   ├── frequency.py          # Vibrational/ZPE calculation workflow
│   ├── md.py                 # MD simulation workflow
│   └── screening.py          # High-throughput parameter screening
├── analysis/
│   ├── energetics.py         # Adsorption/surface energies
│   ├── thermochemistry.py    # Gibbs free energy (ZPE, thermal corrections, entropy)
│   ├── barriers.py           # NEB barrier extraction, BEP relations
│   ├── trajectory.py         # MD analysis (MSD, RDF, diffusion)
│   ├── rds.py                # Rate-determining step identification
│   ├── surface_comparison.py # Surface ranking, volcano plots
│   └── microkinetic.py       # Microkinetic modeling (steady-state, TOF)
├── database/
│   ├── naming.py             # Directory naming conventions
│   └── asedb.py              # ASE database integration
└── cli/
    └── main.py               # Command-line interface
```

## Key Modules

### Structure Building

#### `structure/adsorbates.py` - 6 Adsorbate Placement Methods
1. `add_simple()` - Manual (x,y) + height using `add_adsorbate()`
2. `add_random()` - Random placement with Euler angle rotation
3. `add_grid()` - Grid-based systematic placement
4. `add_on_site()` - Site-specific (above La, V, O atoms)
5. `add_with_collision()` - Random with minimum distance check
6. `add_catkit()` - CatKit integration (top/bridge/hollow sites)

#### `structure/defects.py` - Oxynitride Generation
- `create_oxynitride(nitrogen_fraction, vacancy_concentration)` - O→N substitution + vacancies
- Parameter `x` in LaVON₂₋ₓ controls vacancy concentration

#### `structure/decomposition.py` - NH3 Decomposition Configurations
- `create_NH2_H_configs()` - Remove H, place on various sites
- `create_NH_2H_configs()` - Multiple H-H distance configurations
- `create_N_3H_configs()` - Clustered vs spread H arrangements
- `filter_unique_configs()` - RMSD-based duplicate removal

### Analysis

#### `analysis/thermochemistry.py` - Gibbs Free Energy
- ΔG = ΔE₀ + ΔZPE + ΔH(T) - TΔS(T)
- `calculate_ZPE()` - Zero-point energy from frequencies
- `calculate_thermal_correction()` - ΔH(T) enthalpy
- `calculate_entropy()` - Vibrational, rotational, translational S(T)

#### `analysis/rds.py` - Rate-Determining Step
- `thermodynamic_rds()` - Highest ΔE step (approximation)
- `energy_span_model()` - Campbell's energy span approach
- `brep_barrier()` - Brønsted-Evans-Polanyi estimation

#### `analysis/surface_comparison.py` - Surface Ranking
- `SurfaceComparator` class with metrics
- Energy span, RDS barrier, N binding strength
- `plot_energy_profiles()` - Visual comparison
- `volcano_plot()` - Activity vs descriptor

## Key API Examples

### 1. Create Oxynitride Surface
```python
from nh3sofc.structure import BulkStructure, SurfaceBuilder, DefectBuilder

bulk = BulkStructure.from_cif("LaVO3.cif")
surface = SurfaceBuilder(bulk).create_surface((0,0,1), layers=6, vacuum=15)
oxynitride = DefectBuilder(surface).create_oxynitride(
    nitrogen_fraction=0.67,    # 2/3 O → N
    vacancy_concentration=0.1  # 10% vacancies (parameter x)
)
```

### 2. NH3 Decomposition Pathway
```python
from nh3sofc.structure import DecompositionBuilder

decomp = DecompositionBuilder(optimized_nh3_slab)

# Generate intermediates with multiple H configurations
nh2_h_configs = decomp.create_NH2_H_configs(n_configs=10)
nh_2h_configs = decomp.create_NH_2H_configs(n_configs=10)
n_3h_configs = decomp.create_N_3H_configs(n_configs=5)

# Filter duplicates
unique_configs = decomp.filter_unique_configs(all_configs, threshold=0.5)
```

### 3. Surface Comparison
```python
from nh3sofc.analysis import SurfaceComparator, find_rds

surfaces = {
    'LaO-term': {'NH3*': 0.0, 'NH2*+H*': 0.8, 'NH*+2H*': 1.5, 'N*+3H*': 1.2},
    'VO2-term': {'NH3*': 0.0, 'NH2*+H*': 0.6, 'NH*+2H*': 1.2, 'N*+3H*': 0.9},
}

comparator = SurfaceComparator(surfaces)
ranking = comparator.rank_by_energy_span()
comparator.plot_energy_profiles('comparison.png')

rds = find_rds(surfaces['LaO-term'], method='energy_span')
```

### 4. High-Throughput Screening
```python
from nh3sofc.workflows import ScreeningWorkflow

wf = ScreeningWorkflow(
    base_structure=surface,
    parameter_space={
        "vacancy_concentration": [0.0, 0.05, 0.10, 0.15, 0.20],
        "nitrogen_fraction": [0.5, 0.67],
    }
)
wf.setup()
results = wf.run()
```

## Dependencies

```
ase>=3.22.0
numpy>=1.21.0
scipy>=1.7.0
mace-torch>=0.3.0  # optional
spglib>=2.0.0      # optional
catkit             # optional, for advanced site detection
```

## Implementation Order

**Phase A: Foundation**
1. `core/constants.py` + `core/base.py`
2. `structure/bulk.py`
3. `structure/surface.py`
4. `structure/adsorbates.py`

**Phase B: Structure Generators**
5. `structure/defects.py`
6. `structure/decomposition.py`

**Phase C: VASP Integration**
7. `calculators/vasp/inputs.py`
8. `calculators/vasp/outputs.py`
9. `calculators/vasp/frequency.py`
10. `jobs/pbs.py`

**Phase D: Workflows**
11. `workflows/relaxation.py`
12. `workflows/decomposition.py`
13. `workflows/neb.py`
14. `workflows/frequency.py`

**Phase E: Analysis**
15. `analysis/energetics.py`
16. `analysis/thermochemistry.py`
17. `analysis/rds.py`
18. `analysis/surface_comparison.py`

**Phase F: Advanced**
19. `database/`
20. `calculators/mace/`
21. `workflows/screening.py`
22. `analysis/microkinetic.py`
23. `analysis/trajectory.py`

## Verification Plan

1. **Unit tests**: Test each module independently
   - Structure: verify surface creation, defect counts
   - VASP: validate INCAR/KPOINTS generation
   - Adsorbates: check collision detection, site identification

2. **Integration tests**: Full workflow on Ni/YSZ benchmark
   - Build (111) Ni surface
   - Add NH3 adsorbate
   - Generate VASP inputs
   - Verify directory structure and file contents

3. **Validation**: Compare calculated adsorption energies with literature values for Ni-NH3 system
