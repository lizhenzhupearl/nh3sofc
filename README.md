# nh3sofc

A Python package for automating DFT and ML force field calculations in NH3 solid oxide fuel cell research.

## Project Summary

**nh3sofc** is a simulation automation package designed for computational materials science research on ammonia-based solid oxide fuel cells (SOFCs). The package focuses on:

- **Perovskite oxynitrides** (e.g., LaVO3 → LaVON₂₋ₓ) with tunable vacancy concentration
- **NH3 decomposition catalysis** on electrode surfaces
- **High-throughput screening** of material compositions and defect concentrations

### Key Capabilities

| Module | Description |
|--------|-------------|
| `structure/` | Build surfaces, create oxynitride defects, place adsorbates |
| `calculators/` | Interface with VASP (DFT) and MACE (ML force field) |
| `workflows/` | Automated relaxation, NEB, decomposition pathways |
| `analysis/` | Thermochemistry, rate-determining step, surface comparison |
| `database/` | ASE database integration with standardized naming |

### Target Materials

- **Primary**: Perovskite oxynitrides (LaVON₂₋ₓ, LaVO₃, etc.)
- **Benchmark**: Ni/YSZ cermet (traditional SOFC anode)
- **Future**: Non-stoichiometric oxynitrides from ABO₃, ABO₄ oxides

### Computational Tools

- **DFT**: VASP with Hubbard U, vdW-D3 corrections
- **MLFF**: MACE for accelerated calculations and active learning
- **HPC**: PBS/Torque job management

## Installation

```bash
pip install -e .
```

### Dependencies

```
ase>=3.22.0
numpy>=1.21.0
scipy>=1.7.0
mace-torch>=0.3.0  # optional
spglib>=2.0.0      # optional
```

## Quick Start

### 1. Build an Oxynitride Surface

```python
from nh3sofc.structure import BulkStructure, SurfaceBuilder, DefectBuilder

# Load bulk structure
bulk = BulkStructure.from_cif("LaVO3.cif")

# Create (001) surface
surface = SurfaceBuilder(bulk).create_surface(
    miller_index=(0, 0, 1),
    layers=6,
    vacuum=15.0
)

# Create oxynitride with vacancies
oxynitride = DefectBuilder(surface).create_oxynitride(
    nitrogen_fraction=0.67,      # 2/3 of O → N
    vacancy_concentration=0.1    # 10% vacancies (parameter x)
)
```

### 2. Place NH3 on Surface

```python
from nh3sofc.structure import AdsorbatePlacer

placer = AdsorbatePlacer(surface)

# Site-specific placement
configs = placer.add_on_site("NH3", site_type="ontop", atom_types=["La", "V"])

# Or generate all configurations
all_configs = placer.add_at_all_sites("NH3", site_type="ontop")
```

### 3. Generate VASP Inputs

```python
from nh3sofc.calculators.vasp import VASPInputGenerator

vasp = VASPInputGenerator(
    atoms=slab_with_nh3,
    calc_type="relax",
    work_dir="./calculations/LaVON_001_NH3"
)

vasp.generate_all(
    encut=520,
    hubbard_u={"V": 3.25},
    vdw="D3BJ",
    is_surface=True
)
```

### 4. Analyze NH3 Decomposition

```python
from nh3sofc.analysis import SurfaceComparator, find_rds

# Energy profile from calculations
energies = {
    'NH3*': 0.0,
    'NH2*+H*': 0.8,
    'NH*+2H*': 1.5,
    'N*+3H*': 1.2
}

# Find rate-determining step
rds = find_rds(energies, method='energy_span')
print(f"RDS: {rds['step']}, Barrier: {rds['barrier']:.2f} eV")
```

## NH3 Decomposition Workflow

The package implements a complete workflow for studying NH3 decomposition:

```
NH3(g) → NH3* → NH2* + H* → NH* + 2H* → N* + 3H* → N2(g) + 3H2(g)
```

### Co-adsorption Approach

For each decomposition step, multiple H configurations are generated:

```python
from nh3sofc.structure import DecompositionBuilder

decomp = DecompositionBuilder(optimized_nh3_slab)

# Step 1: NH3* → NH2* + H*
nh2_h_configs = decomp.create_NH2_H_configs(n_configs=10)

# Step 2: NH2*+H* → NH*+2H*
nh_2h_configs = decomp.create_NH_2H_configs(n_configs=10)

# Step 3: NH*+2H* → N*+3H*
n_3h_configs = decomp.create_N_3H_configs(n_configs=5)
```

### Thermochemistry

Go beyond 0K electronic energies:

```python
from nh3sofc.analysis import ThermochemistryCalculator

thermo = ThermochemistryCalculator()

# Calculate Gibbs free energy at reaction temperature
dG = thermo.calculate_gibbs_energy(
    atoms=structure,
    temperature=673,  # K (400°C typical SOFC operating temp)
    include_zpe=True,
    include_entropy=True
)
```

## Documentation

**Full documentation:** [https://lizhenzhupearl.github.io/nh3sofc](https://lizhenzhupearl.github.io/nh3sofc)

- [Installation Guide](docs/installation.md) - Detailed installation instructions
- [Quick Start](docs/quickstart.md) - Get started quickly
- [Tutorials](docs/tutorials/) - Step-by-step guides
- [API Reference](docs/api/) - Complete API documentation

### Build Documentation Locally

```bash
pip install mkdocs mkdocs-material mkdocstrings[python]
mkdocs serve  # Preview at http://127.0.0.1:8000
mkdocs build  # Build static site
```

## Project Structure

```
NH3SOFC/
├── nh3sofc/                 # Main package
│   ├── structure/           # Structure building
│   ├── calculators/         # VASP, MACE interfaces
│   ├── workflows/           # Automated workflows
│   ├── analysis/            # Post-processing
│   └── database/            # Data management
├── examples/                # Usage examples
├── tests/                   # Unit tests
├── plan.md                  # Implementation plan
├── nh3sofc.md              # Workflow documentation
└── README.md               # This file
```

## Contributing

This package is under active development. Key areas for contribution:

1. **Structure module**: Additional defect types, surface reconstructions
2. **MACE integration**: Training workflows, uncertainty quantification
3. **Analysis**: Microkinetic modeling, descriptor development
4. **Testing**: Unit tests, validation against literature

## License

MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This package builds on:
- [ASE](https://wiki.fysik.dtu.dk/ase/) - Atomic Simulation Environment
- [MACE](https://github.com/ACEsuit/mace) - ML force field
- [VASP](https://www.vasp.at/) - Vienna Ab initio Simulation Package
