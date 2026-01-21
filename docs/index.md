# NH3-SOFC Simulation Package

**Automated DFT and ML force field calculations for ammonia solid oxide fuel cell research**

---

## Overview

**nh3sofc** is a Python package that automates computational materials science workflows for studying ammonia (NH3) decomposition catalysis on solid oxide fuel cell (SOFC) electrode materials. The package integrates with VASP for DFT calculations and MACE for machine learning force fields.

<div class="grid cards" markdown>

-   :material-cube-outline: **Structure Building**

    ---

    Build surfaces from bulk crystals, create oxynitride defects, and place adsorbate molecules with multiple methods.

    [:octicons-arrow-right-24: Structure Module](api/structure.md)

-   :material-calculator: **VASP & MACE**

    ---

    Generate VASP input files with proper settings for relaxation, NEB, and frequency calculations. Use MACE for fast ML predictions.

    [:octicons-arrow-right-24: Calculators](api/calculators.md)

-   :material-chart-line: **Analysis Tools**

    ---

    Calculate adsorption energies, Gibbs free energies, identify rate-determining steps, and compare catalyst surfaces.

    [:octicons-arrow-right-24: Analysis](api/analysis.md)

-   :material-database: **Data Management**

    ---

    Standardized naming conventions and ASE database integration for organizing calculations.

    [:octicons-arrow-right-24: Database](api/database.md)

-   :material-atom: **Exsolution Simulation**

    ---

    Model metal nanoparticle exsolution from perovskites, including defect formation, segregation, and catalysis on exsolved particles.

    [:octicons-arrow-right-24: Exsolution Tutorial](tutorials/exsolution.md)

</div>

## Key Features

### Target Materials

| Material Class | Examples | Focus |
|---------------|----------|-------|
| Perovskite Oxynitrides | LaVON₂₋ₓ, LaVO₃, Ni/GDC | Primary research target |
| Exsolution Perovskites | La₀.₄Sr₀.₄Ti₀.₉Ni₀.₁O₃ | Metal nanoparticle catalysis |
| Mixed Ion Conductors | LSGM, BCZY | Proton/oxide conductors |
| Cermet Anodes | Ni/YSZ, Ni/GDC | Benchmark materials |

### NH3 Decomposition Pathway

The package implements automated workflows for studying the complete NH3 decomposition mechanism:

```
NH3(g) → NH3* → NH2* + H* → NH* + 2H* → N* + 3H* → ½N2(g) + 3/2H2(g)
```

### Exsolution Simulation

Model the complete exsolution pathway for perovskite materials:

```
Pristine → Defective → Segregated → Exsolved (with nanoparticle)
```

- Support for Ni, Co, Fe exsolution from B-site
- Socketed nanoparticle modeling (realistic anchoring)
- Oxygen vacancy coupling with metal reduction
- Multiple adsorption sites: metal top, interface edge, vacancy sites

### Computational Methods

- **DFT**: VASP with Hubbard U corrections, vdW-D3/D3BJ
- **ML Force Fields**: MACE foundation models and custom training
- **Thermochemistry**: ZPE, entropy, Gibbs free energy corrections
- **Kinetics**: BEP relations, energy span model, microkinetic modeling
- **Exsolution**: Vacancy formation, segregation, and particle binding energetics

## Quick Example

```python
from nh3sofc.structure import BulkStructure, SurfaceBuilder, DefectBuilder
from nh3sofc.calculators.vasp import VASPInputGenerator

# Build oxynitride surface
bulk = BulkStructure.from_cif("LaVO3.cif")
surface = SurfaceBuilder(bulk).create_surface((0,0,1), layers=6, vacuum=15)
oxynitride = DefectBuilder(surface).create_oxynitride(
    nitrogen_fraction=0.67,
    vacancy_concentration=0.1
)

# Generate VASP inputs
vasp = VASPInputGenerator(oxynitride, calc_type="relax", work_dir="./calc")
vasp.generate_all(encut=520, hubbard_u={"V": 3.25}, vdw="D3BJ")
```

## Installation

```bash
# Clone the repository
git clone https://github.com/lizhenzhupearl/nh3sofc.git
cd nh3sofc

# Install in development mode
pip install -e .
```

See [Installation Guide](installation.md) for detailed instructions.

## Documentation

<div class="grid" markdown>

[:material-book-open-variant: **Tutorials**](tutorials/index.md){ .md-button }

Step-by-step guides for common workflows

[:material-api: **API Reference**](api/index.md){ .md-button }

Complete API documentation

[:material-frequently-asked-questions: **Examples**](examples/index.md){ .md-button }

Code examples and use cases

</div>

## Package Structure

```
nh3sofc/
├── core/           # Constants, base classes
├── structure/      # Bulk, surface, defects, adsorbates
├── calculators/
│   ├── vasp/       # VASP input/output handling
│   └── mace/       # MACE ML force field interface
├── jobs/           # PBS/HPC job management
├── workflows/      # Automated calculation workflows
├── analysis/       # Energetics, thermochemistry, kinetics
└── database/       # Naming conventions, ASE DB integration
```

## Citation

If you use this package, please cite:

```bibtex
@software{nh3sofc,
  title = {NH3-SOFC: Automated Simulations for Ammonia Fuel Cell Research},
  author = {Li, Zhenzhu},
  year = {2026},
  url = {https://github.com/lizhenzhupearl/nh3sofc}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/lizhenzhupearl/nh3sofc/blob/main/LICENSE) file for details.
