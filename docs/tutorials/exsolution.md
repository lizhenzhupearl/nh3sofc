# Exsolution Simulation

This tutorial covers simulating exsolution processes in perovskite materials, where transition metal cations migrate from the bulk to the surface under reducing conditions to form metallic nanoparticles.

## Background

**Exsolution** is a phenomenon where B-site transition metals (Ni, Co, Fe) in perovskite oxides (e.g., La₀.₄Sr₀.₄Ti₀.₉Ni₀.₁O₃) emerge from the lattice to form metallic nanoparticles on the surface. Key features:

- Occurs under reducing atmospheres (H₂, NH₃)
- Creates **socketed** nanoparticles anchored into the surface
- Generates oxygen vacancies that enhance catalytic activity
- Produces highly stable, sinter-resistant metal particles

## Exsolution Pathway

The exsolution process involves four stages:

```
1. Pristine      → 2. Defective     → 3. Segregated   → 4. Exsolved
   Perovskite       + Vacancies        + Surface B       + Nanoparticle
```

## Quick Start

### Create Exsolution Structures

```python
from nh3sofc.structure import BulkStructure, SurfaceBuilder, ExsolutionBuilder

# Load perovskite structure
bulk = BulkStructure.from_cif("LaSrTiNiO3.cif")

# Create (001) surface
surface = SurfaceBuilder(bulk).create_surface(
    miller_index=(0, 0, 1),
    layers=6,
    vacuum=15.0
)
surface.fix_bottom_layers(2)

# Create exsolution builder
builder = ExsolutionBuilder(surface)

# Generate full exsolution pathway
pathway = builder.create_exsolution_pathway(
    metal="Ni",
    particle_size=13,
    vacancy_fraction=0.1
)

# Save structures
for step in pathway:
    atoms = step["atoms"]
    stage = step["stage"]
    atoms.write(f"{stage}.vasp", format="vasp")
    print(f"{stage}: {step['description']}")
```

### Individual Steps

#### 1. Identify Perovskite Sites

```python
sites = builder.identify_perovskite_sites()
print(f"A-site atoms: {len(sites['A_site'])}")  # La, Sr
print(f"B-site atoms: {len(sites['B_site'])}")  # Ti, Ni
print(f"O-site atoms: {len(sites['O_site'])}")  # O
```

#### 2. Create Defective Perovskite

```python
defective = builder.create_defective_perovskite(
    a_site_vacancy_fraction=0.05,   # 5% A-site vacancies
    b_site_vacancy_fraction=0.1,    # 10% B-site (Ni) vacancies
    oxygen_vacancy_fraction=0.08,   # 8% O vacancies
    b_site_element="Ni",            # Specifically remove Ni
    random_seed=42
)
defective.write("defective.vasp", format="vasp")
```

#### 3. Create Nanoparticle on Surface

```python
exsolved = builder.create_nanoparticle(
    metal="Ni",
    n_atoms=13,                  # Ni13 cluster (magic number)
    shape="hemispherical",       # or "icosahedral"
    position="hollow",           # or "ontop", "bridge", "random"
    interface_distance=2.0,      # Metal-oxide distance
    socketed=True,               # Remove surface atoms (real exsolution)
    random_seed=42
)
exsolved.write("exsolved_Ni13.vasp", format="vasp")
```

## Running Calculations

### Single Exsolution Study

```python
from nh3sofc.workflows import ExsolutionWorkflow

wf = ExsolutionWorkflow(
    atoms=surface.atoms,
    work_dir="./exsolution_Ni13",
    metal="Ni",
    particle_size=13,
    vacancy_fraction=0.1,
    calculator="vasp",
    encut=520,
    hubbard_u={"Ni": 6.2, "Ti": 3.0},
    vdw="D3BJ",
)

# Generate all structures
wf.generate_pathway_structures()

# Setup VASP calculations
wf.setup()

# Submit jobs
# cd exsolution_Ni13 && bash submit_all.sh
```

### Parse Results

```python
# After VASP calculations complete
results = wf.parse_results()

print(f"Pristine energy: {results['pristine']['energy']:.2f} eV")
print(f"Exsolved energy: {results['exsolved']['energy']:.2f} eV")
print(f"Exsolution energy: {results['exsolution_energy']:.2f} eV")
print(f"Favorable: {results['summary']['favorable']}")
```

## Energetics Analysis

### Calculate Exsolution Driving Force

```python
from nh3sofc.analysis import ExsolutionEnergetics

energetics = ExsolutionEnergetics(
    metal="Ni",
    temperature=873.0,  # 600°C
    p_o2=1e-20          # Reducing atmosphere
)

# Set reference
energetics.set_reference_energies(e_pristine=-250.0)

# Add DFT results
energetics.add_stage("pristine", -250.0)
energetics.add_stage("defective", -235.0, n_o_vacancies=4)
energetics.add_stage("segregated", -238.0)
energetics.add_stage("exsolved", -265.0, n_particle_atoms=13, n_o_vacancies=4)

# Calculate thermodynamics
result = energetics.get_exsolution_driving_force()

print(f"ΔG_exsolution = {result['delta_G']:.2f} eV")
print(f"μ_O = {result['mu_O']:.2f} eV")
print(f"Favorable: {result['favorable']}")

# Print summary
energetics.print_summary()
```

### Effect of Temperature and Pressure

```python
import numpy as np

# Vary temperature
for T in [673, 773, 873, 973]:
    result = energetics.get_exsolution_driving_force(temperature=T)
    print(f"T={T}K: ΔG={result['delta_G']:.2f} eV")

# Vary oxygen partial pressure
for log_p in [-25, -20, -15, -10]:
    result = energetics.get_exsolution_driving_force(p_o2=10**log_p)
    print(f"log(p_O2)={log_p}: ΔG={result['delta_G']:.2f} eV")
```

## High-Throughput Screening

### Screen Multiple Metals and Sizes

```python
from nh3sofc.workflows import ExsolutionScreeningWorkflow

screening = ExsolutionScreeningWorkflow(
    base_structure=surface.atoms,
    parameter_space={
        "metal": ["Ni", "Co", "Fe"],
        "particle_size": [1, 4, 13],
        "vacancy_fraction": [0.05, 0.1, 0.15],
    },
    work_dir="./exsolution_screening",
    calculator="vasp",
    encut=520,
)

# Generate all combinations
configs = screening.generate_all()
print(f"Total configurations: {len(configs)}")

# Setup calculations
screening.setup_all()

# After completion:
results = screening.parse_all()
best = screening.get_best_result(metric="exsolution_energy", minimize=True)
print(f"Best configuration: {best['config']}")
```

## NH3 Catalysis on Exsolved Particles

### Identify Adsorption Sites

Exsolved particles have multiple unique adsorption site types:

```python
# Get adsorption sites on exsolved structure
sites = builder.get_adsorption_sites(exsolved)

print(f"Metal top sites: {len(sites['metal_top'])}")
print(f"Interface edge sites: {len(sites['interface_edge'])}")
print(f"Vacancy sites: {len(sites['vacancy_site'])}")
print(f"Oxide surface sites: {len(sites['oxide_surface'])}")
```

### Couple with NH3 Decomposition

```python
# Continue from exsolution workflow
decomp_wf = wf.couple_with_decomposition()
decomp_wf.setup()

# This creates NH3 decomposition pathway on the exsolved particle
```

### Compare with Clean Surface

```python
# Energies from NH3 decomposition on exsolved particle
exsolved_energies = {
    "NH3*": -0.8,
    "NH2*+H*": 0.2,
    "NH*+2H*": 0.5,
    "N*+3H*": 0.3,
}

# Energies from clean perovskite surface
clean_energies = {
    "NH3*": -0.5,
    "NH2*+H*": 0.6,
    "NH*+2H*": 1.2,
    "N*+3H*": 0.9,
}

comparison = energetics.compare_with_clean_surface(
    exsolved_energies,
    clean_energies
)

print("More favorable on exsolved:", comparison["more_favorable_on_exsolved"])
print("More favorable on clean:", comparison["more_favorable_on_clean"])
```

## Key Differences: Exsolved vs. Deposited Catalysts

| Feature | Exsolved | Deposited |
|---------|----------|-----------|
| Anchoring | Socketed into support | Weakly bound |
| Sintering | Resistant | Prone to coalescence |
| Active sites | Metal + interface + vacancies | Metal surface only |
| Coking | Resistant (small particles) | Variable |

## Best Practices

1. **Particle sizes**: Use magic numbers (1, 4, 7, 13, 19) for stable clusters
2. **Vacancy coupling**: ~2 O vacancies per reduced B-site cation
3. **Socket modeling**: Enable `socketed=True` for realistic exsolution
4. **Temperature**: Typical exsolution occurs at 600-900°C
5. **Atmosphere**: Use very low p(O₂) (~10⁻²⁰ atm) for reducing conditions

## Example Directory Structure

```
exsolution_study/
├── pristine/
│   ├── INCAR, POSCAR, KPOINTS, POTCAR
│   └── run.pbs
├── defective/
│   └── ...
├── segregated/
│   └── ...
├── exsolved/
│   └── ...
├── submit_all.sh
└── nh3_decomposition/  # From couple_with_decomposition()
    ├── NH3/
    ├── NH2_H/
    └── ...
```

## Next Steps

- [Thermochemistry](thermochemistry.md) - Calculate Gibbs free energies
- [Microkinetics](microkinetics.md) - Model reaction rates on exsolved particles
- [Surface Comparison](surface_comparison.md) - Compare different catalyst surfaces
