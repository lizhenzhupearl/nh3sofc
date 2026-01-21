# Analysis Module

The `nh3sofc.analysis` module provides tools for analyzing calculation results.

## AdsorptionEnergyCalculator

Calculate adsorption energies.

```python
from nh3sofc.analysis import AdsorptionEnergyCalculator
```

### Constructor

```python
AdsorptionEnergyCalculator(
    e_surface: float = None,
    gas_references: dict = None
)
```

### Methods

#### calculate

```python
calculate(
    e_total: float,
    adsorbate: str = "NH3",
    n_adsorbates: int = 1
) -> float
```

Calculate adsorption energy: `E_ads = E(ads/surf) - E(surf) - n*E(gas)`

**Example:**

```python
calc = AdsorptionEnergyCalculator()
calc.set_surface_energy(-150.0)
calc.set_gas_reference("NH3", -19.54)

E_ads = calc.calculate(-170.0, "NH3")
print(f"Adsorption energy: {E_ads:.3f} eV")
```

---

## HarmonicThermo

Harmonic thermodynamics for adsorbates.

```python
from nh3sofc.analysis import HarmonicThermo
```

### Constructor

```python
HarmonicThermo(
    frequencies: List[float],   # cm^-1
    electronic_energy: float = 0.0
)
```

### Methods

| Method | Description |
|--------|-------------|
| `get_zpe()` | Zero-point energy (eV) |
| `get_vibrational_energy(T)` | U_vib at temperature T (eV) |
| `get_vibrational_entropy(T)` | S_vib at temperature T (eV/K) |
| `get_gibbs_energy(T)` | Gibbs free energy (eV) |
| `get_heat_capacity(T)` | Cv at temperature T (eV/K) |

**Example:**

```python
# Frequencies from VASP frequency calculation
frequencies = [3400, 3300, 1600, 1100, 500, 400]  # cm^-1

thermo = HarmonicThermo(frequencies, electronic_energy=-170.0)

print(f"ZPE: {thermo.get_zpe():.4f} eV")
print(f"G(673K): {thermo.get_gibbs_energy(673):.4f} eV")
```

---

## IdealGasThermo

Ideal gas thermodynamics including translational, rotational, and vibrational contributions.

```python
from nh3sofc.analysis import IdealGasThermo
```

### Constructor

```python
IdealGasThermo(
    frequencies: List[float],
    electronic_energy: float,
    mass: float,              # amu
    geometry: str = "nonlinear",
    symmetry_number: int = 1,
    spin: int = 0
)
```

### Methods

| Method | Description |
|--------|-------------|
| `get_translational_entropy(T, p)` | S_trans (eV/K) |
| `get_rotational_entropy(T)` | S_rot (eV/K) |
| `get_total_entropy(T, p)` | Total entropy (eV/K) |
| `get_enthalpy_correction(T)` | H(T) - H(0) (eV) |
| `get_gibbs_energy(T, p)` | Gibbs free energy (eV) |

**Example:**

```python
# H2 gas at 673 K, 1 bar
h2_thermo = IdealGasThermo(
    frequencies=[4400],  # H-H stretch
    electronic_energy=-6.77,
    mass=2.016,
    geometry="linear",
    symmetry_number=2
)

G_H2 = h2_thermo.get_gibbs_energy(T=673, p=1.0)
```

---

## RDSAnalyzer

Identify the rate-determining step.

```python
from nh3sofc.analysis import RDSAnalyzer
```

### Constructor

```python
RDSAnalyzer()
```

### Methods

#### set_pathway

```python
set_pathway(pathway: Dict[str, float])
```

Set reaction pathway energies.

#### find_rds_thermodynamic

```python
find_rds_thermodynamic() -> Tuple[str, float]
```

Find RDS using thermodynamic approximation (highest ΔE step).

#### find_rds_bep

```python
find_rds_bep() -> Tuple[str, float]
```

Find RDS using BEP-estimated barriers.

#### find_rds_energy_span

```python
find_rds_energy_span() -> Tuple[float, str, str]
```

Find RDS using energy span model.

**Returns:** `(energy_span, TDTS, TDI)`

**Example:**

```python
analyzer = RDSAnalyzer()
analyzer.set_pathway({
    'NH3*': 0.0,
    'NH2*+H*': 0.8,
    'NH*+2H*': 1.5,
    'N*+3H*': 1.2
})

# Thermodynamic RDS
rds, dE = analyzer.find_rds_thermodynamic()
print(f"RDS: {rds}, ΔE = {dE:.3f} eV")

# Energy span
span, tdts, tdi = analyzer.find_rds_energy_span()
print(f"Energy span: {span:.3f} eV")
```

---

## BEPRelation

Brønsted-Evans-Polanyi relations for barrier estimation.

```python
from nh3sofc.analysis import BEPRelation
```

### Constructor

```python
BEPRelation(
    alpha: float = None,
    beta: float = None,
    reaction_type: str = "default"
)
```

### Pre-defined Reaction Types

| Type | α | β | Description |
|------|---|---|-------------|
| `"N-H_dissociation"` | 0.87 | 0.95 | NH3 → NH2 + H |
| `"C-H_dissociation"` | 0.75 | 0.90 | C-H breaking |
| `"O-H_dissociation"` | 0.65 | 0.85 | O-H breaking |
| `"N-N_formation"` | 0.50 | 1.20 | N2 formation |
| `"H-H_formation"` | 0.60 | 0.80 | H2 formation |
| `"default"` | 0.80 | 1.00 | General |

### Methods

#### estimate_barrier

```python
estimate_barrier(reaction_energy: float) -> float
```

Estimate activation barrier: `E_a = α*ΔE + β`

**Example:**

```python
bep = BEPRelation(reaction_type="N-H_dissociation")
barrier = bep.estimate_barrier(0.5)  # ΔE = 0.5 eV
print(f"Estimated barrier: {barrier:.3f} eV")
```

---

## SurfaceComparator

Compare and rank catalyst surfaces.

```python
from nh3sofc.analysis import SurfaceComparator
```

### Constructor

```python
SurfaceComparator(
    surfaces: Dict[str, Dict[str, float]],
    barriers: Dict[str, Dict[str, float]] = None
)
```

### Methods

#### rank_by_energy_span

```python
rank_by_energy_span() -> List[Tuple[str, float]]
```

Rank surfaces by energy span (lower is better).

#### rank_by_max_barrier

```python
rank_by_max_barrier() -> List[Tuple[str, float]]
```

Rank surfaces by maximum barrier.

#### get_best_surface

```python
get_best_surface(method: str = "energy_span") -> str
```

Get the best performing surface.

#### plot_energy_profiles

```python
plot_energy_profiles(filename: str = None)
```

Plot comparison of energy profiles.

**Example:**

```python
surfaces = {
    'LaO-term': {'NH3*': 0.0, 'NH2*+H*': 0.8, 'NH*+2H*': 1.5, 'N*+3H*': 1.2},
    'VO2-term': {'NH3*': 0.0, 'NH2*+H*': 0.6, 'NH*+2H*': 1.2, 'N*+3H*': 0.9},
}

comparator = SurfaceComparator(surfaces)

# Ranking
ranking = comparator.rank_by_energy_span()
for i, (name, span) in enumerate(ranking, 1):
    print(f"{i}. {name}: {span:.3f} eV")

# Best surface
best = comparator.get_best_surface()
print(f"Best: {best}")

# Plot
comparator.plot_energy_profiles("comparison.png")
```

---

## MicroKineticModel

Microkinetic modeling for steady-state analysis.

```python
from nh3sofc.analysis import MicroKineticModel
```

### Constructor

```python
MicroKineticModel(
    temperature: float = 673.0,
    total_sites: float = 1.0
)
```

### Methods

#### add_species

```python
add_species(
    name: str,
    is_site: bool = False,
    initial_coverage: float = 0.0
)
```

Add a surface species.

#### add_gas

```python
add_gas(name: str, pressure: float = 1.0)
```

Add a gas phase species with partial pressure.

#### add_reaction

```python
add_reaction(
    name: str,
    reactants: Dict[str, int],
    products: Dict[str, int],
    Ea_fwd: float,
    Ea_rev: float = None,
    dG: float = None
)
```

Add an elementary reaction.

#### solve_steady_state

```python
solve_steady_state() -> Dict[str, float]
```

Solve for steady-state coverages.

#### get_tof

```python
get_tof(
    product_reaction: str = None,
    coverages: Dict[str, float] = None
) -> float
```

Calculate turnover frequency.

**Example:**

```python
model = MicroKineticModel(temperature=673)

# Add species
model.add_species("*", is_site=True)
model.add_species("NH3*")
model.add_species("NH2*")
model.add_species("H*")

# Add gas
model.add_gas("NH3(g)", pressure=0.1)
model.add_gas("H2(g)", pressure=0.01)

# Add reactions
model.add_reaction(
    "NH3_ads",
    reactants={"NH3(g)": 1, "*": 1},
    products={"NH3*": 1},
    Ea_fwd=0.0,
    dG=-0.5
)

# Solve
coverages = model.solve_steady_state()
tof = model.get_tof()

print(f"TOF: {tof:.4e} s^-1")
```

---

## NH3DecompositionModel

Pre-built microkinetic model for NH3 decomposition.

```python
from nh3sofc.analysis import NH3DecompositionModel
```

### Constructor

```python
NH3DecompositionModel(
    temperature: float = 673.0,
    barriers: Dict[str, float] = None,
    reaction_energies: Dict[str, float] = None
)
```

**Default barriers:** NH3_ads=0.0, NH3_NH2=1.2, NH2_NH=1.0, NH_N=0.8, N_N2=1.5, H_H2=0.6 eV

**Example:**

```python
# With custom barriers from DFT
barriers = {
    "NH3_ads": 0.0,
    "NH3_NH2": 1.0,  # From NEB
    "NH2_NH": 0.9,
    "NH_N": 0.7,
    "N_N2": 1.3,
    "H_H2": 0.5
}

model = NH3DecompositionModel(temperature=673, barriers=barriers)
tof = model.get_tof()
print(f"TOF: {tof:.4e} s^-1")
```

---

## Convenience Functions

### calculate_adsorption_energy

```python
calculate_adsorption_energy(
    e_adsorbate_surface: float,
    e_surface: float,
    e_adsorbate_gas: float
) -> float
```

### calculate_zpe

```python
calculate_zpe(frequencies: List[float]) -> float
```

### find_thermodynamic_rds

```python
find_thermodynamic_rds(pathway: Dict[str, float]) -> Tuple[str, float]
```

### estimate_barrier_bep

```python
estimate_barrier_bep(
    reaction_energy: float,
    reaction_type: str = "default"
) -> float
```

### calculate_rate_constant

```python
calculate_rate_constant(
    barrier: float,
    temperature: float,
    prefactor: float = 1e13
) -> float
```

### calculate_tof

```python
calculate_tof(
    barriers: Dict[str, float],
    reaction_energies: Dict[str, float],
    temperature: float = 673.0
) -> float
```

---

## ExsolutionEnergetics

Analyze energetics of exsolution processes.

```python
from nh3sofc.analysis import ExsolutionEnergetics
```

### Constructor

```python
ExsolutionEnergetics(
    metal: str = "Ni",
    temperature: float = 873.0,
    p_o2: float = 1e-20
)
```

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `metal` | `str` | `"Ni"` | Exsolution metal |
| `temperature` | `float` | `873.0` | Temperature (K) |
| `p_o2` | `float` | `1e-20` | O2 partial pressure (atm) |

### Methods

#### set_reference_energies

```python
set_reference_energies(
    e_pristine: float = None,
    e_bulk_metal_per_atom: float = None,
    e_o2: float = None
)
```

Set reference energies for calculations.

#### add_stage

```python
add_stage(
    stage: str,
    energy: float,
    n_particle_atoms: int = 0,
    n_o_vacancies: int = 0
)
```

Add energy for a calculation stage.

#### get_oxygen_chemical_potential

```python
get_oxygen_chemical_potential(
    temperature: float = None,
    p_o2: float = None
) -> float
```

Calculate μ_O(T, p) = 0.5 * [E(O₂) + μ°(T) + kT*ln(p/p°)]

#### calculate_vacancy_formation_energy

```python
calculate_vacancy_formation_energy(
    e_defective: float = None,
    e_pristine: float = None,
    n_vacancies: int = 1,
    include_chemical_potential: bool = True
) -> float
```

E_vac = [E(defective) - E(pristine) + n * μ_O] / n

#### calculate_segregation_energy

```python
calculate_segregation_energy(
    e_segregated: float = None,
    e_bulk_distributed: float = None
) -> float
```

E_seg = E(surface_segregated) - E(bulk_distributed)

#### calculate_exsolution_energy

```python
calculate_exsolution_energy(
    e_exsolved: float = None,
    e_substrate: float = None,
    n_atoms: int = None,
    e_bulk_metal: float = None,
    include_vacancy_correction: bool = True
) -> float
```

Calculate exsolution driving force.

#### calculate_particle_binding_energy

```python
calculate_particle_binding_energy(
    e_system: float,
    e_surface: float,
    e_isolated_particle: float
) -> float
```

E_bind = E(particle/surface) - E(surface) - E(isolated_particle)

#### get_exsolution_driving_force

```python
get_exsolution_driving_force(
    temperature: float = None,
    p_o2: float = None
) -> Dict[str, float]
```

Calculate T and p dependent driving force.

**Returns:**

```python
{
    "delta_G": float,
    "delta_E": float,
    "vacancy_term": float,
    "entropy_term": float,
    "mu_O": float,
    "favorable": bool
}
```

#### compare_with_clean_surface

```python
compare_with_clean_surface(
    exsolved_energies: Dict[str, float],
    clean_surface_energies: Dict[str, float]
) -> Dict[str, Any]
```

Compare catalytic activity of exsolved vs clean surface.

#### print_summary

```python
print_summary()
```

Print formatted energetics summary.

**Example:**

```python
energetics = ExsolutionEnergetics(metal="Ni", temperature=873)

# Set references
energetics.set_reference_energies(e_pristine=-250.0)

# Add calculation results
energetics.add_stage("pristine", -250.0)
energetics.add_stage("defective", -235.0, n_o_vacancies=4)
energetics.add_stage("exsolved", -265.0, n_particle_atoms=13, n_o_vacancies=4)

# Calculate
result = energetics.get_exsolution_driving_force()
print(f"ΔG_exsolution = {result['delta_G']:.2f} eV")
print(f"Favorable: {result['favorable']}")

energetics.print_summary()
```

---

## calculate_exsolution_driving_force

Convenience function for quick calculation.

```python
from nh3sofc.analysis import calculate_exsolution_driving_force

result = calculate_exsolution_driving_force(
    e_pristine=-250.0,
    e_exsolved=-265.0,
    n_particle_atoms=13,
    metal="Ni",
    n_o_vacancies=4,
    temperature=873.0,
    p_o2=1e-20
)
print(f"ΔG = {result['delta_G']:.2f} eV")
```
