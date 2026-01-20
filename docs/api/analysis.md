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
