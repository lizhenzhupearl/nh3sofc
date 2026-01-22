# Theoretical Framework

This page describes the theoretical methods implemented in NH3SOFC for catalysis modeling and analysis.

## Overview

NH3SOFC implements a comprehensive theoretical framework for computational catalysis, covering the complete workflow from electronic structure analysis to microkinetic modeling.

```
Electronic Structure → Descriptors → Energetics → Thermochemistry → Kinetics → Activity
```

## Electronic Structure Descriptors

### D-Band Center Analysis

The d-band model (Hammer-Norskov) is fundamental for understanding adsorption on transition metals. The d-band center is defined as:

$$\varepsilon_d = \frac{\int \varepsilon \cdot \rho_d(\varepsilon) \, d\varepsilon}{\int \rho_d(\varepsilon) \, d\varepsilon}$$

where $\rho_d(\varepsilon)$ is the d-band density of states.

```python
from nh3sofc.analysis import DBandAnalyzer, calculate_d_band_center

# Parse VASP DOSCAR and analyze d-band
analyzer = DBandAnalyzer.from_doscar("DOSCAR")

# Get d-band center for surface atoms
d_center = analyzer.get_d_band_center(atom_index=0)
print(f"D-band center: {d_center:.3f} eV (relative to Fermi)")

# Get comprehensive d-band properties
props = analyzer.get_d_band_properties([0, 1, 2, 3])
for atom, data in props.items():
    print(f"Atom {atom}: center={data['center']:.2f} eV, "
          f"width={data['width']:.2f} eV, filling={data['filling']:.2f}")

# Average d-band center for surface atoms
avg_center = analyzer.get_surface_average_d_band_center([0, 1, 2, 3])
```

**Available descriptors:**
- D-band center (1st moment)
- D-band width (2nd moment)
- D-band filling

## Adsorption Energetics

### Adsorption Energy

$$E_{ads} = E_{ads/surf} - E_{surf} - E_{gas}$$

```python
from nh3sofc.analysis import AdsorptionEnergyCalculator

calc = AdsorptionEnergyCalculator()
calc.set_surface_energy(-250.0)
calc.set_gas_reference("NH3", -19.54)

E_ads = calc.calculate(e_total=-270.5, adsorbate="NH3")
print(f"Adsorption energy: {E_ads:.3f} eV")
```

### Coverage-Dependent Energies

Differential adsorption energy as a function of coverage:

$$E_{ads}^{diff}(\theta) = \frac{\partial E_{total}}{\partial n_{ads}}$$

```python
from nh3sofc.analysis import calculate_coverage_dependent_energy

# energies at different coverages
energies = [-270.5, -290.8, -310.5]  # 1, 2, 3 adsorbates
coverages = [0.25, 0.5, 0.75]

diff_energies = calculate_coverage_dependent_energy(energies, coverages)
```

## Thermochemistry

### Harmonic Approximation (Adsorbates)

For surface-bound species, we use the harmonic approximation:

$$G = E_{DFT} + ZPE + \int_0^T C_v \, dT - TS_{vib}$$

where the zero-point energy and vibrational entropy are calculated from vibrational frequencies.

```python
from nh3sofc.analysis import HarmonicThermo

# Frequencies from VASP frequency calculation (cm⁻¹)
frequencies = [3450, 3350, 1620, 1100, 580, 450]

thermo = HarmonicThermo(frequencies, e_electronic=-5.23)
G = thermo.get_gibbs_energy(temperature=700)  # K
print(f"Gibbs free energy: {G:.3f} eV")
```

### Ideal Gas Thermochemistry

For gas-phase species:

$$G = E_{DFT} + ZPE + H(T) - TS(T, p)$$

Including translational, rotational, and vibrational contributions.

```python
from nh3sofc.analysis import IdealGasThermo

thermo = IdealGasThermo(
    frequencies=[3450, 3350, 1620, 1100, 580, 450],
    e_electronic=-19.54,
    geometry="nonlinear",
    spin=0,
    mass=17.03  # NH3
)

G = thermo.get_gibbs_energy(temperature=700, pressure=1e5)
```

## Reaction Kinetics

### BEP Relations

Brønsted-Evans-Polanyi relations connect activation barriers to reaction energies:

$$E_a = \alpha \cdot \Delta E + \beta$$

```python
from nh3sofc.analysis import BEPRelation

# Built-in BEP parameters for N-H dissociation
bep = BEPRelation.for_reaction("N-H_dissociation")

# Estimate barrier from reaction energy
delta_E = 0.5  # eV (endothermic)
barrier = bep.estimate_barrier(delta_E)
print(f"Estimated barrier: {barrier:.2f} eV")

# Or fit from your own data
custom_bep = BEPRelation.fit_from_data(
    reaction_energies=[0.1, 0.3, 0.5, 0.7],
    barriers=[0.8, 0.9, 1.0, 1.1]
)
```

### Energy Span Model (Kozuch-Shaik)

The energy span model identifies the TOF-determining transition state (TDTS) and TOF-determining intermediate (TDI):

$$\delta E = \begin{cases}
T_{TDTS} - I_{TDI} & \text{if TDTS after TDI} \\
T_{TDTS} - I_{TDI} + \Delta G_r & \text{if TDTS before TDI}
\end{cases}$$

```python
from nh3sofc.analysis import EnergySpanModel

# Define intermediates and transition states
intermediates = {
    "NH3*": 0.0,
    "NH2*+H*": 0.3,
    "NH*+2H*": 0.5,
    "N*+3H*": 0.2,
}

transition_states = {
    "TS1": 0.8,   # NH3* → NH2*+H*
    "TS2": 1.1,   # NH2*+H* → NH*+2H*
    "TS3": 0.9,   # NH*+2H* → N*+3H*
}

model = EnergySpanModel(intermediates, transition_states, delta_G_rxn=-0.5)

# Get energy span and TOF-determining species
span = model.get_energy_span()
tdi, tdts = model.get_determining_states()
print(f"Energy span: {span:.2f} eV")
print(f"TOF-determining intermediate: {tdi}")
print(f"TOF-determining transition state: {tdts}")
```

### Rate-Determining Step Analysis

Multiple methods for RDS identification:

```python
from nh3sofc.analysis import RDSAnalyzer

analyzer = RDSAnalyzer(
    intermediates=intermediates,
    transition_states=transition_states
)

# Method 1: Thermodynamic (highest reaction energy)
rds_thermo = analyzer.find_rds_thermodynamic()

# Method 2: BEP-estimated barriers
rds_bep = analyzer.find_rds_bep()

# Method 3: Energy span model
rds_span = analyzer.find_rds_energy_span()

# Summary
analyzer.print_summary()
```

## Microkinetic Modeling

### Steady-State Surface Kinetics

Solve for steady-state coverages using the mean-field approximation:

$$\frac{d\theta_i}{dt} = \sum_j \nu_{ij} r_j = 0$$

```python
from nh3sofc.analysis import MicroKineticModel

model = MicroKineticModel()

# Define elementary steps
model.add_step("NH3_ads", barrier_fwd=0.0, barrier_rev=0.5,
               reactants=["*"], products=["NH3*"])
model.add_step("NH3_diss", barrier_fwd=0.8, barrier_rev=0.5,
               reactants=["NH3*", "*"], products=["NH2*", "H*"])
# ... more steps

# Solve for coverages at reaction conditions
coverages = model.solve_steady_state(temperature=700, pressures={"NH3": 0.1})

# Calculate turnover frequency
tof = model.calculate_tof(coverages)
```

### NH3 Decomposition Model

Pre-configured model for NH3 decomposition with default parameters:

```python
from nh3sofc.analysis import NH3DecompositionModel

model = NH3DecompositionModel()

# Use default barriers or customize
model.set_barrier("NH3_diss", 0.85)

# Solve and analyze
coverages = model.solve(temperature=700)
tof = model.get_tof(coverages)

# Sensitivity analysis
sensitivity = model.degree_of_rate_control()
```

## Catalyst Screening

### Volcano Plots

Visualize activity as a function of binding energy descriptor:

```python
from nh3sofc.analysis import SurfaceComparator, ActivityDescriptor

# Compare multiple catalyst surfaces
comparator = SurfaceComparator()

comparator.add_surface("Ni(111)",
    intermediates={"NH3*": -0.8, "NH2*": -0.5, ...},
    barriers={"TS1": 0.9, ...})
comparator.add_surface("Pt(111)", ...)
comparator.add_surface("Ru(0001)", ...)

# Generate volcano plot
comparator.plot_volcano(
    descriptor="NH3*",  # Use NH3 binding as descriptor
    activity_metric="energy_span"
)

# Find optimal catalyst
best = comparator.get_best_catalyst(metric="energy_span")
print(f"Best catalyst: {best}")
```

### Activity Prediction from Descriptors

```python
from nh3sofc.analysis import ActivityDescriptor

descriptor = ActivityDescriptor()

# Fit scaling relation from known data
descriptor.fit_scaling_relation(
    descriptor_values=[-0.8, -0.6, -0.4, -0.2],  # NH3 binding
    activities=[1e5, 1e6, 1e7, 1e6]  # TOF
)

# Predict activity for new catalyst
predicted_activity = descriptor.predict_activity(descriptor_value=-0.5)
```

## Exsolution Thermodynamics

Specialized analysis for metal exsolution from perovskites:

```python
from nh3sofc.analysis import ExsolutionEnergetics

exsol = ExsolutionEnergetics()

# Calculate exsolution driving force
driving_force = exsol.calculate_driving_force(
    e_pristine=-300.0,
    e_defective=-298.5,
    e_segregated=-299.0,
    temperature=1073,  # K
    p_O2=1e-20  # atm
)

# Vacancy formation energy
e_vac = exsol.vacancy_formation_energy(
    e_defective=-298.5,
    e_pristine=-300.0,
    n_vacancies=1
)
```

## Complete Workflow Example

```python
from nh3sofc.analysis import (
    DBandAnalyzer,
    AdsorptionEnergyCalculator,
    HarmonicThermo,
    NH3DecompositionModel,
    SurfaceComparator
)

# 1. Electronic structure analysis
d_analyzer = DBandAnalyzer.from_doscar("DOSCAR")
d_center = d_analyzer.get_surface_average_d_band_center([0, 1, 2, 3])

# 2. Adsorption energies
calc = AdsorptionEnergyCalculator()
calc.set_surface_energy(e_surface)
E_NH3 = calc.calculate(e_nh3_surface, "NH3")

# 3. Free energy corrections
thermo = HarmonicThermo(frequencies, E_NH3)
G_NH3 = thermo.get_gibbs_energy(temperature=700)

# 4. Microkinetic modeling
model = NH3DecompositionModel()
model.set_barrier("NH3_diss", barrier_from_dft)
tof = model.get_tof(model.solve(temperature=700))

# 5. Catalyst comparison
comparator = SurfaceComparator()
comparator.add_surface("My_catalyst", intermediates, barriers)
ranking = comparator.rank_by("energy_span")
```

## References

The theoretical methods implemented in this package are based on:

1. **D-band model**: Hammer, B. & Norskov, J.K. *Adv. Catal.* 45, 71-129 (2000)
2. **BEP relations**: Michaelides, A. et al. *J. Am. Chem. Soc.* 125, 3704 (2003)
3. **Energy span model**: Kozuch, S. & Shaik, S. *Acc. Chem. Res.* 44, 101 (2011)
4. **Microkinetic modeling**: Dumesic, J.A. et al. *The Microkinetics of Heterogeneous Catalysis* (1993)
5. **Computational catalysis**: Norskov, J.K. et al. *Nat. Chem.* 1, 37 (2009)

## Next Steps

- [Microkinetics Tutorial](microkinetics.md) - Detailed microkinetic modeling examples
- [Surface Comparison](surface_comparison.md) - Catalyst screening workflows
- [Exsolution Tutorial](exsolution.md) - Exsolution energetics analysis
