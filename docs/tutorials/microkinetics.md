# Tutorial: Microkinetic Modeling

This tutorial covers microkinetic modeling to predict reaction rates and identify rate-limiting steps.

## Learning Objectives

- Build microkinetic models from DFT data
- Solve for steady-state coverages
- Calculate turnover frequencies (TOF)
- Identify rate-determining steps

## Overview

Microkinetic modeling connects DFT energetics to observable kinetics:

```
DFT Energies → Rate Constants → Coverages → TOF
```

## Step 1: Calculate Rate Constants

```python
from nh3sofc.analysis import RateConstantCalculator

calc = RateConstantCalculator(temperature=673)  # K

# From barrier (Eyring equation)
barrier = 1.0  # eV
k = calc.from_barrier(barrier)
print(f"Rate constant: {k:.2e} s^-1")

# Include entropy (more accurate)
k = calc.from_barrier(
    barrier=1.0,
    delta_S=-0.0005  # eV/K (entropy change)
)
```

## Step 2: Build Microkinetic Model

### Manual Model Construction

```python
from nh3sofc.analysis import MicroKineticModel

model = MicroKineticModel(temperature=673, total_sites=1.0)

# Add surface species
model.add_species("*", is_site=True)       # Empty site
model.add_species("NH3*")
model.add_species("NH2*")
model.add_species("H*")
model.add_species("N*")

# Add gas species with partial pressures (bar)
model.add_gas("NH3(g)", pressure=0.1)
model.add_gas("H2(g)", pressure=0.01)
model.add_gas("N2(g)", pressure=0.01)

# Add reactions with barriers
model.add_reaction(
    name="NH3_adsorption",
    reactants={"NH3(g)": 1, "*": 1},
    products={"NH3*": 1},
    Ea_fwd=0.0,   # Adsorption barrier
    dG=-0.5       # Free energy change
)

model.add_reaction(
    name="NH3_dissociation",
    reactants={"NH3*": 1, "*": 1},
    products={"NH2*": 1, "H*": 1},
    Ea_fwd=1.2,   # Forward barrier
    Ea_rev=0.4    # Reverse barrier
)

# ... add more reactions
```

### Using Pre-built NH3 Decomposition Model

```python
from nh3sofc.analysis import NH3DecompositionModel

# Default barriers
model = NH3DecompositionModel(temperature=673)

# Or custom barriers from your DFT
barriers = {
    "NH3_ads": 0.0,
    "NH3_NH2": 1.0,
    "NH2_NH": 0.9,
    "NH_N": 0.7,
    "N_N2": 1.3,
    "H_H2": 0.5,
}

model = NH3DecompositionModel(
    temperature=673,
    barriers=barriers
)
```

## Step 3: Solve Steady State

```python
coverages = model.solve_steady_state()

print("Steady-state coverages:")
for species, theta in coverages.items():
    print(f"  {species}: {theta:.4f}")
```

## Step 4: Calculate TOF

```python
tof = model.get_tof()
print(f"Turnover frequency: {tof:.4e} s^-1")

# Per reaction
reaction_rates = model.get_reaction_rates(coverages)
for rxn, rate in reaction_rates.items():
    print(f"  {rxn}: {rate:.4e} s^-1")
```

## Step 5: Sensitivity Analysis

```python
# Degree of rate control
drc = model.degree_of_rate_control()

print("Degree of Rate Control:")
for rxn, chi in sorted(drc.items(), key=lambda x: -abs(x[1])):
    print(f"  {rxn}: {chi:.3f}")
```

## Step 6: Temperature Dependence

```python
import numpy as np
import matplotlib.pyplot as plt

temperatures = np.linspace(500, 900, 50)
tofs = []

for T in temperatures:
    model.set_temperature(T)
    coverages = model.solve_steady_state()
    tofs.append(model.get_tof(coverages))

plt.figure(figsize=(8, 5))
plt.semilogy(1000/temperatures, tofs, 'b-', linewidth=2)
plt.xlabel('1000/T (K$^{-1}$)')
plt.ylabel('TOF (s$^{-1}$)')
plt.title('Arrhenius Plot')
plt.savefig('arrhenius.png', dpi=150)
```

## Complete Example

```python
from nh3sofc.analysis import NH3DecompositionModel
import numpy as np
import matplotlib.pyplot as plt

# 1. Create model with custom barriers
barriers = {
    "NH3_ads": 0.0,
    "NH3_NH2": 1.0,  # From NEB
    "NH2_NH": 0.85,
    "NH_N": 0.70,
    "N_N2": 1.25,
    "H_H2": 0.50,
}

model = NH3DecompositionModel(
    temperature=673,
    barriers=barriers
)

# 2. Solve steady state
coverages = model.solve_steady_state()

print("Steady-state coverages:")
for species, theta in coverages.items():
    if theta > 1e-6:
        print(f"  {species}: {theta:.4f}")

# 3. Calculate TOF
tof = model.get_tof()
print(f"\nTOF: {tof:.4e} s^-1")

# 4. Rate control analysis
drc = model.degree_of_rate_control()
rds = max(drc, key=drc.get)
print(f"\nRate-determining step: {rds}")

# 5. Apparent activation energy
from nh3sofc.core.constants import KB_EV

T1, T2 = 650, 700
model.set_temperature(T1)
tof1 = model.get_tof()
model.set_temperature(T2)
tof2 = model.get_tof()

Ea_app = -KB_EV * np.log(tof2/tof1) / (1/T2 - 1/T1)
print(f"Apparent Ea: {Ea_app:.2f} eV")
```

## Advanced: Pressure Dependence

```python
import numpy as np

p_NH3_range = np.logspace(-3, 0, 20)  # 0.001 to 1 bar
tofs = []

for p in p_NH3_range:
    model.set_pressure("NH3(g)", p)
    coverages = model.solve_steady_state()
    tofs.append(model.get_tof(coverages))

# Reaction order
from scipy.stats import linregress
slope, _, _, _, _ = linregress(np.log(p_NH3_range), np.log(tofs))
print(f"Reaction order in NH3: {slope:.2f}")
```

## Best Practices

1. **Consistent energies** - Use same DFT settings for all barriers
2. **Include ZPE** - Correct barriers for zero-point energy
3. **Check convergence** - Verify steady-state is reached
4. **Validate** - Compare with experimental data when available

## Next Steps

- [Surface Comparison](surface_comparison.md) - Compare TOF across surfaces
- [Thermochemistry](thermochemistry.md) - Improve barriers with thermodynamic corrections
