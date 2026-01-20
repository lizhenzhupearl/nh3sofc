# Tutorial: Frequency & Thermochemistry

This tutorial covers vibrational frequency calculations and thermodynamic corrections.

## Learning Objectives

- Calculate vibrational frequencies
- Compute zero-point energy (ZPE)
- Calculate Gibbs free energy at reaction temperature

## Background

At finite temperature, electronic energy alone is insufficient. We need:

```
G(T) = E_elec + ZPE + H(T) - TS(T)
```

## Step 1: Set Up Frequency Calculation

```python
from ase.io import read
from nh3sofc.workflows import FrequencyWorkflow

# Load optimized structure
atoms = read("nh3_on_surface.traj")

# Identify adsorbate atoms (e.g., NH3 indices)
adsorbate_indices = [48, 49, 50, 51]

# Set up frequency calculation
freq = FrequencyWorkflow(
    atoms,
    work_dir="./frequency",
    adsorbate_indices=adsorbate_indices,
    encut=520,
)

freq.setup(nfree=2)  # Central difference
```

## Step 2: Parse Frequencies

After VASP completes:

```python
results = freq.parse_results()

print("Frequencies (cm^-1):")
for f in results["real_frequencies"]:
    print(f"  {f:.1f}")

if results["imaginary"]:
    print(f"WARNING: Imaginary frequencies: {results['imaginary']}")
    print("Structure may not be a minimum!")

print(f"\nZPE: {results['zpe']:.4f} eV")
```

## Step 3: Thermodynamic Properties

### For Adsorbates (Harmonic Approximation)

```python
from nh3sofc.analysis import HarmonicThermo

# Create thermo object
thermo = HarmonicThermo(
    frequencies=results["real_frequencies"],
    electronic_energy=-170.0  # From relaxation
)

# At reaction temperature (673 K = 400°C)
T = 673

print(f"ZPE: {thermo.get_zpe():.4f} eV")
print(f"U_vib(T): {thermo.get_vibrational_energy(T):.4f} eV")
print(f"S_vib(T): {thermo.get_vibrational_entropy(T):.6f} eV/K")
print(f"G(T): {thermo.get_gibbs_energy(T):.4f} eV")
```

### For Gas Molecules (Ideal Gas)

```python
from nh3sofc.analysis import IdealGasThermo

# NH3 gas
nh3_thermo = IdealGasThermo(
    frequencies=[3400, 3337, 1627, 1627, 968, 968],  # cm^-1
    electronic_energy=-19.54,
    mass=17.031,  # amu
    geometry="nonlinear",
    symmetry_number=3,
    spin=0
)

T = 673
p = 1.0  # bar

G_NH3 = nh3_thermo.get_gibbs_energy(T, p)
print(f"G(NH3, 673K, 1bar): {G_NH3:.4f} eV")
```

### H2 Gas Example

```python
h2_thermo = IdealGasThermo(
    frequencies=[4400],  # H-H stretch
    electronic_energy=-6.77,
    mass=2.016,
    geometry="linear",
    symmetry_number=2,
    spin=0
)

G_H2 = h2_thermo.get_gibbs_energy(673, 1.0)
```

## Step 4: Gibbs Energy Corrections

### Adsorption Free Energy

```python
# ΔG_ads = G(ads/surf) - G(surf) - G(gas)
dG_ads = G_nh3_adsorbed - G_surface - G_NH3_gas
```

### Reaction Free Energy

```python
# NH3* → NH2* + H*
dG_rxn = (G_NH2_H_adsorbed - G_NH3_adsorbed)

# Temperature effect on equilibrium
from nh3sofc.core.constants import KB_EV
K_eq = np.exp(-dG_rxn / (KB_EV * T))
```

## Complete Example

```python
from ase.io import read
from nh3sofc.workflows import FrequencyWorkflow
from nh3sofc.analysis import HarmonicThermo, IdealGasThermo

# 1. Calculate frequencies for adsorbed NH3
atoms = read("nh3_on_surface.traj")

freq = FrequencyWorkflow(
    atoms,
    work_dir="./freq_nh3",
    adsorbate_indices=[48, 49, 50, 51],
)
freq.setup()

# 2. After VASP completes
results = freq.parse_results()
E_elec = -170.0  # From relaxation

# 3. Get Gibbs energy at 673 K
thermo = HarmonicThermo(
    frequencies=results["real_frequencies"],
    electronic_energy=E_elec
)

G_ads = thermo.get_gibbs_energy(673)
print(f"G(NH3*, 673K): {G_ads:.4f} eV")

# 4. Compare with gas phase
nh3_gas = IdealGasThermo(
    frequencies=[3400, 3337, 1627, 1627, 968, 968],
    electronic_energy=-19.54,
    mass=17.031,
    geometry="nonlinear",
    symmetry_number=3,
)

G_gas = nh3_gas.get_gibbs_energy(673, 1.0)

# 5. Free energy of adsorption
G_surface = -150.0  # From surface calculation
dG_ads = G_ads - G_surface - G_gas
print(f"ΔG_ads(673K): {dG_ads:.4f} eV")
```

## Best Practices

1. **Only vibrate adsorbate** - Fix surface atoms to save computation
2. **Check for imaginary frequencies** - Indicates structure is not a minimum
3. **Use same functional** - Frequencies should match relaxation settings
4. **Temperature consistency** - Use same T throughout analysis

## Common Issues

| Issue | Solution |
|-------|----------|
| Imaginary frequencies | Re-optimize structure with tighter convergence |
| Very low frequencies | May indicate weak binding or translation modes |
| Missing modes | Check all adsorbate atoms are included |

## Next Steps

- [Microkinetic Modeling](microkinetics.md) - Use thermochemistry for rate modeling
- [Surface Comparison](surface_comparison.md) - Compare surfaces at reaction conditions
