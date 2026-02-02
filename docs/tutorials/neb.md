# Tutorial: NEB Transition States

This tutorial covers finding transition states using the Nudged Elastic Band (NEB) method.

## Learning Objectives

- Set up NEB calculations
- Use climbing image NEB
- Extract barrier heights
- Validate transition states

## Overview

NEB finds the minimum energy path between two states:

```
Initial State → TS → Final State
     ↓           ↓        ↓
   NH3*         TS      NH2* + H*
```

## Step 1: Prepare Endpoints

```python
from ase.io import read

# Load optimized initial and final structures (POSCAR format)
initial = read("work/decomposition/LaVO3_001/NH3/relaxed.vasp", format="vasp")  # NH3*
final = read("work/decomposition/LaVO3_001/NH2_H/relaxed.vasp", format="vasp")  # NH2* + H*
```

**Important:** Both endpoints must be fully optimized (converged relaxation).

## Step 2: Set Up NEB

### Using NEBWorkflow

```python
from nh3sofc.workflows import NEBWorkflow

neb = NEBWorkflow(
    initial=initial,
    final=final,
    work_dir="./neb_nh3_dissociation",
    n_images=7,
    calculator="vasp",
    encut=520,
    kspacing=0.03,
)

# Setup with IDPP interpolation
neb.setup(
    interpolation="idpp",  # Image Dependent Pair Potential
    climb=True,            # Climbing image NEB
)
```

### Directory Structure

```
neb_nh3_dissociation/
├── 00/          # Initial image
│   ├── POSCAR
│   └── ...
├── 01/          # Intermediate image 1
├── 02/          # ...
├── ...
├── 07/          # Final image
└── run.pbs
```

## Step 3: Run with MACE (Quick Test)

Before expensive VASP calculations, test with MACE:

```python
# Run NEB with MACE ML force field
images = neb.run_mace(
    fmax=0.03,   # Force convergence (eV/Å)
    steps=500,   # Max optimization steps
)

# Check results
results = neb.parse_results()
print(f"Forward barrier: {results['barrier_forward']:.3f} eV")
print(f"Reverse barrier: {results['barrier_reverse']:.3f} eV")
```

## Step 4: Parse VASP Results

After VASP NEB completes:

```python
results = neb.parse_results()

print(f"Forward barrier: {results['barrier_forward']:.3f} eV")
print(f"Reverse barrier: {results['barrier_reverse']:.3f} eV")
print(f"TS image index: {results['ts_image']}")

# Get transition state structure and save as POSCAR
from nh3sofc import write_poscar
ts = results["ts_structure"]
write_poscar(ts, "work/neb/NH3_to_NH2_H/ts.vasp")
```

## Step 5: Visualize

```python
neb.plot_energy_profile("neb_profile.png")
```

```python
import matplotlib.pyplot as plt
import numpy as np

energies = results["energies"]
x = np.arange(len(energies))

plt.figure(figsize=(8, 5))
plt.plot(x, energies, 'o-', markersize=10)
plt.xlabel("Reaction Coordinate")
plt.ylabel("Energy (eV)")
plt.title("NEB Energy Profile")
plt.axhline(0, color='gray', linestyle='--')
plt.savefig("neb_profile.png", dpi=150)
```

## Step 6: Validate TS (Frequency)

A true transition state has exactly one imaginary frequency:

```python
from nh3sofc.workflows import FrequencyWorkflow

ts = results["ts_structure"]

freq = FrequencyWorkflow(
    ts,
    work_dir="./ts_freq",
    adsorbate_indices=ts_atom_indices,
)

freq.setup()

# After calculation...
freq_results = freq.parse_results()

imaginary = freq_results["imaginary"]
if len(imaginary) == 1:
    print(f"Valid TS! Imaginary frequency: {imaginary[0]:.1f}i cm^-1")
else:
    print(f"Warning: {len(imaginary)} imaginary frequencies")
```

## Best Practices

1. **Optimize endpoints** - Ensure both initial and final states are fully relaxed
2. **Number of images** - Use 5-9 images (7 is typical)
3. **IDPP interpolation** - Better than linear for complex paths
4. **Climbing image** - Always use for accurate barrier
5. **Check forces** - Final forces should be < 0.05 eV/Å
6. **Validate with frequency** - Confirm exactly one imaginary mode

## Complete Example

```python
from ase.io import read
from nh3sofc.workflows import NEBWorkflow, FrequencyWorkflow
from nh3sofc import write_poscar

# 1. Load endpoints (POSCAR format)
initial = read("work/decomposition/LaVO3_001/NH3/relaxed.vasp", format="vasp")
final = read("work/decomposition/LaVO3_001/NH2_H/relaxed.vasp", format="vasp")

# 2. Setup NEB with meaningful work directory
neb = NEBWorkflow(
    initial=initial,
    final=final,
    work_dir="work/neb/LaVO3_001_NH3_to_NH2_H",
    n_images=7,
    encut=520,
    hubbard_u={"V": 3.25},
)

neb.setup(interpolation="idpp", climb=True)

# 3. Quick test with MACE
images = neb.run_mace(fmax=0.05)

# 4. Parse results
results = neb.parse_results()
print(f"Barrier: {results['barrier_forward']:.3f} eV")

# 5. Save transition state
ts = results["ts_structure"]
write_poscar(ts, "work/neb/LaVO3_001_NH3_to_NH2_H/ts.vasp")

# 6. Plot
neb.plot_energy_profile("work/neb/LaVO3_001_NH3_to_NH2_H/barrier.png")
```

## Next Steps

- [Thermochemistry](thermochemistry.md) - Add ZPE corrections to barriers
- [Microkinetic Modeling](microkinetics.md) - Use barriers for rate predictions
