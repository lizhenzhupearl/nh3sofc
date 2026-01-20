# Tutorial: NH3 Decomposition Pathway

This tutorial walks through a complete workflow for studying the NH3 decomposition reaction on a catalyst surface.

## Learning Objectives

By the end of this tutorial, you will be able to:

- Generate decomposition intermediate configurations
- Set up batch calculations for all steps
- Calculate energy profiles
- Identify the rate-determining step

## Background

NH3 decomposition follows the co-adsorption pathway:

```
NH3* → NH2* + H* → NH* + 2H* → N* + 3H* → ½N2(g) + 3/2H2(g)
```

Each step involves dissociating one N-H bond and placing H on the surface.

## Step 1: Prepare the Initial Structure

First, we need an optimized NH3 on surface structure:

```python
from ase.io import read
from nh3sofc.structure import AdsorbatePlacer

# Load your relaxed surface
surface = read("relaxed_surface.traj")

# Add NH3
placer = AdsorbatePlacer(surface)
nh3_configs = placer.add_on_site("NH3", site_type="ontop", height=2.0)

# Use the first configuration
nh3_on_slab = nh3_configs[0]
```

**Important**: For accurate decomposition energetics, the initial NH3* structure should be fully optimized.

## Step 2: Generate Decomposition Intermediates

### Using DecompositionBuilder

```python
from nh3sofc.structure import DecompositionBuilder

# Initialize builder with relaxed NH3/surface
decomp = DecompositionBuilder(nh3_on_slab)

# Step 1: NH3* → NH2* + H*
nh2_h_configs = decomp.create_NH2_H_configs(
    n_configs=10,  # Generate 10 configurations
    h_sites="auto" # Automatically find H binding sites
)
print(f"NH2+H configurations: {len(nh2_h_configs)}")

# Step 2: NH2* + H* → NH* + 2H*
nh_2h_configs = decomp.create_NH_2H_configs(n_configs=10)
print(f"NH+2H configurations: {len(nh_2h_configs)}")

# Step 3: NH* + 2H* → N* + 3H*
n_3h_configs = decomp.create_N_3H_configs(n_configs=5)
print(f"N+3H configurations: {len(n_3h_configs)}")
```

### Understanding Configuration Generation

For NH2+H step:
- One H is removed from NH3
- The H is placed at various surface sites (ontop, bridge, hollow)
- Distance constraints ensure physical configurations

```python
# Customize H placement
nh2_h_custom = decomp.create_NH2_H_configs(
    n_configs=10,
    h_height=1.0,           # Height above surface
    min_h_distance=1.5,     # Min distance from N
    max_h_distance=4.0,     # Max distance from N
    random_seed=42          # Reproducibility
)
```

## Step 3: Set Up Calculations

### Using DecompositionWorkflow

```python
from nh3sofc.workflows import DecompositionWorkflow

workflow = DecompositionWorkflow(
    nh3_on_slab=nh3_on_slab,
    work_dir="./decomposition_study",
    n_configs_per_step=5,
    calculator="vasp",
    # VASP parameters
    encut=520,
    kspacing=0.03,
    hubbard_u={"V": 3.25},
    vdw="D3BJ",
    # PBS parameters
    nodes=1,
    ppn=24,
    walltime="24:00:00",
)

# Generate configurations and set up all calculations
file_paths = workflow.setup()

print("Calculation directories created:")
for step, paths in file_paths.items():
    print(f"  {step}: {len(paths)} configurations")
```

### Directory Structure

```
decomposition_study/
├── initial_configs/          # Initial structures
│   ├── NH3.xyz
│   ├── NH2_H_000.xyz
│   └── ...
├── NH3/                      # NH3* calculations
│   └── config_000/
│       ├── INCAR
│       ├── KPOINTS
│       ├── POSCAR
│       ├── POTCAR
│       └── run.pbs
├── NH2_H/                    # NH2* + H* calculations
│   ├── config_000/
│   ├── config_001/
│   └── ...
├── NH_2H/                    # NH* + 2H* calculations
│   └── ...
└── N_3H/                     # N* + 3H* calculations
    └── ...
```

### Submit Jobs

```bash
# Submit all jobs
cd decomposition_study
for dir in */config_*/; do
    cd $dir
    qsub run.pbs
    cd ../..
done
```

## Step 4: Parse Results

After calculations complete:

```python
# Parse all results
results = workflow.parse_results()

# Check status
status = workflow.get_status()
for step, s in status.items():
    print(f"{step}: {s['completed']}/{sum(s.values())} completed")

# Get lowest energies
lowest = workflow.get_lowest_energies()
for step, energy in lowest.items():
    print(f"{step}: {energy:.4f} eV")
```

## Step 5: Calculate Energy Profile

### Relative Energies

```python
# Get energy profile relative to NH3*
profile = workflow.get_energy_profile(reference="NH3")

print("\nEnergy Profile (vs NH3*):")
for step, dE in profile.items():
    print(f"  {step:10s}: {dE:+.3f} eV")
```

### Reaction Energies

```python
# Step-by-step reaction energies
reactions = workflow.get_reaction_energies()

print("\nReaction Energies:")
for name, dE in reactions.items():
    sign = "exothermic" if dE < 0 else "endothermic"
    print(f"  {name:20s}: {dE:+.3f} eV ({sign})")
```

### Plotting

```python
import matplotlib.pyplot as plt
import numpy as np

# Prepare data
steps = list(profile.keys())
energies = [profile[s] for s in steps]

# Create energy diagram
fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(steps))
ax.step(x, energies, where='mid', linewidth=2, color='blue')
ax.scatter(x, energies, s=100, zorder=5)

# Add labels
ax.set_xticks(x)
ax.set_xticklabels(steps, rotation=45, ha='right')
ax.set_ylabel('Energy (eV)', fontsize=12)
ax.set_title('NH3 Decomposition Energy Profile', fontsize=14)
ax.axhline(0, color='gray', linestyle='--', alpha=0.5)

# Annotate reaction energies
for i in range(len(steps)-1):
    dE = energies[i+1] - energies[i]
    color = 'green' if dE < 0 else 'red'
    ax.annotate(f'{dE:+.2f}', xy=(i+0.5, (energies[i]+energies[i+1])/2),
                fontsize=10, color=color, ha='center')

plt.tight_layout()
plt.savefig('energy_profile.png', dpi=150)
```

## Step 6: Identify Rate-Determining Step

### Using RDS Analyzer

```python
from nh3sofc.analysis import RDSAnalyzer

analyzer = RDSAnalyzer()
analyzer.set_pathway(profile)

# Thermodynamic RDS (highest ΔE step)
rds_step, rds_energy = analyzer.find_rds_thermodynamic()
print(f"Thermodynamic RDS: {rds_step} (ΔE = {rds_energy:.3f} eV)")

# BEP-estimated barriers
bep_rds, bep_barrier = analyzer.find_rds_bep()
print(f"BEP-estimated RDS: {bep_rds} (E_a ≈ {bep_barrier:.3f} eV)")

# Energy span model
span, tdts, tdi = analyzer.find_rds_energy_span()
print(f"Energy span: {span:.3f} eV")
print(f"  TDTS: {tdts}")
print(f"  TDI: {tdi}")
```

### Full Analysis Report

```python
analyzer.print_analysis()
```

Output:
```
============================================================
Rate-Determining Step Analysis
============================================================

1. Thermodynamic RDS (highest ΔE):
   Step: NH3*→NH2*+H*
   ΔE = 0.800 eV

2. BEP-estimated RDS (highest E_a):
   Step: NH3*→NH2*+H*
   E_a ≈ 1.646 eV

3. Energy Span Model:
   Energy span: 1.500 eV
   TDTS: TS_NH3*_NH2*+H*
   TDI: NH3*

4. Step-by-step barriers:
   NH3*→NH2*+H*: E_a = 1.646 eV (BEP), ΔE = 0.800 eV
   NH2*+H*→NH*+2H*: E_a = 1.558 eV (BEP), ΔE = 0.700 eV
   NH*+2H*→N*+3H*: E_a = 0.515 eV (BEP), ΔE = -0.300 eV
============================================================
```

## Step 7: Compare with Other Surfaces

```python
from nh3sofc.analysis import SurfaceComparator

# Collect results from multiple surfaces
surfaces = {
    'LaO-term': workflow.get_energy_profile(),
    # Add more surfaces...
}

comparator = SurfaceComparator(surfaces)
ranking = comparator.rank_by_energy_span()

print("\nSurface Ranking (by activity):")
for i, (name, span) in enumerate(ranking, 1):
    print(f"  {i}. {name}: δE = {span:.3f} eV")
```

## Best Practices

### 1. Configuration Sampling
- Generate more configs than needed, then filter
- Use RMSD filtering to remove duplicates
- Include diverse H positions

### 2. Convergence Checks
- Verify all calculations converged
- Check for imaginary frequencies (optimize further if found)
- Ensure forces are below threshold

### 3. Thermochemistry
- Include ZPE corrections for accurate energetics
- Consider temperature effects for reaction conditions

## Next Steps

- [NEB Calculations](neb.md) - Calculate actual transition state barriers
- [Frequency Calculations](thermochemistry.md) - Add thermodynamic corrections
- [Microkinetic Modeling](microkinetics.md) - Predict reaction rates

## Complete Script

```python
"""Complete NH3 decomposition workflow."""

from ase.io import read
from nh3sofc.structure import AdsorbatePlacer
from nh3sofc.workflows import DecompositionWorkflow
from nh3sofc.analysis import RDSAnalyzer, SurfaceComparator

# 1. Load relaxed NH3/surface (from previous calculation)
nh3_on_slab = read("relaxed_nh3_on_surface.traj")

# 2. Set up workflow
workflow = DecompositionWorkflow(
    nh3_on_slab=nh3_on_slab,
    work_dir="./decomposition",
    n_configs_per_step=5,
    encut=520,
    hubbard_u={"V": 3.25},
    vdw="D3BJ",
)

# 3. Setup calculations
workflow.setup()

# 4. After calculations complete...
# results = workflow.parse_results()

# 5. Analyze
# profile = workflow.get_energy_profile()
# workflow.print_summary()
```
