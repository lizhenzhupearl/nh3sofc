# Tutorial: MACE ML Force Fields

This tutorial covers using MACE machine learning force fields for fast calculations.

## Learning Objectives

- Use MACE foundation models
- Train custom MACE models
- Implement active learning workflows

## Overview

MACE provides ML force fields that are:
- 100-1000x faster than DFT
- Accurate for similar systems
- Useful for screening and exploration

## Step 1: Using Foundation Models

```python
from nh3sofc.calculators.mace import MACECalculatorWrapper

# Load foundation model
mace = MACECalculatorWrapper(
    foundation_model="medium",  # "small", "medium", or "large"
    device="auto"  # "cpu", "cuda", or "auto"
)

# Get ASE calculator
calc = mace.get_calculator()

# Attach to atoms
atoms.calc = calc
energy = atoms.get_potential_energy()
forces = atoms.get_forces()

print(f"Energy: {energy:.4f} eV")
```

## Step 2: Fast Relaxation

```python
from ase.optimize import BFGS

atoms.calc = mace.get_calculator()

opt = BFGS(atoms, trajectory="relax.traj")
opt.run(fmax=0.03)

print(f"Final energy: {atoms.get_potential_energy():.4f} eV")
```

## Step 3: Fast NEB

```python
from nh3sofc.workflows import NEBWorkflow

neb = NEBWorkflow(
    initial=initial,
    final=final,
    work_dir="./mace_neb",
    n_images=7,
)

# Run with MACE (much faster than VASP)
images = neb.run_mace(fmax=0.03, steps=500)

results = neb.parse_results()
print(f"Barrier: {results['barrier_forward']:.3f} eV")
```

## Step 4: Training Custom Models

### Extract Training Data

```python
from nh3sofc.calculators.mace import TrainingDataExtractor

extractor = TrainingDataExtractor()

# Extract from VASP calculations
for calc_dir in ["./calc1", "./calc2", "./calc3"]:
    extractor.extract_from_directory(
        calc_dir,
        include_trajectory=True  # Include ionic steps
    )

training_data = extractor.extract_all()
print(f"Extracted {len(training_data)} configurations")

# Filter high-energy configs
filtered = extractor.filter_by_energy(
    training_data,
    max_energy_per_atom=-2.0
)

# Save
extractor.write_xyz("training_data.xyz", filtered)
```

### Generate Training Config

```python
from nh3sofc.calculators.mace import MACETrainingConfig

config = MACETrainingConfig(
    train_file="train.xyz",
    valid_file="valid.xyz",
    model_name="LaVON_MACE"
)

config.generate_config(
    r_max=6.0,
    max_epochs=500,
    batch_size=5,
)

config.write_config("mace_config.yaml")
config.generate_training_script("train.sh")
```

### Train Model

```bash
# Run training script
bash train.sh

# Or manually
mace_run_train --config mace_config.yaml
```

## Step 5: Uncertainty Estimation

```python
from nh3sofc.calculators.mace import MACEEnsemble

# Ensemble of models
ensemble = MACEEnsemble([
    "./model_1.model",
    "./model_2.model",
    "./model_3.model",
])

results = ensemble.calculate_with_uncertainty(atoms)

print(f"Energy: {results['energy']:.4f} ± {results['energy_std']:.4f} eV")
print(f"Max force std: {results['max_force_std']:.4f} eV/Å")

# Flag high uncertainty for DFT recalculation
if results['max_force_std'] > 0.1:
    print("High uncertainty - add to training data!")
```

## Step 6: Active Learning

```python
from nh3sofc.calculators.mace import ActiveLearningWorkflow

al = ActiveLearningWorkflow(
    initial_training_data="train.xyz",
    ensemble_size=3,
    uncertainty_threshold=0.1,  # eV/Å
    work_dir="./active_learning",
)

# Run active learning loop
for iteration in range(5):
    # 1. Train ensemble
    al.train_ensemble()

    # 2. Explore with MACE
    new_configs = al.explore(
        base_structure=surface,
        n_configs=100,
        method="md",  # or "random"
        temperature=673,
    )

    # 3. Select high-uncertainty configs
    selected = al.select_for_dft(new_configs, n_select=10)

    # 4. Run DFT on selected
    al.run_dft_calculations(selected)

    # 5. Update training data
    al.update_training_data()

    print(f"Iteration {iteration}: Added {len(selected)} new configs")
```

## Complete Example

```python
from ase.io import read
from nh3sofc.calculators.mace import MACECalculatorWrapper, TrainingDataExtractor
from nh3sofc.workflows import NEBWorkflow
from ase.optimize import BFGS

# 1. Quick relaxation with foundation model
atoms = read("structure.xyz")
mace = MACECalculatorWrapper(foundation_model="medium")

atoms.calc = mace.get_calculator()
opt = BFGS(atoms)
opt.run(fmax=0.05)

# 2. Use as starting point for VASP
atoms.write("mace_relaxed.xyz")

# 3. Or run NEB with MACE
initial = read("initial.traj")
final = read("final.traj")

neb = NEBWorkflow(initial, final, n_images=5)
images = neb.run_mace(fmax=0.05)

print(f"Quick barrier estimate: {neb.parse_results()['barrier_forward']:.2f} eV")
```

## Best Practices

1. **Foundation models** - Good for similar chemistries
2. **Custom training** - Needed for unusual systems
3. **Validate with DFT** - Always check key results
4. **Active learning** - Most efficient for new systems
5. **Ensemble uncertainty** - Detect extrapolation

## Next Steps

- [Microkinetic Modeling](microkinetics.md) - Use ML-computed data for kinetics
- [Screening](screening.md) - ML-accelerated screening
