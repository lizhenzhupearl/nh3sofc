# Tutorial: Surface Comparison

This tutorial covers comparing and ranking catalyst surfaces for NH3 decomposition.

## Learning Objectives

- Compare energy profiles across surfaces
- Rank surfaces by activity
- Create volcano plots
- Identify optimal catalysts

## Step 1: Collect Energy Profiles

```python
from nh3sofc.analysis import SurfaceComparator

# Energy profiles from DFT calculations
surfaces = {
    'LaO-termination': {
        'NH3*': 0.0,
        'NH2*+H*': 0.8,
        'NH*+2H*': 1.5,
        'N*+3H*': 1.2,
    },
    'VO2-termination': {
        'NH3*': 0.0,
        'NH2*+H*': 0.6,
        'NH*+2H*': 1.2,
        'N*+3H*': 0.9,
    },
    'LaO-vac10%': {
        'NH3*': 0.0,
        'NH2*+H*': 0.7,
        'NH*+2H*': 1.3,
        'N*+3H*': 1.0,
    },
}

comparator = SurfaceComparator(surfaces)
```

## Step 2: Rank by Energy Span

The energy span model identifies the thermodynamic driving force:

```python
ranking = comparator.rank_by_energy_span()

print("Surface Ranking (by activity):")
for i, (name, span) in enumerate(ranking, 1):
    print(f"  {i}. {name}: δE = {span:.3f} eV")
```

Lower energy span = higher activity.

## Step 3: Rank by Maximum Barrier

Using BEP-estimated barriers:

```python
ranking = comparator.rank_by_max_barrier()

print("\nRanking by Maximum Barrier:")
for i, (name, barrier) in enumerate(ranking, 1):
    print(f"  {i}. {name}: E_a,max = {barrier:.3f} eV")
```

## Step 4: Get Best Surface

```python
# By energy span
best = comparator.get_best_surface(method="energy_span")
print(f"\nBest surface (energy span): {best}")

# By max barrier
best = comparator.get_best_surface(method="max_barrier")
print(f"Best surface (max barrier): {best}")
```

## Step 5: Visualize Comparison

### Energy Profile Comparison

```python
comparator.plot_energy_profiles("energy_profiles.png")
```

### Custom Plot

```python
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(10, 6))

steps = list(list(surfaces.values())[0].keys())
x = np.arange(len(steps))
colors = ['blue', 'red', 'green']

for (name, profile), color in zip(surfaces.items(), colors):
    energies = [profile[s] for s in steps]
    ax.plot(x, energies, 'o-', label=name, color=color, markersize=8)

ax.set_xticks(x)
ax.set_xticklabels(steps, rotation=45, ha='right')
ax.set_ylabel('Energy (eV)')
ax.set_title('NH3 Decomposition Energy Profiles')
ax.legend()
ax.axhline(0, color='gray', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('comparison.png', dpi=150)
```

## Step 6: Volcano Plot

```python
from nh3sofc.analysis import create_volcano_plot

# Descriptor (e.g., N binding energy)
descriptors = {
    'LaO-termination': -0.5,
    'VO2-termination': -0.8,
    'LaO-vac10%': -0.6,
}

# Activity metric (e.g., TOF or -energy_span)
activities = {
    name: -span for name, span in comparator.rank_by_energy_span()
}

create_volcano_plot(
    descriptors,
    activities,
    xlabel='N Binding Energy (eV)',
    ylabel='-Energy Span (eV)',
    filename='volcano.png'
)
```

### Manual Volcano Plot

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 6))

for name in surfaces:
    x = descriptors[name]
    y = activities[name]
    ax.scatter(x, y, s=100, label=name)
    ax.annotate(name, (x, y), xytext=(5, 5), textcoords='offset points')

ax.set_xlabel('N Binding Energy (eV)')
ax.set_ylabel('-Energy Span (eV)')
ax.set_title('Volcano Plot')

plt.tight_layout()
plt.savefig('volcano.png', dpi=150)
```

## Step 7: Statistical Analysis

```python
import pandas as pd

# Create DataFrame
data = []
for name, profile in surfaces.items():
    row = {'surface': name}
    row.update(profile)
    row['energy_span'] = comparator.get_energy_span(name)
    row['max_barrier'] = comparator.get_max_barrier(name)
    data.append(row)

df = pd.DataFrame(data)
print(df.to_string(index=False))

# Save
df.to_csv('surface_comparison.csv', index=False)
```

## Complete Example

```python
from nh3sofc.analysis import SurfaceComparator, RDSAnalyzer
from nh3sofc.workflows import DecompositionWorkflow

# 1. Collect results from multiple workflows
surfaces = {}

for termination in ['LaO', 'VO2']:
    workflow = DecompositionWorkflow(
        nh3_on_slab=slabs[termination],
        work_dir=f"./{termination}",
    )

    # After calculations complete
    profile = workflow.get_energy_profile(reference="NH3")
    surfaces[termination] = profile

# 2. Compare
comparator = SurfaceComparator(surfaces)

# 3. Rank
ranking = comparator.rank_by_energy_span()
print("Ranking:")
for i, (name, span) in enumerate(ranking, 1):
    print(f"  {i}. {name}: {span:.3f} eV")

# 4. RDS for each surface
for name, profile in surfaces.items():
    analyzer = RDSAnalyzer()
    analyzer.set_pathway(profile)
    rds, dE = analyzer.find_rds_thermodynamic()
    print(f"{name} RDS: {rds}")

# 5. Plot comparison
comparator.plot_energy_profiles("comparison.png")

# 6. Best surface
best = comparator.get_best_surface()
print(f"\nRecommended surface: {best}")
```

## Scoring Function

Custom multi-criteria ranking:

```python
def score_surface(profile, barriers, weights=None):
    """Score surface by multiple criteria."""
    if weights is None:
        weights = {
            'energy_span': 0.4,
            'max_barrier': 0.4,
            'n_binding': 0.2,
        }

    scores = {
        'energy_span': -max(profile.values()),
        'max_barrier': -max(barriers.values()),
        'n_binding': -profile.get('N*+3H*', 0),
    }

    total = sum(w * scores[k] for k, w in weights.items())
    return total

# Apply
surface_scores = {}
for name in surfaces:
    profile = surfaces[name]
    barriers = comparator.get_bep_barriers(name)
    surface_scores[name] = score_surface(profile, barriers)

# Rank
ranked = sorted(surface_scores.items(), key=lambda x: -x[1])
print("\nCustom Ranking:")
for name, score in ranked:
    print(f"  {name}: {score:.3f}")
```

## Best Practices

1. **Consistent DFT** - Use same settings for all surfaces
2. **Multiple metrics** - Don't rely on single criterion
3. **Include barriers** - Energy span alone may miss kinetics
4. **Validate** - Compare with experimental data

## Next Steps

- [Microkinetic Modeling](microkinetics.md) - Predict TOF for each surface
- [High-Throughput Screening](screening.md) - Screen more compositions
