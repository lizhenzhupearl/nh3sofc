# Computation inventory (To be toolnized)

***We can use ase to perform necessary management of the lattice, cut the surface, put molecule on the surface, run simulations**

## **STEP 1: build the slabs**

### **STEP 1.1 build the surface**

`from ase.io import read, write
from ase.build import surface`

`*# Read bulk structure*`

`bulk = read('LaVO3.cif')`

`*# Create (001) surface with LaO termination*`

`*# You need to specify the 'layers' parameter carefully*`

`slab_LaO = surface(bulk, (0, 0, 1), layers=5, vacuum=15)`

### **STEP 1.2 place the molecules**

Yes! ASE has several ways to place molecules on surfaces, from simple manual placement to more sophisticated random/systematic placement. Let me show you the options:

## Method 1: Simple `add_adsorbate()` (Manual Placement)

```python
from ase.io import read, write
from ase.build import add_adsorbate, molecule

# Load your surface slab
slab = read('LaVO3_001.cif')

# Create NH3 molecule
nh3 = molecule('NH3')

# Add at specific position
add_adsorbate(slab, nh3,
              height=2.0,        # Distance above surface (Å)
              position=(5, 5),   # xy coordinates on surface
              mol_index=0)       # Which atom in molecule to anchor

write('slab_with_NH3.cif', slab)

```

**Limitations**: You need to manually specify positions - not truly random.

## Method 2: Random Placement (Custom Function)

```python
import numpy as np
from ase.io import read, write
from ase.build import molecule
from ase import Atoms

def add_molecule_randomly(slab, molecule_name, height=2.0, n_attempts=10):
    """
    Randomly place a molecule on the surface

    Parameters:
    -----------
    slab : ASE Atoms object
        The surface slab
    molecule_name : str
        Name of molecule (e.g., 'NH3', 'H2O')
    height : float
        Height above surface (Å)
    n_attempts : int
        Number of random configurations to generate
    """

    # Get surface dimensions
    cell = slab.get_cell()
    a = cell[0, 0]  # x-dimension
    b = cell[1, 1]  # y-dimension

    # Find top of slab (highest z-coordinate)
    z_top = slab.positions[:, 2].max()

    configs = []

    for i in range(n_attempts):
        # Create fresh copy of slab
        slab_copy = slab.copy()

        # Create molecule
        mol = molecule(molecule_name)

        # Random position on surface
        x_rand = np.random.uniform(0, a)
        y_rand = np.random.uniform(0, b)
        z_pos = z_top + height

        # Random rotation (Euler angles)
        alpha = np.random.uniform(0, 2*np.pi)
        beta = np.random.uniform(0, np.pi)
        gamma = np.random.uniform(0, 2*np.pi)

        # Rotate molecule
        mol.rotate(alpha, 'z')
        mol.rotate(beta, 'y')
        mol.rotate(gamma, 'x')

        # Translate molecule to position
        mol.translate([x_rand, y_rand, z_pos])

        # Combine slab + molecule
        combined = slab_copy + mol

        configs.append(combined)

        # Save
        write(f'config_{i:03d}.cif', combined)

    return configs

# Usage
slab = read('LaVO3_001_slab.cif')
configs = add_molecule_randomly(slab, 'NH3', height=2.5, n_attempts=20)

```

## Method 3: Grid-Based Systematic Placement

```python
import numpy as np
from ase.io import read, write
from ase.build import molecule

def place_molecule_on_grid(slab, molecule_name, height=2.0,
                           grid_density=3, orientations=4):
    """
    Place molecule at grid points with different orientations

    Parameters:
    -----------
    grid_density : int
        Number of grid points along each axis
    orientations : int
        Number of random orientations per grid point
    """

    cell = slab.get_cell()
    a, b = cell[0, 0], cell[1, 1]
    z_top = slab.positions[:, 2].max()

    # Create grid
    x_grid = np.linspace(0, a, grid_density)
    y_grid = np.linspace(0, b, grid_density)

    configs = []
    count = 0

    for x in x_grid:
        for y in y_grid:
            for _ in range(orientations):
                slab_copy = slab.copy()
                mol = molecule(molecule_name)

                # Random orientation
                mol.rotate(np.random.uniform(0, 2*np.pi), 'z')
                mol.rotate(np.random.uniform(0, np.pi), 'y')
                mol.rotate(np.random.uniform(0, 2*np.pi), 'x')

                # Position
                mol.translate([x, y, z_top + height])

                combined = slab_copy + mol
                configs.append(combined)

                write(f'grid_config_{count:03d}.cif', combined)
                count += 1

    return configs

# Usage
slab = read('LaVO3_001_slab.cif')
configs = place_molecule_on_grid(slab, 'NH3', height=2.5,
                                  grid_density=4, orientations=3)
print(f"Generated {len(configs)} configurations")

```

## Method 4: Site-Specific Placement (Most Physically Meaningful)

```python
from ase.io import read, write
from ase.build import molecule
import numpy as np

def place_on_surface_atoms(slab, molecule_name, height=2.0,
                           atom_types=['La', 'V'], n_orient=5):
    """
    Place molecule above specific atom types with random orientations

    Parameters:
    -----------
    atom_types : list
        Types of surface atoms to place molecule above
    n_orient : int
        Number of orientations per site
    """

    # Find surface atoms (top 20% of slab)
    z_positions = slab.positions[:, 2]
    z_cutoff = z_positions.max() - 0.2 * (z_positions.max() - z_positions.min())

    surface_indices = [i for i, z in enumerate(z_positions) if z > z_cutoff]

    # Filter by atom type
    target_indices = [i for i in surface_indices
                     if slab[i].symbol in atom_types]

    print(f"Found {len(target_indices)} surface sites: {atom_types}")

    configs = []
    count = 0

    for idx in target_indices:
        site_pos = slab.positions[idx]

        for _ in range(n_orient):
            slab_copy = slab.copy()
            mol = molecule(molecule_name)

            # Random orientation
            mol.rotate(np.random.uniform(0, 2*np.pi), 'z')
            mol.rotate(np.random.uniform(0, np.pi), 'y')
            mol.rotate(np.random.uniform(0, 2*np.pi), 'x')

            # Place above the site
            mol.translate([site_pos[0], site_pos[1], site_pos[2] + height])

            combined = slab_copy + mol
            configs.append(combined)

            write(f'site_{slab[idx].symbol}_{count:03d}.cif', combined)
            count += 1

    return configs

# Usage - place NH3 above La and V sites
slab = read('LaVO3_001_slab.cif')
configs = place_on_surface_atoms(slab, 'NH3', height=2.5,
                                  atom_types=['La', 'V'], n_orient=3)

```

## Method 5: With Collision Detection

```python
from ase.io import read, write
from ase.build import molecule
import numpy as np

def add_molecule_smart(slab, molecule_name, height=2.0,
                       min_distance=2.0, max_attempts=100):
    """
    Add molecule with collision detection

    Parameters:
    -----------
    min_distance : float
        Minimum distance between molecule and surface atoms (Å)
    max_attempts : int
        Maximum tries to find valid position
    """

    cell = slab.get_cell()
    a, b = cell[0, 0], cell[1, 1]
    z_top = slab.positions[:, 2].max()

    for attempt in range(max_attempts):
        # Random position
        x = np.random.uniform(0, a)
        y = np.random.uniform(0, b)
        z = z_top + height

        mol = molecule(molecule_name)

        # Random orientation
        mol.rotate(np.random.uniform(0, 2*np.pi), 'z')
        mol.rotate(np.random.uniform(0, np.pi), 'y')

        mol.translate([x, y, z])

        # Check distances to surface atoms
        distances = []
        for slab_atom in slab:
            for mol_atom in mol:
                dist = np.linalg.norm(slab_atom.position - mol_atom.position)
                distances.append(dist)

        # If no collision, accept
        if min(distances) > min_distance:
            combined = slab + mol
            print(f"Success on attempt {attempt + 1}")
            return combined

    print(f"Failed to place molecule after {max_attempts} attempts")
    return None

# Usage
slab = read('LaVO3_001_slab.cif')
result = add_molecule_smart(slab, 'NH3', height=2.5, min_distance=1.8)
if result:
    write('NH3_on_surface.cif', result)

```

## Method 6: Using External Package (Catkit)

```python
# Catkit has more sophisticated tools for this
# pip install catkit

from catkit.gen.adsorption import AdsorptionSites
from ase.io import read

slab = read('LaVO3_001_slab.cif')

# Find adsorption sites automatically
sites = AdsorptionSites(slab)
coordinates = sites.get_coordinates()

print(f"Found {len(coordinates)} adsorption sites:")
print(f"  Top sites: {len(sites.get_sites('top'))}")
print(f"  Bridge sites: {len(sites.get_sites('bridge'))}")
print(f"  Hollow sites: {len(sites.get_sites('hollow'))}")

# Add molecule to each site
from catkit.gen.adsorption import Builder
builder = Builder(slab)

for i, coord in enumerate(coordinates[:5]):  # First 5 sites
    config = builder.add_adsorbate('NH3', coordinate=coord, height=2.0)
    write(f'catkit_site_{i}.cif', config)

```

## Recommended Workflow for Your Study

### **For initial exploration**:

```python
# 1. Place on specific sites (La, V, O)
configs = place_on_surface_atoms(slab, 'NH3', atom_types=['La', 'V'])

# 2. Quick optimization (cheap calculator)
from ase.optimize import BFGS
# ... optimize each config

# 3. Select best 3-5 configurations
# 4. Do full DFT optimization

```

### **For thorough screening**:

```python
# Combine grid + orientations
configs = place_molecule_on_grid(slab, 'NH3',
                                  grid_density=5,
                                  orientations=4)
# This gives 5×5×4 = 100 configurations
# Pre-optimize with cheap method (e.g., DFTB, force field)
# Then DFT on top 10

```

## Post-Processing: Finding Unique Configurations

After generating many structures, you'll want to identify duplicates:

```python
from ase.io import read
from ase.geometry import get_distances
import numpy as np

def are_similar(atoms1, atoms2, threshold=0.5):
    """Check if two configurations are similar"""
    if len(atoms1) != len(atoms2):
        return False

    # Compare positions (accounting for periodic boundaries)
    pos1 = atoms1.get_positions()
    pos2 = atoms2.get_positions()

    # Simple RMSD check
    rmsd = np.sqrt(np.mean((pos1 - pos2)**2))
    return rmsd < threshold

def filter_unique_configs(filenames, threshold=0.5):
    """Remove duplicate configurations"""
    unique = []

    for fname in filenames:
        atoms = read(fname)
        is_unique = True

        for existing in unique:
            if are_similar(atoms, existing, threshold):
                is_unique = False
                break

        if is_unique:
            unique.append(atoms)

    print(f"Reduced {len(filenames)} configs to {len(unique)} unique")
    return unique

```

## My Recommendation

For **NH3 on LaVO3**:

1. **Start with Method 4** (site-specific placement):
    - Target La sites (Lewis acid for NH3)
    - Try 3-5 orientations per site
    - This is physically motivated
2. **Add a few random configurations** (Method 2):
    - Just to be safe you haven't missed something
    - 10-20 random configs
3. **Pre-screen with cheap method**:
    - Use EMT, DFTB, or single-point DFT
    - Rank by energy
4. **Full DFT on top candidates**:
    - Optimize top 5-10 configurations
    - These become your starting points for reaction pathway

Would you like me to write a complete script tailored to your LaVO3/NH3 system?

### **STEP 1.3 Beyond single configuration 0K energy**

Great question! Let's break down what you have, what you need for Gibbs free energy, and what conclusions you can draw from your current data.

## What You Currently Have

From your DFT calculations:

- **Adsorption energies (E_ads)** at T = 0 K for different configurations:
    - NH₃* (molecularly adsorbed)
    - NH₂* + H*
    - NH* + 2H*
    - N* + 3H*
    - etc.

These are **electronic energies** (ΔE₀) including zero-point energy (ZPE) at best.

## What's Missing for Gibbs Free Energy

To get **ΔG(T, P)**, you need to add **thermal corrections**:

### **ΔG = ΔE₀ + ΔZPE + ΔH(T) - TΔS(T)**

Where:

- **ΔE₀**: Electronic energy (you have this)
- **ΔZPE**: Zero-point energy correction
- **ΔH(T)**: Enthalpy correction (vibrational, rotational, translational)
- **TΔS(T)**: Entropy term (temperature-dependent)

### What You Need to Calculate:

### 1. **Vibrational Frequencies** (Phonon/Hessian calculations)

**For surface species** (NH₃*, NH₂*, NH*, N*, H*):

```
- Run frequency calculations on optimized geometries
- Fix substrate atoms (or use partial Hessian)
- Get vibrational modes of adsorbates only

```

**Extract**:

- ZPE = (1/2)Σ ℏωᵢ
- Vibrational entropy: S_vib(T)
- Vibrational enthalpy: H_vib(T)

### 2. **Gas Phase Reference States**

For **NH₃(g)** and **H₂(g)**:

```
- Full frequency calculation in a box
- Or use experimental/tabulated data

```

**Extract**:

- Translational entropy: S_trans (from Sackur-Tetrode)
- Rotational entropy: S_rot
- Vibrational contributions
- PV work term (≈ RT for ideal gas)

### 3. **Temperature and Pressure**

Specify reaction conditions:

- **T**: Typical NH₃ decomposition: 400-800 K
- **P**: Standard (1 bar) or operating pressure

## Minimal Data Needed

### **Option 1: Zero-Point Energy Only** (Quick approximation)

Calculate **only ZPE** for all species:

**ΔE_rxn → ΔE₀ = ΔE_rxn + ΔZPE**

This corrects for quantum effects but ignores entropy.

**Computational cost**: Moderate (frequency calculations)

### **Option 2: Full Gibbs Free Energy** (Complete thermodynamics)

Calculate frequencies + apply statistical mechanics:

**ΔG(T) = ΔE₀ + ΔG_thermal**

Where ΔG_thermal includes all thermal corrections.

**Computational cost**: Same as Option 1, but more post-processing

## Practical Workflow

### Using VASP (example):

```bash
# 1. Optimize geometry
IBRION = 2
NSW = 100

# 2. Calculate frequencies
IBRION = 5  # or 6 or 7
NFREE = 2
POTIM = 0.015
# Selectively displace only adsorbate atoms (SELECTIVE DYNAMICS)

```

### Using scripts (ASE/pymatgen):

```python
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from ase.vibrations import Vibrations

# For adsorbed species
vib = Vibrations(atoms, indices=[...])  # indices of adsorbate
vib.run()
vib_energies = vib.get_energies()

# Calculate thermodynamic corrections
thermo = HarmonicThermo(vib_energies=vib_energies)
ZPE = thermo.get_ZPE_correction()
G_correction = thermo.get_helmholtz_energy(T, verbose=False)

# For gas phase molecules
thermo_gas = IdealGasThermo(
    vib_energies=vib_energies,
    geometry='nonlinear',  # for NH3
    atoms=molecule,
    symmetrynumber=3,  # for NH3
    spin=0
)
G_gas = thermo_gas.get_gibbs_energy(T, P)

```

## What You CAN Conclude Without Gibbs Free Energy

### ✅ **Valid Conclusions from ΔE_ads alone**:

### 1. **Relative Stability at T ≈ 0 K**

- Which configuration is most stable?
- **Example**: If E_ads(NH₃*) > E_ads(NH₂* + H*), then NH₃ prefers to dissociate at 0 K

### 2. **Reaction Energy (ΔE_rxn)**

Calculate electronic energy changes:

```
NH₃* → NH₂* + H*          ΔE₁ = E(NH₂* + H*) - E(NH₃*)
NH₂* → NH* + H*           ΔE₂ = E(NH* + H*) - E(NH₂*)
NH* → N* + H*             ΔE₃ = E(N* + H*) - E(NH*)

```

**Interpretation**:

- ΔE < 0: Exothermic (thermodynamically favorable at 0 K)
- ΔE > 0: Endothermic (requires energy input)

### 3. **Binding Strength Comparison**

- Which intermediate binds most strongly?
- **Example**: If |E_ads(N*)| >> |E_ads(NH₃*)|, nitrogen may poison the surface

### 4. **Qualitative Reaction Pathway**

Sketch an energy profile:

```
       NH₂*+H*
         /\
        /  \
NH₃*   /    \   NH*+2H*
  |   /      \    /\
  |  /        \  /  \
  | /          \/    \
  |/                  N*+3H*

```

**Identify**:

- Most stable intermediate
- Potential rate-determining step (largest uphill step)

### 5. **Structure-Reactivity Relationships**

- Compare different surface terminations (LaO vs VO₂)
- Compare different adsorption sites
- **Example**: "NH₃ binds 0.5 eV stronger on LaO than VO₂"

### ⚠️ **What You CANNOT Conclude**:

❌ **Actual operating temperature behavior**

- Entropy becomes dominant at high T
- Stable at 0 K ≠ stable at 700 K

❌ **Reaction spontaneity** (ΔG > 0 or < 0)

- Need entropy to determine this

❌ **Equilibrium constants**

- K_eq = exp(-ΔG/RT)

❌ **Coverage effects at realistic conditions**

- Requires chemical potential (μ = G)

❌ **Gas-surface equilibrium**

- NH₃(g) ⇌ NH₃* requires ΔG_ads(T,P)

## Practical Recommendations

### **If you have limited time/resources**:

**Minimum viable analysis**:

1. ✅ **Add ZPE corrections** (relatively cheap)
    - Changes binding energies by ~0.1-0.3 eV typically
    - Important for N-H bond breaking
2. ✅ **Use ΔE₀ to identify trends**
    - Most stable intermediate
    - Most endothermic step
    - Compare surfaces/sites
3. ✅ **Discuss limitations**
    - "These are 0 K energies; entropy effects at operating T may change rankings"

### **For a complete study**:

1. ✅ Calculate vibrational frequencies for all species
2. ✅ Apply thermal corrections at relevant T (e.g., 600 K)
3. ✅ Construct full ΔG reaction profile
4. ✅ Calculate transition states (if needed) for activation barriers
5. ✅ Compare ΔG values to determine rate-determining step

## Example Interpretation

### **With only ΔE_ads**:

> "NH₃ dissociation to NH₂* + H* is endothermic by 0.8 eV at the electronic level, suggesting this step may be challenging. The N* species is very strongly bound (-4.2 eV), which could lead to surface poisoning."
> 

### **With ΔG(T)**:

> "At 600 K, NH₃ dissociation has ΔG = +0.3 eV, making it slightly unfavorable. However, H₂ desorption (ΔG = -0.5 eV) drives the overall reaction forward. The rate-determining step is N* + H* → NH* with ΔG‡ = 1.2 eV."
> 

## Bottom Line

- **You can publish** meaningful results with just ΔE_ads (many papers do), but acknowledge limitations
- **For NH₃ decomposition** (high T reaction), Gibbs free energy is highly recommended
- **Minimum addition**: Calculate ZPE (adds ~1 day of computation)
- **Complete analysis**: Full vibrational analysis + thermal corrections (adds ~1 week)

What's your timeline and computational budget? I can help you decide which level of theory is appropriate for your study.

## **STEP 2: perform the calculation**

*Q. When I calculate the decomposition of NH3, should I put NH3 on the surface, then NH2, and H; the NH, H, H on the surface?*

## The Key Question: Sequential vs Co-adsorbed

### **Option 1: Co-adsorbed species (RECOMMENDED)**

Put all fragments on the surface simultaneously:

- **NH3**→ *NH2 + H** → *NH + 2H** → *N + 3H**

### **Option 2: Sequential removal**

Put only the current fragment:

- **NH3**→ **NH2**→ **NH**→ **N**(H atoms removed to gas phase as H2)

## Understanding the Mechanisms

### **Mechanism A: H remains on surface** (Co-adsorption)

`Step 1: NH3* → NH2* + H*
Step 2: NH2* + H* → NH* + 2H*
Step 3: NH* + 2H* → N* + 3H*
Step 4: 2H* → H2(g)  [separate step]`

**Physical picture**:

- H atoms accumulate on surface during decomposition
- H-H recombination and desorption happens separately
- This is the **more realistic mechanism** for catalytic conditions

**When to use**:

- ✅ Studying the decomposition pathway itself
- ✅ Realistic catalytic conditions (H coverage matters)
- ✅ Understanding surface chemistry
- ✅ Finding rate-determining step

### **Step-by-step what to calculate c**o-adsorption:

### **1. Initial state: NH3 adsorbed**

`Configuration: NH3* on clean surface
Calculate: E(NH3*/slab)`

### **2. First dehydrogenation: NH2 + H**

`Configuration: NH2* and H* both on surface
Calculate: E(NH2* + H*/slab)

Important: Try multiple configurations!
- H near NH2 (associated)
- H far from NH2 (dissociated)
- Different H adsorption sites`

### **3. Second dehydrogenation: NH + 2H**

`Configuration: NH* and 2H* on surface
Calculate: E(NH* + 2H*/slab)

Important: 
- Try different H-H distances
- Different H sites
- NH position variations`

### **4. Third dehydrogenation: N + 3H**

`Configuration: N* and 3H* on surface
Calculate: E(N* + 3H*/slab)

Important:
- 3H atoms can be spread out or clustered
- Check if H-H pairing occurs (toward H2 formation)`

### **5. H2 formation and desorption** (separate step)

`Configuration 1: 2H* → H2* (molecular on surface)
Configuration 2: H2* → H2(g) (desorption)`

## Practical Workflow

### **Step 1: Find stable adsorption sites**

`from ase.io import read, write
from ase.build import molecule

slab = read('LaVO3_slab.cif')

# Place NH3 on La site (Lewis acid)
configs_NH3 = place_on_surface_atoms(slab, 'NH3', atom_types=['La'])

# Optimize all, find lowest energy
# This is your reference: NH3*`

### **Step 2: Generate NH2 + H configurations**

python

`# Start from optimized NH3* geometry
nh3_optimized = read('NH3_optimized.traj')

# Remove one H from NH3
# Method 1: Manual
# - Delete H atom from NH3
# - This leaves NH2 in same position

# Method 2: Using ASE
from ase import Atoms
import numpy as np

def create_NH2_H_configs(nh3_on_slab, n_configs=10):
    """
    Create multiple NH2 + H configurations
    """
    configs = []
    
    # Find NH3 atoms (last 4 atoms if you added NH3 to slab)
    n_slab = len(nh3_on_slab) - 4
    nh3_indices = list(range(n_slab, len(nh3_on_slab)))
    
    # Identify N and H atoms
    n_idx = nh3_indices[0]  # Assuming N is first
    h_indices = nh3_indices[1:]  # Three H atoms` 

`*# Remove one H, place it on different surface sites*`

`for i, h_to_remove in enumerate(h_indices):
        for config_num in range(n_configs // 3):
            
            new_slab = nh3_on_slab.copy()`

`*# Get H position*`

`h_pos = new_slab[h_to_remove].position.copy()`

`*# Find surface sites for H*`

`*# Option 1: Random displacement*`

`x_disp = np.random.uniform(-3, 3)
            y_disp = np.random.uniform(-3, 3)
            new_h_pos = h_pos + [x_disp, y_disp, 0]`

`*# Option 2: Place near specific surface atoms*`

`*# (La, V, O sites)*`

`*# Update H position*`

`new_slab[h_to_remove].position = new_h_pos
            
            configs.append(new_slab)
            write(f'NH2_H_config_{len(configs)}.cif', new_slab)
    
    return configs`

`*# Generate configurations*`

`configs = create_NH2_H_configs(nh3_optimized, n_configs=15)`

### **Step 3: Systematic H placement strategy**

For **NH2 + H**, **NH + 2H**, **N + 3H**, I recommend:

python

`def place_H_atoms_on_sites(slab_with_NHx, n_H_atoms, site_types=['La', 'V', 'O']):
    """
    Place H atoms on different surface sites
    
    Strategy:
    - Try each H on different types of sites
    - Try different H-H separations
    - Try H near vs far from NHx fragment
    """
    
    # Find surface sites
    surface_sites = get_surface_sites(slab_with_NHx, site_types)
    
    configs = []
    
    if n_H_atoms == 1:
        # Simple: place H on each site type
        for site in surface_sites[:10]:  # Top 10 sites
            config = slab_with_NHx.copy()
            # Add H at site
            # ... 
            configs.append(config)
    
    elif n_H_atoms == 2:
        # More complex: try different H-H distances
        for site1 in surface_sites[:5]:
            for site2 in surface_sites[:5]:`

`if site1 != site2:
                    config = slab_with_NHx.copy()`

`*# Add 2H atoms*`

`*# ...*`

`configs.append(config)
    
    elif n_H_atoms == 3:`

`*# Even more complex*`

`*# Try: (1) all H clustered, (2) spread out, (3) pairs*`

`*# ...*`

`pass
    
    return configs
````

`*## What Energies to Extract*`

`For each configuration, calculate:`

`*### **Adsorption Energy** (relative to gas phase):*`

````
E_ads(NH3*) = E(NH3*/slab) - E(slab) - E(NH3_gas)

E_ads(NH2* + H*) = E(NH2* + H*/slab) - E(slab) - E(NH3_gas)
                 = E(NH2* + H*/slab) - E(slab) - [E(NH2_gas) + 1/2 E(H2_gas)]
````

`*### **Reaction Energy** (step-by-step):*`

````
ΔE1 = E(NH2* + H*/slab) - E(NH3*/slab)
ΔE2 = E(NH* + 2H*/slab) - E(NH2* + H*/slab)
ΔE3 = E(N* + 3H*/slab) - E(NH* + 2H*/slab)`

`````

`*### **H2 desorption energy**:*`

````
ΔE_des = E(slab + H2_gas) - E(2H*/slab)
       = E(slab) + E(H2_gas) - E(2H*/slab)
````

`*## Coverage Effects*`

`**Important consideration**: How many H atoms are realistic?`

`*### **Low H coverage** (what you're calculating):*`

`- 1 NH3, produces up to 3 H atoms
- Coverage: ~0.1-0.25 ML (monolayer)`

`*### **High H coverage** (might be more realistic):*`

`- Multiple NH3 molecules decomposing
- H coverage: 0.5-1.0 ML
- **This changes energetics!**

**Recommendation**:
1. Start with **single NH3** (low coverage) for mechanism
2. Later, check effect of **H pre-coverage** if you have time
   - Add extra H atoms to surface
   - Recalculate NH3 → NH2 + H with H already present`

`*## Example Energy Diagram*`

`````

`Energy
^
|
|     TS2         TS3
|      /\          /\
|     /  \        /  \`

`|    /    \      /    \
| TS1      NH*+2H*     N*+3H*
|  /\      /            /
| /  \    /            /
|/    NH2*+H*         /
|      /             /
NH3*   /             /
|   /             /
|  /    H2 desorption
| /              |
|/               v
+-------------------> Reaction coordinate`

`**What you calculate**:
- Points: NH3*, NH2*+H*, NH*+2H*, N*+3H*
- Transition states (TS): If you do NEB or CI-NEB later`

`*## Common Mistakes to Avoid*`

`❌ **Mistake 1**: Removing H to gas phase immediately
- Wrong: NH3* → NH2* + 1/2 H2(g)
- Right: NH3* → NH2* + H*

❌ **Mistake 2**: Using only one configuration per step
- You MUST try multiple H positions
- Lowest energy configuration might not be obvious

❌ **Mistake 3**: Ignoring H-H interactions
- At NH*+2H* and N*+3H*, H atoms can cluster
- This affects energetics

❌ **Mistake 4**: **Not relaxing the entire system**
- Always relax slab + adsorbates together
- Fix bottom layers only`

`*## Summary: What to Calculate*`

`*### **Minimum viable study**:*`

````
1. NH3* (1 config - lowest energy)
2. NH2* + H* (5-10 configs - different H sites)
3. NH* + 2H* (5-10 configs - different H arrangements)
4. N* + 3H* (3-5 configs - different H arrangements)
5. 2H* → H2* → H2(g) (2-3 configs)`

### **Comprehensive study**:

- All of the above
- Plus: Transition states (NEB/CI-NEB)
- Plus: Coverage dependence
- Plus: Different surface terminations (LaO vs VO2)

## My Recommendation for You

**Phase 1**: Co-adsorption (NH2+H, NH+2H, N+3H)

- This gives you the decomposition pathway
- Most physically meaningful

**Phase 2** (if needed): H2 formation/desorption

- Understand if H removal is limiting

**Phase 3** (optional): Compare to sequential mechanism

- Shows importance of H co-adsorption

STEP 3: Understanding results

Excellent question! Comparing surfaces for catalytic performance requires analyzing multiple descriptors. Let me show you how to determine which surface is "better" for NH3 decomposition.

## Key Metrics to Compare Surfaces

### 1. **Thermodynamic Criteria** (Feasibility)

### A. **Overall Reaction Energy**

For the complete decomposition: *NH3 → N* + 3H***

```python
# For each surface
ΔE_total = E(N* + 3H*/slab) - E(NH3*/slab)

```

**Interpretation**:

- More negative = more exothermic = thermodynamically favorable
- But **too negative** can mean N* poisoning (see below)

### B. **Step-by-step Reaction Energies**

```python
ΔE1 = E(NH2* + H*/slab) - E(NH3*/slab)       # First N-H break
ΔE2 = E(NH* + 2H*/slab) - E(NH2* + H*/slab)  # Second N-H break
ΔE3 = E(N* + 3H*/slab) - E(NH* + 2H*/slab)   # Third N-H break

```

**What to look for**:

- ✅ **All steps should be feasible** (not too endothermic)
- ✅ **Balanced energetics** (no single step >>1 eV endothermic)
- ⚠️ Very endothermic step (>1.5 eV) = likely bottleneck

### C. **Hydrogen Desorption Energy**

```python
ΔE_H2_des = E(slab) + E(H2_gas) - E(2H*/slab)

```

**Critical for catalysis**:

- ✅ **Should be exothermic or slightly endothermic** (ΔE ≈ 0 to +0.5 eV)
- ❌ Too negative (<-1 eV): H binds too weakly, may not dissociate NH3
- ❌ Too positive (>+1.5 eV): H sticks too strongly, poisons surface

### 2. **Kinetic Criteria** (Activity)

### A. **Rate-Determining Step (RDS)**

The step with the **highest activation barrier** or **most endothermic**:

```python
# Without transition states (approximation)
steps = {
    'NH3 → NH2+H': ΔE1,
    'NH2+H → NH+2H': ΔE2,
    'NH+2H → N+3H': ΔE3,
    'H2 desorption': ΔE_H2_des
}

RDS = max(steps.items(), key=lambda x: x[1])

```

**Better surface**:

- ✅ Lower RDS barrier
- ✅ RDS is a different step than other surfaces (complementary)

### B. **Activation Barriers** (if you calculate with NEB)

```python
Ea1 = E(TS_NH3→NH2) - E(NH3*)
Ea2 = E(TS_NH2→NH) - E(NH2* + H*)
Ea3 = E(TS_NH→N) - E(NH* + 2H*)

```

**Better surface**:

- ✅ Lower activation barriers overall
- ✅ Particularly for the RDS

### 3. **Electronic Structure Criteria**

### A. **Binding Strength Hierarchy**

Calculate adsorption energies:

```python
E_ads(NH3) = E(NH3*/slab) - E(slab) - E(NH3_gas)
E_ads(NH2) = E(NH2*/slab) - E(slab) - E(NH2_gas)
E_ads(NH) = E(NH*/slab) - E(slab) - E(NH_gas)
E_ads(N) = E(N*/slab) - E(slab) - E(N_gas)
E_ads(H) = E(H*/slab) - E(slab) - 0.5*E(H2_gas)

```

**Ideal trends** (Sabatier principle):

- ✅ **Moderate binding**: Not too weak, not too strong
- ✅ **Order**: |E_ads(N)| > |E_ads(NH)| > |E_ads(NH2)| > |E_ads(NH3)|
- ⚠️ **N poisoning check**: If |E_ads(N)| >> 4 eV, N may not desorb

### B. **Scaling Relations**

Plot binding energies against each other:

```python
import matplotlib.pyplot as plt

# For all three surfaces
surfaces = ['LaO-term', 'VO2-term', 'Surface3']
E_N = [-4.2, -3.8, -3.5]  # Example values
E_H = [-2.8, -2.5, -2.2]

plt.scatter(E_N, E_H)
for i, surf in enumerate(surfaces):
    plt.annotate(surf, (E_N[i], E_H[i]))

plt.xlabel('E_ads(N) [eV]')
plt.ylabel('E_ads(H) [eV]')
plt.title('Scaling Relation')
plt.show()

```

**Look for**:

- Surface that breaks scaling relations (outlier can be better)
- Balance between N and H binding

### 4. **Catalyst Descriptor Approach**

### **Volcano Plot** (if you have enough data)

```python
import numpy as np
import matplotlib.pyplot as plt

# Theoretical: Activity vs. binding energy descriptor
# Typically use E_ads(N) or E_ads(NH) as descriptor

def activity_model(E_descriptor):
    """
    Simple BEP-type model
    Activity ∝ exp(-ΔG_apparent/kT)
    """
    # Simplified: assume linear scaling
    if E_descriptor < -4.0:  # Too strong binding
        return np.exp(-abs(E_descriptor + 4.0))
    elif E_descriptor > -2.5:  # Too weak binding
        return np.exp(-abs(E_descriptor + 2.5))
    else:  # Optimal range
        return 1.0

# Your surfaces
surfaces = {
    'LaO-term': -3.8,
    'VO2-term': -3.5,
    'Surface3': -4.2
}

activities = {surf: activity_model(E) for surf, E in surfaces.items()}

# Plot
E_range = np.linspace(-5, -2, 100)
activity_range = [activity_model(E) for E in E_range]

plt.plot(E_range, activity_range, 'k-', label='Volcano')
for surf, E in surfaces.items():
    plt.scatter(E, activities[surf], s=100, label=surf)

plt.xlabel('E_ads(N) [eV]')
plt.ylabel('Relative Activity')
plt.legend()
plt.show()

```

**Best surface**: Closest to volcano peak

### 5. **Surface Stability Criteria**

### A. **Surface Energy**

```python
γ = (E_slab - N_units * E_bulk) / (2 * Area)

```

**Better surface**:

- ✅ Lower surface energy (more stable, more abundant)
- Balance with activity (stable but inactive surface is useless)

### B. **Reconstruction Check**

After NH3 adsorption/decomposition:

- ❌ Significant surface reconstruction = unstable
- ✅ Minimal reconstruction = robust catalyst

### 6. **Poisoning Resistance**

### A. *N Coverage Effects*

If N binds too strongly:

```python
# Check if N can be removed
ΔE_N_removal = E(slab) + E(N_gas) - E(N*/slab)
# or
ΔE_N2_formation = E(slab) + E(N2_gas) - E(2N*/slab)

```

**Better surface**:

- ✅ N removal is feasible (not >2 eV endothermic)

### B. **H Coverage Effects**

Test stability with multiple H atoms:

```python
# Add extra H to surface, recalculate NH3 adsorption
E_ads(NH3, with_H_precov) vs E_ads(NH3, clean)

```

**Better surface**:

- ✅ Activity doesn't drop significantly with H coverage

## Practical Comparison Framework

### **Method 1: Energy Span Model**

```python
def energy_span(energies):
    """
    Calculate energy span (apparent barrier)
    Based on: Campbell, J. Catal. 2001
    """
    intermediates = list(energies.values())
    E_span = max(intermediates) - min(intermediates)
    return E_span

# For each surface
surface1_energies = {
    'NH3*': 0.0,
    'NH2*+H*': 0.8,
    'NH*+2H*': 1.5,
    'N*+3H*': 1.2
}

surface2_energies = {
    'NH3*': 0.0,
    'NH2*+H*': 0.6,
    'NH*+2H*': 1.2,
    'N*+3H*': 0.9
}

span1 = energy_span(surface1_energies)  # 1.5 eV
span2 = energy_span(surface2_energies)  # 1.2 eV

print(f"Surface 1 span: {span1:.2f} eV")
print(f"Surface 2 span: {span2:.2f} eV")
# Lower span = better

```

### **Method 2: Weighted Scoring**

```python
def score_surface(metrics):
    """
    Score surface based on multiple criteria
    """
    score = 0

    # Thermodynamics (30%)
    if -1.0 < metrics['ΔE_total'] < 0:
        score += 30

    # RDS barrier (40%)
    if metrics['RDS_barrier'] < 1.0:
        score += 40
    elif metrics['RDS_barrier'] < 1.5:
        score += 20

    # H2 desorption (20%)
    if -0.5 < metrics['ΔE_H2_des'] < 0.5:
        score += 20

    # N binding (10%)
    if -4.5 < metrics['E_ads_N'] < -3.0:
        score += 10

    return score

# Example
surface_A = {
    'ΔE_total': -0.5,
    'RDS_barrier': 1.2,
    'ΔE_H2_des': 0.2,
    'E_ads_N': -3.8
}

surface_B = {
    'ΔE_total': -0.3,
    'RDS_barrier': 0.9,
    'ΔE_H2_des': 0.4,
    'E_ads_N': -3.5
}

print(f"Surface A score: {score_surface(surface_A)}")
print(f"Surface B score: {score_surface(surface_B)}")

```

### **Method 3: Visual Comparison**

```python
import matplotlib.pyplot as plt
import numpy as np

def plot_energy_profiles(surfaces_data):
    """
    Plot energy profiles for all surfaces
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    steps = ['NH3*', 'NH2*+H*', 'NH*+2H*', 'N*+3H*']
    x = np.arange(len(steps))

    for surf_name, energies in surfaces_data.items():
        ax.plot(x, energies, 'o-', label=surf_name, linewidth=2, markersize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(steps)
    ax.set_ylabel('Relative Energy [eV]')
    ax.set_xlabel('Reaction Coordinate')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('surface_comparison.png', dpi=300)
    plt.show()

# Usage
surfaces_data = {
    'LaO-term': [0.0, 0.8, 1.5, 1.2],
    'VO2-term': [0.0, 0.6, 1.2, 0.9],
    'Surface3': [0.0, 0.9, 1.7, 1.5]
}

plot_energy_profiles(surfaces_data)

```

## Complete Analysis Script

```python
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read

class SurfaceComparator:
    def __init__(self, surface_name):
        self.name = surface_name
        self.energies = {}

    def add_energy(self, state, energy):
        """Add energy for a given state"""
        self.energies[state] = energy

    def calculate_reaction_energies(self):
        """Calculate step-wise reaction energies"""
        steps = {
            'NH3 → NH2+H': self.energies['NH2*+H*'] - self.energies['NH3*'],
            'NH2+H → NH+2H': self.energies['NH*+2H*'] - self.energies['NH2*+H*'],
            'NH+2H → N+3H': self.energies['N*+3H*'] - self.energies['NH*+2H*'],
        }
        return steps

    def get_RDS(self):
        """Identify rate-determining step"""
        steps = self.calculate_reaction_energies()
        rds = max(steps.items(), key=lambda x: x[1])
        return rds

    def energy_span(self):
        """Calculate energy span"""
        energies = list(self.energies.values())
        return max(energies) - min(energies)

    def summary(self):
        """Print summary"""
        print(f"\n{'='*50}")
        print(f"Surface: {self.name}")
        print(f"{'='*50}")

        steps = self.calculate_reaction_energies()
        for step, energy in steps.items():
            print(f"{step:20s}: ΔE = {energy:+.2f} eV")

        rds_step, rds_energy = self.get_RDS()
        print(f"\nRate-determining step: {rds_step}")
        print(f"RDS barrier: {rds_energy:.2f} eV")
        print(f"Energy span: {self.energy_span():.2f} eV")

# Usage
surfaces = []

# Surface 1: LaO-terminated
surf1 = SurfaceComparator('LaO-terminated (001)')
surf1.add_energy('NH3*', 0.0)
surf1.add_energy('NH2*+H*', 0.8)
surf1.add_energy('NH*+2H*', 1.5)
surf1.add_energy('N*+3H*', 1.2)
surfaces.append(surf1)

# Surface 2: VO2-terminated
surf2 = SurfaceComparator('VO2-terminated (001)')
surf2.add_energy('NH3*', 0.0)
surf2.add_energy('NH2*+H*', 0.6)
surf2.add_energy('NH*+2H*', 1.2)
surf2.add_energy('N*+3H*', 0.9)
surfaces.append(surf2)

# Surface 3: Another surface
surf3 = SurfaceComparator('(110) surface')
surf3.add_energy('NH3*', 0.0)
surf3.add_energy('NH2*+H*', 1.0)
surf3.add_energy('NH*+2H*', 1.8)
surf3.add_energy('N*+3H*', 1.6)
surfaces.append(surf3)

# Compare all surfaces
for surf in surfaces:
    surf.summary()

# Rank surfaces
print(f"\n{'='*50}")
print("RANKING (by energy span):")
print(f"{'='*50}")
ranked = sorted(surfaces, key=lambda s: s.energy_span())
for i, surf in enumerate(ranked, 1):
    print(f"{i}. {surf.name}: {surf.energy_span():.2f} eV")

```

## Decision Criteria Summary

**Best surface should have**:

1. ✅ **Low energy span** (<1.5 eV)
2. ✅ **Balanced step energies** (no single step >1.2 eV)
3. ✅ **Moderate N binding** (-4.0 to -3.0 eV)
4. ✅ **Favorable H2 desorption** (-0.5 to +0.5 eV)
5. ✅ **Low surface energy** (stable facet)
6. ✅ **No severe poisoning** (N or H)

**Red flags** ❌:

- Any step >1.5 eV endothermic
- |E_ads(N)| > 4.5 eV (poisoning)
- Energy span > 2.0 eV
- Significant surface reconstruction

## My Recommendation

For your **three surfaces**, create:

1. **Energy profile plot** (like above)
2. **Table of key metrics**:
    
    ```
    | Metric              | Surf 1 | Surf 2 | Surf 3 ||---------------------|--------|--------|--------|| Energy span         |        |        |        || RDS                 |        |        |        || E_ads(N)            |        |        |        || ΔE(H2 desorption)   |        |        |        |
    
    ```
    
3. **Discussion**: Which surface balances activity and stability?

Would you like help setting up this comparison with your actual data?

STEP 4: rate limiting step

Excellent question! Determining the **rate-determining step (RDS)** is crucial for understanding catalytic activity. Let me explain the different methods and their limitations.

## Methods to Identify the Rate-Determining Step

### **Method 1: Thermodynamic Approximation** (What You Can Do Now)

Using only **reaction energies** (no transition states):

### **A. Highest Endothermic Step**

```python
# Calculate reaction energies for each step
steps = {
    'NH3* → NH2*+H*': E_NH2_H - E_NH3,
    'NH2*+H* → NH*+2H*': E_NH_2H - E_NH2_H,
    'NH*+2H* → N*+3H*': E_N_3H - E_NH_2H,
    'H2 desorption': E_clean + E_H2_gas - E_2H
}

# Find most endothermic
rds_approx = max(steps.items(), key=lambda x: x[1])
print(f"Approximate RDS: {rds_approx[0]}")
print(f"ΔE = {rds_approx[1]:.2f} eV")

```

**Assumption**:

- Step with highest ΔE is likely the slowest
- Based on **Brønsted-Evans-Polanyi (BEP) relation**: ΔE‡ ∝ ΔE_rxn

**Limitations**:

- ❌ Ignores entropic effects
- ❌ Assumes similar pre-factors for all steps
- ❌ No actual activation barriers
- ⚠️ Only **qualitative guidance**

### **B. Energy Span Model**

More sophisticated thermodynamic approach:

```python
import numpy as np

def find_rds_energy_span(energies, T=673):
    """
    Energy span model (Campbell, 2001)
    Identifies apparent activation barrier

    Parameters:
    -----------
    energies : dict
        {state: relative_energy} in eV
    T : float
        Temperature in K
    """
    kB = 8.617e-5  # eV/K

    states = list(energies.keys())
    E = np.array(list(energies.values()))

    # Find TOF-determining transition and intermediate
    # Simplified: look for max span
    E_max = E.max()
    E_min = E.min()

    energy_span = E_max - E_min

    # Identify states
    TDI = states[np.argmax(E)]  # TOF-determining intermediate
    TDTS = states[np.argmin(E)]  # TOF-determining transition state

    return {
        'energy_span': energy_span,
        'TDI': TDI,
        'TDTS': TDTS,
        'apparent_barrier': energy_span
    }

# Usage
energies = {
    'NH3*': 0.0,
    'NH2*+H*': 0.8,
    'NH*+2H*': 1.5,  # Highest
    'N*+3H*': 1.2
}

result = find_rds_energy_span(energies, T=673)
print(f"Energy span: {result['energy_span']:.2f} eV")
print(f"Rate-limiting: between {result['TDTS']} and {result['TDI']}")

```

**Better than simple ΔE, but still approximate**

---

### **Method 2: Calculate Transition States** (ACCURATE, but expensive)

Use **NEB (Nudged Elastic Band)** or **CI-NEB** to find activation barriers:

### **A. NEB Calculation Setup**

```python
from ase.neb import NEB
from ase.optimize import BFGS
from ase.io import read, write

# For each reaction step: NH3* → NH2*+H*

# 1. Initial state (reactant)
initial = read('NH3_on_surface.traj')

# 2. Final state (product)
final = read('NH2_H_on_surface.traj')

# 3. Create intermediate images
images = [initial]
images += [initial.copy() for i in range(7)]  # 7 intermediate images
images += [final]

# 4. Interpolate
neb = NEB(images)
neb.interpolate()

# 5. Attach calculator to intermediate images only
for image in images[1:-1]:
    image.calc = your_calculator  # MLIP or DFT

# 6. Optimize the band
optimizer = BFGS(neb, trajectory='neb.traj')
optimizer.run(fmax=0.05)

# 7. Get activation barrier
energies = [image.get_potential_energy() for image in images]
Ea_forward = max(energies) - energies[0]
Ea_reverse = max(energies) - energies[-1]

print(f"Forward barrier: {Ea_forward:.2f} eV")
print(f"Reverse barrier: {Ea_reverse:.2f} eV")

```

### **B. Complete Pathway Analysis**

```python
def analyze_full_pathway():
    """
    Calculate activation barriers for all steps
    """
    steps = [
        ('NH3* → NH2*+H*', 'nh3_nh2h'),
        ('NH2*+H* → NH*+2H*', 'nh2h_nh2h'),
        ('NH*+2H* → N*+3H*', 'nh2h_n3h'),
    ]

    barriers = {}

    for step_name, file_prefix in steps:
        # Run NEB
        Ea = run_neb(f'{file_prefix}_initial.traj',
                     f'{file_prefix}_final.traj')
        barriers[step_name] = Ea

    # Find RDS
    rds = max(barriers.items(), key=lambda x: x[1])

    return barriers, rds

barriers, rds = analyze_full_pathway()
print(f"\nActivation barriers:")
for step, Ea in barriers.items():
    marker = " ← RDS" if step == rds[0] else ""
    print(f"{step}: {Ea:.2f} eV{marker}")

```

**Advantages**:

- ✅ Gives actual activation barriers
- ✅ Most accurate method
- ✅ Gold standard for publications

**Disadvantages**:

- ❌ Very expensive (7-15 optimizations per step)
- ❌ Requires good initial guess for TS
- ❌ Can fail for complex mechanisms

---

### **Method 3: Free Energy Analysis** (COMPLETE, includes entropy)

Calculate **Gibbs free energy** barriers at reaction temperature:

```python
from ase.thermochemistry import HarmonicThermo
from ase.vibrations import Vibrations

def calculate_free_energy_barrier(initial_state, final_state, ts_state, T=673):
    """
    Calculate ΔG‡ including entropy

    Steps:
    1. Calculate frequencies for initial, TS, final states
    2. Get thermal corrections
    3. Calculate ΔG‡ = G(TS) - G(initial)
    """

    states = {
        'initial': initial_state,
        'TS': ts_state,
        'final': final_state
    }

    G = {}

    for name, atoms in states.items():
        # 1. Electronic energy
        atoms.calc = your_calculator
        E = atoms.get_potential_energy()

        # 2. Vibrational analysis (only adsorbates)
        vib = Vibrations(atoms, indices=adsorbate_indices)
        vib.run()
        vib_energies = vib.get_energies()

        # Remove imaginary modes for TS (keep only real modes for G calculation)
        if name == 'TS':
            vib_energies = [e for e in vib_energies if e > 0]

        # 3. Thermochemistry
        thermo = HarmonicThermo(vib_energies=vib_energies,
                                electronicenergy=E)

        # Get Helmholtz free energy (for adsorbates)
        G[name] = thermo.get_helmholtz_energy(T)

    # Calculate barrier
    dG_forward = G['TS'] - G['initial']
    dG_reverse = G['TS'] - G['final']

    return dG_forward, dG_reverse

# For each step
dG_barrier_step1 = calculate_free_energy_barrier(nh3_state, nh2h_state, ts1_state, T=673)
# ... repeat for all steps

# Find RDS
rds_step = max(all_barriers, key=lambda x: x[1])

```

**Most complete analysis**:

- ✅ Includes entropy (ΔS)
- ✅ Temperature-dependent
- ✅ Most realistic

**Cost**:

- Very expensive (frequencies + NEB for all steps)

---

### **Method 4: Microkinetic Modeling** (MOST RIGOROUS)

Solve rate equations for the entire mechanism:

```python
import numpy as np
from scipy.integrate import odeint

def microkinetic_model(energies, barriers, T=673, P_NH3=1.0):
    """
    Solve microkinetic model to determine RDS

    Returns:
    --------
    rates : dict
        Reaction rate for each step
    coverages : dict
        Surface coverage of each species
    """

    kB = 8.617e-5  # eV/K
    h = 4.136e-15  # eV·s

    # Calculate rate constants (Eyring equation)
    def rate_constant(Ea, dG_rxn):
        k_forward = (kB * T / h) * np.exp(-Ea / (kB * T))
        k_reverse = (kB * T / h) * np.exp(-(Ea - dG_rxn) / (kB * T))
        return k_forward, k_reverse

    # For each step
    k = {}
    k['ads'] = rate_constant(0, energies['NH3*'] - 0)  # Adsorption
    k['step1'] = rate_constant(barriers['step1'], energies['NH2*+H*'] - energies['NH3*'])
    k['step2'] = rate_constant(barriers['step2'], energies['NH*+2H*'] - energies['NH2*+H*'])
    k['step3'] = rate_constant(barriers['step3'], energies['N*+3H*'] - energies['NH*+2H*'])
    k['des'] = rate_constant(barriers['H2_des'], 0)  # Desorption

    # Solve rate equations (steady-state or dynamic)
    def rate_equations(coverages, t):
        theta_NH3, theta_NH2, theta_NH, theta_N, theta_H = coverages
        theta_free = 1 - sum(coverages)

        d_theta = np.zeros(5)

        # NH3 adsorption/dissociation
        d_theta[0] = k['ads'][0] * P_NH3 * theta_free - k['step1'][0] * theta_NH3

        # NH2+H formation/consumption
        d_theta[1] = k['step1'][0] * theta_NH3 - k['step2'][0] * theta_NH2 * theta_H

        # NH+2H formation/consumption
        d_theta[2] = k['step2'][0] * theta_NH2 * theta_H - k['step3'][0] * theta_NH * theta_H**2

        # N+3H formation
        d_theta[3] = k['step3'][0] * theta_NH * theta_H**2

        # H balance
        d_theta[4] = (k['step1'][0] * theta_NH3 + k['step2'][0] * theta_NH2 * theta_H
                     + k['step3'][0] * theta_NH * theta_H**2
                     - 2 * k['des'][0] * theta_H**2)

        return d_theta

    # Solve
    coverages_init = [0.1, 0, 0, 0, 0]
    t = np.linspace(0, 1000, 10000)
    solution = odeint(rate_equations, coverages_init, t)

    # Steady-state coverages
    coverages_ss = solution[-1]

    # Calculate step rates
    rates = {
        'step1': k['step1'][0] * coverages_ss[0],
        'step2': k['step2'][0] * coverages_ss[1] * coverages_ss[4],
        'step3': k['step3'][0] * coverages_ss[2] * coverages_ss[4]**2,
    }

    # RDS is the slowest step
    rds = min(rates.items(), key=lambda x: x[1])

    return rates, coverages_ss, rds

# Run microkinetic model
rates, coverages, rds = microkinetic_model(energies, barriers, T=673)

print(f"\nSteady-state rates:")
for step, rate in rates.items():
    print(f"{step}: {rate:.2e} s⁻¹")
print(f"\nRDS: {rds[0]} (rate = {rds[1]:.2e} s⁻¹)")

```

**Most powerful**:

- ✅ Accounts for all steps simultaneously
- ✅ Includes coverage effects
- ✅ Pressure/temperature dependence
- ✅ Can predict turnover frequency (TOF)

**Challenges**:

- Requires activation barriers for all steps
- Requires pre-exponential factors
- More complex to set up

---

## Practical Recommendations

### **Your Current Situation** (only have energies):

### **Quick Analysis** (5 minutes):

```python
# Simple approach
energies = {
    'NH3*': 0.0,
    'NH2*+H*': 0.8,
    'NH*+2H*': 1.5,
    'N*+3H*': 1.2,
    '2H*': -0.5,  # If you have H2 on surface
}

# Calculate step energies
steps = {
    'Step 1: NH3* → NH2*+H*': 0.8 - 0.0,
    'Step 2: NH2*+H* → NH*+2H*': 1.5 - 0.8,
    'Step 3: NH*+2H* → N*+3H*': 1.2 - 1.5,
}

print("Reaction energies:")
for step, dE in steps.items():
    print(f"{step}: ΔE = {dE:+.2f} eV")

# Approximate RDS
rds_approx = max(steps.items(), key=lambda x: x[1])
print(f"\nApproximate RDS: {rds_approx[0]}")
print(f"(Most endothermic step: ΔE = {rds_approx[1]:.2f} eV)")

print("\n⚠️ Warning: This is only an approximation!")
print("For accurate RDS, calculate transition states with NEB.")

```

### **Energy Profile Plot**:

```python
import matplotlib.pyplot as plt
import numpy as np

def plot_energy_profile(energies):
    """
    Visualize reaction coordinate
    """
    states = list(energies.keys())
    E = list(energies.values())

    x = np.arange(len(states))

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot energy profile
    ax.plot(x, E, 'o-', linewidth=2, markersize=10)

    # Highlight highest point
    max_idx = np.argmax(E)
    ax.plot(x[max_idx], E[max_idx], 'ro', markersize=15,
            label=f'Highest energy: {states[max_idx]}')

    # Labels
    ax.set_xticks(x)
    ax.set_xticklabels(states, rotation=45, ha='right')
    ax.set_ylabel('Relative Energy (eV)', fontsize=12)
    ax.set_xlabel('Reaction Coordinate', fontsize=12)
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.legend()
    ax.grid(alpha=0.3)

    # Annotate energy differences
    for i in range(len(x)-1):
        dE = E[i+1] - E[i]
        mid_x = (x[i] + x[i+1]) / 2
        mid_y = (E[i] + E[i+1]) / 2
        ax.annotate(f'ΔE={dE:+.2f}', xy=(mid_x, mid_y),
                   fontsize=9, ha='center',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig('energy_profile.png', dpi=300)
    plt.show()

# Usage
energies = {
    'NH₃*': 0.0,
    'NH₂*+H*': 0.8,
    'NH*+2H*': 1.5,
    'N*+3H*': 1.2
}

plot_energy_profile(energies)

```

---

### **For Publication-Quality Work**:

**Minimum requirement**:

1. ✅ Report reaction energies (ΔE)
2. ✅ Identify approximate RDS (highest ΔE step)
3. ⚠️ **Caveat**: "Approximate, based on thermodynamics"
4. ✅ State: "Activation barriers not calculated"

**Better**:

1. ✅ Calculate NEB for suspected RDS step only
2. ✅ Report one activation barrier
3. ✅ Validates thermodynamic approximation

**Best**:

1. ✅ NEB for all steps
2. ✅ Vibrational analysis
3. ✅ Free energy barriers at T
4. ✅ (Optional) Microkinetic modeling

---

## Decision Tree

```
Do you have computational resources for NEB?
│
├─ NO → Use thermodynamic approximation
│        ├─ Find most endothermic step
│        ├─ Report with caveat
│        └─ Suggest future work: TS calculations
│
└─ YES → Calculate transition states
         │
         ├─ Limited time?
         │  └─ Calculate TS for suspected RDS only
         │
         └─ Full study?
            ├─ NEB for all steps
            ├─ Add vibrational corrections
            └─ Microkinetic model

```

---

## Summary Table

| Method | Accuracy | Cost | When to Use |
| --- | --- | --- | --- |
| **Highest ΔE** | Low | Free | Initial screening |
| **Energy span** | Medium | Free | Better than ΔE alone |
| **NEB (one step)** | High | Moderate | Validate RDS |
| **NEB (all steps)** | High | High | Complete study |
| **Free energy (ΔG‡)** | Very High | Very High | Accurate kinetics |
| **Microkinetic** | Highest | Highest | TOF prediction |

---

## What I Recommend for You

**Right now** (with only energies):

```python
1. Calculate all step energies (ΔE)
2. Identify highest ΔE → "approximate RDS"
3. Plot energy profile
4. Report: "Based on reaction energetics, Step X appears rate-limiting"
5. Add: "Transition state calculations needed for confirmation"

```

**Next step** (if possible):

```python
1. Run NEB on the suspected RDS
2. Get actual activation barrier
3. Validates your thermodynamic analysis

```

Would you like help setting up NEB calculations, or should we work with the thermodynamic approximation for now?