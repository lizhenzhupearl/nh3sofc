# Plan: Surface Builder Enhancements for Perovskite Materials

**Date:** 2026-01-22
**Status:** Completed

## Summary

Enhance `SurfaceBuilder` and `SlabStructure` classes to properly handle perovskite surfaces (LaVO3, etc.) including:
- Polarity analysis and compensation
- Named termination selection (LaO vs VO2)
- Symmetric slab creation
- Stoichiometry validation
- Surface dipole calculation

## Current State

### What's Already Implemented
- Basic surface creation from bulk (`create_surface()`)
- Supercell expansion (`repeat_xy()`, `create_surface_with_size()`)
- All terminations enumeration (`get_all_terminations()`)
- Bottom layer fixing (`fix_bottom_layers()`)
- Surface atom detection (`get_surface_atoms()`)
- Vacuum management (`add_vacuum()`, `get_vacuum_size()`)

### What's Missing
1. Polarity check and compensation
2. Symmetric slab creation
3. Perovskite-specific termination naming (LaO, VO2)
4. Surface dipole moment calculation
5. Stoichiometry validation
6. Surface reconstruction support

---

## Phase 1: Polarity Analysis and Compensation

### 1.1 Add `calculate_layer_charges()` method to SlabStructure

Calculate formal charges per layer to detect polarity.

```python
def calculate_layer_charges(
    self,
    formal_charges: Optional[Dict[str, float]] = None,
    layer_tolerance: float = 0.5
) -> List[Dict]:
    """
    Calculate formal charges per atomic layer.

    Parameters
    ----------
    formal_charges : dict, optional
        Formal charges by element (e.g., {"La": +3, "V": +5, "O": -2}).
        If None, uses common oxidation states.
    layer_tolerance : float
        Z-distance tolerance for layer detection (Angstrom)

    Returns
    -------
    list of dict
        List of {"z": z_position, "atoms": [...], "charge": total_charge}
    """
```

Default formal charges for common elements:
```python
DEFAULT_FORMAL_CHARGES = {
    "La": +3, "Sr": +2, "Ba": +2, "Ca": +2,  # A-site
    "V": +5, "Ti": +4, "Mn": +4, "Fe": +3, "Ni": +2, "Co": +3,  # B-site
    "O": -2, "N": -3,  # Anions
    "H": +1,  # Adsorbate
}
```

### 1.2 Add `check_polarity()` method

```python
def check_polarity(
    self,
    formal_charges: Optional[Dict[str, float]] = None
) -> Dict:
    """
    Check if surface is polar (has net dipole).

    Returns
    -------
    dict
        {
            "is_polar": bool,
            "dipole_moment": float (eÅ),
            "layer_charges": list,
            "top_layer_charge": float,
            "bottom_layer_charge": float,
            "recommendation": str
        }
    """
```

### 1.3 Add `calculate_dipole_moment()` method

```python
def calculate_dipole_moment(
    self,
    formal_charges: Optional[Dict[str, float]] = None
) -> Tuple[float, np.ndarray]:
    """
    Calculate electric dipole moment of the slab.

    Returns
    -------
    tuple
        (magnitude in eÅ, direction vector)
    """
```

---

## Phase 2: Symmetric Slab Creation

### 2.1 Add `create_symmetric_slab()` to SurfaceBuilder

```python
def create_symmetric_slab(
    self,
    miller_index: Tuple[int, int, int],
    layers: int = 6,
    vacuum: float = 15.0,
    termination: str = "AO",  # or "BO2" for perovskites
) -> SlabStructure:
    """
    Create a symmetric slab with same termination on both surfaces.

    This cancels the dipole moment for polar surfaces.

    Parameters
    ----------
    miller_index : tuple
        Miller indices
    layers : int
        Number of layers (will be adjusted to ensure symmetry)
    vacuum : float
        Vacuum thickness
    termination : str
        Desired termination: "AO" or "BO2" for perovskites

    Returns
    -------
    SlabStructure
        Symmetric slab with zero net dipole

    Notes
    -----
    For ABO3 perovskites along (001):
    - "AO" termination: ...AO-BO2-AO-BO2-AO... (odd number of AO layers)
    - "BO2" termination: ...BO2-AO-BO2-AO-BO2... (odd number of BO2 layers)
    """
```

### 2.2 Add `make_symmetric()` to SlabStructure

```python
def make_symmetric(self, method: str = "mirror") -> "SlabStructure":
    """
    Make an existing slab symmetric.

    Parameters
    ----------
    method : str
        "mirror": Mirror the top half to create symmetric slab
        "remove": Remove atoms from one side to balance charges

    Returns
    -------
    SlabStructure
        Symmetric slab
    """
```

---

## Phase 3: Perovskite-Specific Termination Selection

### 3.1 Add `PerovskiteSurfaceBuilder` class

```python
class PerovskiteSurfaceBuilder(SurfaceBuilder):
    """
    Specialized surface builder for ABO3 perovskite structures.

    Handles:
    - Named terminations (AO, BO2)
    - Automatic polarity compensation
    - Stoichiometry validation

    Examples
    --------
    >>> bulk = BulkStructure.from_cif("LaVO3.cif")
    >>> builder = PerovskiteSurfaceBuilder(bulk, A_site="La", B_site="V")
    >>> slab = builder.create_surface(
    ...     miller_index=(0, 0, 1),
    ...     termination="LaO",  # or "VO2"
    ...     layers=6,
    ...     symmetric=True
    ... )
    """

    def __init__(
        self,
        bulk: Union[BulkStructure, Atoms],
        A_site: Optional[str] = None,
        B_site: Optional[str] = None,
        anion: str = "O"
    ):
        """
        Parameters
        ----------
        bulk : BulkStructure or Atoms
            ABO3 perovskite bulk structure
        A_site : str, optional
            A-site cation (e.g., "La", "Sr"). Auto-detected if None.
        B_site : str, optional
            B-site cation (e.g., "V", "Ti"). Auto-detected if None.
        anion : str
            Anion species (default "O", can be "N" for oxynitrides)
        """
```

### 3.2 Key methods for PerovskiteSurfaceBuilder

```python
def identify_sites(self) -> Dict[str, List[int]]:
    """Identify A-site, B-site, and anion positions in bulk."""

def get_termination_options(
    self,
    miller_index: Tuple[int, int, int]
) -> List[str]:
    """Get available termination names for given Miller index."""

def create_surface(
    self,
    miller_index: Tuple[int, int, int],
    termination: str,  # "LaO", "VO2", "AO", "BO2"
    layers: int = 6,
    vacuum: float = 15.0,
    symmetric: bool = True,
    fix_bottom: int = 2,
) -> SlabStructure:
    """Create surface with specified termination."""

def select_termination(
    self,
    slab: SlabStructure,
    termination: str
) -> SlabStructure:
    """
    Adjust slab to have specified termination.

    May remove/add atoms from top or bottom.
    """
```

---

## Phase 4: Stoichiometry Validation

### 4.1 Add stoichiometry methods to SlabStructure

```python
def get_stoichiometry(self) -> Dict[str, float]:
    """Get normalized stoichiometry (e.g., La:V:O = 1:1:3)."""

def check_stoichiometry(
    self,
    expected: Optional[Dict[str, float]] = None
) -> Dict:
    """
    Check if stoichiometry is reasonable.

    Returns
    -------
    dict
        {
            "is_stoichiometric": bool,
            "actual": dict,
            "expected": dict,
            "deviation": dict,
            "warnings": list
        }
    """

def get_layer_stoichiometry(self) -> List[Dict]:
    """Get stoichiometry per layer."""
```

---

## Phase 5: Layer Analysis Improvements

### 5.1 Enhance layer detection

```python
def identify_layers(
    self,
    tolerance: float = 0.5,
    method: str = "z_clustering"
) -> List[Dict]:
    """
    Identify atomic layers in the slab.

    Parameters
    ----------
    tolerance : float
        Z-distance tolerance for clustering
    method : str
        "z_clustering": Simple z-position clustering
        "coordination": Use coordination environment

    Returns
    -------
    list of dict
        [{
            "z": z_position,
            "indices": [atom indices],
            "composition": {"La": 2, "O": 2},
            "layer_type": "AO" or "BO2" (for perovskites)
        }, ...]
    """

def get_layer_spacing(self) -> List[float]:
    """Get interlayer distances."""

def get_layer_by_type(self, layer_type: str) -> List[int]:
    """Get indices of all atoms in layers of given type (e.g., 'AO')."""
```

---

## Phase 6: Surface Reconstruction Support (Future)

### 6.1 Basic reconstruction patterns

```python
def apply_reconstruction(
    self,
    pattern: str = "sqrt2_x_sqrt2_R45"
) -> SlabStructure:
    """
    Apply common surface reconstruction.

    Parameters
    ----------
    pattern : str
        "sqrt2_x_sqrt2_R45": (√2×√2)R45° reconstruction
        "2x1": (2×1) reconstruction
        "2x2": (2×2) reconstruction
    """
```

---

## Implementation Order

1. **Phase 1**: Polarity analysis (foundation for everything else)
2. **Phase 2**: Symmetric slab creation (uses Phase 1)
3. **Phase 3**: Perovskite-specific builder (uses Phases 1-2)
4. **Phase 4**: Stoichiometry validation (independent, can parallel)
5. **Phase 5**: Layer analysis improvements (enhances all above)
6. **Phase 6**: Reconstruction (future, lower priority)

## Files to Modify

1. `nh3sofc/structure/surface.py` - Main enhancements
2. `nh3sofc/core/constants.py` - Add DEFAULT_FORMAL_CHARGES
3. `docs/tutorials/surface_building.md` - Update documentation
4. `docs/api/structure.md` - API documentation

## Testing Plan

1. Test polarity calculation on known polar surfaces (LaVO3 001)
2. Test symmetric slab creation - verify dipole is ~0
3. Test perovskite termination selection - verify surface composition
4. Test with real LaVO3 and SrTiO3 CIF files
5. Compare against literature surface energies (qualitative)

## Questions for User

1. **Priority**: Should I focus on perovskites only, or also handle other structures (rocksalt, fluorite)?
2. **Reconstruction**: Is (√2×√2)R45° reconstruction important for your NH3 work, or can it wait?
3. **Charge model**: Should we support Bader charges from VASP, or just formal charges?
