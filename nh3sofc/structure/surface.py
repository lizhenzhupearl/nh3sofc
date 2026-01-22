"""Surface generation from bulk structures."""

from pathlib import Path
from typing import Optional, List, Union, Tuple, Dict
from ase.data import covalent_radii, atomic_numbers
import numpy as np
from ase import Atoms
from ase.build import surface, add_vacuum
from ase.constraints import FixAtoms
from ase.io import write

from ..core.base import BaseStructure, get_bottom_atoms, get_surface_atoms
from ..core.constants import (
    DEFAULT_SURFACE_PARAMS,
    DEFAULT_FORMAL_CHARGES,
    PEROVSKITE_SITES,
    SURFACE_POLARITY,
)
from .bulk import BulkStructure


class SlabStructure(BaseStructure):
    """
    Class representing a surface slab.

    Stores information about the slab including Miller indices,
    termination, and fixed atom constraints.
    """

    def __init__(
        self,
        atoms: Atoms,
        miller_index: Optional[Tuple[int, int, int]] = None,
        termination: Optional[str] = None,
        n_layers: Optional[int] = None,
    ):
        """
        Initialize SlabStructure.

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object for the slab
        miller_index : tuple, optional
            Miller indices (h, k, l)
        termination : str, optional
            Surface termination identifier
        n_layers : int, optional
            Number of layers in the slab
        """
        super().__init__(atoms)
        self.miller_index = miller_index
        self.termination = termination
        self.n_layers = n_layers
        self._fixed_indices = []

    @property
    def fixed_indices(self) -> List[int]:
        """Get indices of fixed atoms."""
        return self._fixed_indices

    def fix_bottom_layers(self, n_layers: int, layer_threshold: float = 0.1) -> None:
        """
        Fix the bottom n layers of the slab.

        Parameters
        ----------
        n_layers : int
            Number of bottom layers to fix
        layer_threshold : float
            Tolerance for layer detection (fraction of interlayer spacing)
        """
        z_positions = self.atoms.positions[:, 2]
        unique_z = np.sort(np.unique(np.round(z_positions, decimals=1)))

        if len(unique_z) < n_layers:
            # If fewer unique z than layers, just fix bottom portion
            z_cutoff = unique_z[min(n_layers - 1, len(unique_z) - 1)]
        else:
            z_cutoff = unique_z[n_layers - 1] + layer_threshold

        fix_indices = [i for i, z in enumerate(z_positions) if z <= z_cutoff]
        self._fixed_indices = fix_indices

        # Apply constraint
        constraint = FixAtoms(indices=fix_indices)
        self.atoms.set_constraint(constraint)

    def fix_atoms(self, indices: List[int]) -> None:
        """
        Fix specific atoms by index.

        Parameters
        ----------
        indices : list
            Indices of atoms to fix
        """
        self._fixed_indices = list(indices)
        constraint = FixAtoms(indices=indices)
        self.atoms.set_constraint(constraint)

    def unfix_atoms(self) -> None:
        """Remove all constraints."""
        self._fixed_indices = []
        self.atoms.set_constraint([])

    def get_surface_atoms(self, z_threshold: float = 0.2) -> List[int]:
        """
        Get indices of surface atoms.

        Parameters
        ----------
        z_threshold : float
            Fraction of slab height considered as surface

        Returns
        -------
        list
            Indices of surface atoms
        """
        indices, _ = get_surface_atoms(self.atoms, z_threshold)
        return indices

    def get_surface_area(self) -> float:
        """
        Get the surface area of the slab.

        Returns
        -------
        float
            Surface area in Angstrom^2
        """
        cell = self.atoms.get_cell()
        # Cross product of a and b vectors
        area = np.linalg.norm(np.cross(cell[0], cell[1]))
        return area

    def get_vacuum_size(self) -> float:
        """
        Estimate the vacuum size in the slab.

        Returns
        -------
        float
            Vacuum size in Angstrom
        """
        z_positions = self.atoms.positions[:, 2]
        z_min, z_max = z_positions.min(), z_positions.max()
        slab_height = z_max - z_min
        cell_z = self.atoms.get_cell()[2, 2]
        vacuum = cell_z - slab_height
        return vacuum

    # ========== Layer Analysis Methods ==========

    def estimate_layer_tolerance(self, bond_scale: float = 0.5) -> float:
        """
        Estimate layer tolerance from covalent radii.

        Calculates a tolerance based on the minimum expected bond length
        in the structure. This is useful for complex materials like perovskites
        where atoms in the same layer may have different z-positions.

        Parameters
        ----------
        bond_scale : float
            Fraction of minimum bond length to use as tolerance (default: 0.5)

        Returns
        -------
        float
            Estimated tolerance in Angstrom

        Examples
        --------
        >>> slab.estimate_layer_tolerance()
        1.1  # For LaVO3 (V-O bond ~2.19 A, tolerance ~1.1 A)
        """
        symbols = self.atoms.get_chemical_symbols()
        unique_symbols = list(set(symbols))

        # Get covalent radii for all elements
        radii = {}
        for sym in unique_symbols:
            radii[sym] = covalent_radii[atomic_numbers[sym]]

        # Find minimum bond length (sum of two smallest radii)
        radii_list = list(radii.values())
        if len(radii_list) >= 2:
            sorted_radii = sorted(radii_list)
            min_bond = sorted_radii[0] + sorted_radii[1]
        else:
            min_bond = 2 * radii_list[0] if radii_list else 2.0

        return bond_scale * min_bond

    def identify_layers(
        self,
        tolerance: Union[float, str] = "auto",
    ) -> List[Dict]:
        """
        Identify atomic layers in the slab based on z-positions.

        Parameters
        ----------
        tolerance : float or "auto"
            Z-distance tolerance for layer clustering (Angstrom).
            If "auto", calculates tolerance from covalent radii using
            estimate_layer_tolerance(). Default is "auto" to correctly
            group atoms in complex materials like perovskites.

        Returns
        -------
        list of dict
            List of layer info dicts sorted by z-position (bottom to top):
            [{
                "z": z_position,
                "indices": [atom indices],
                "symbols": [chemical symbols],
                "composition": {"La": 2, "O": 2},
            }, ...]

        Examples
        --------
        >>> layers = slab.identify_layers()  # auto tolerance
        >>> layers = slab.identify_layers(tolerance=0.5)  # explicit tolerance
        """
        # Handle auto tolerance
        if tolerance == "auto":
            tolerance = self.estimate_layer_tolerance()

        z_positions = self.atoms.positions[:, 2]
        symbols = self.atoms.get_chemical_symbols()

        # Cluster atoms by z-position
        sorted_indices = np.argsort(z_positions)
        layers = []
        current_layer_indices = [sorted_indices[0]]
        current_z = z_positions[sorted_indices[0]]

        for idx in sorted_indices[1:]:
            z = z_positions[idx]
            if abs(z - current_z) <= tolerance:
                current_layer_indices.append(idx)
            else:
                # Finish current layer
                layer_symbols = [symbols[i] for i in current_layer_indices]
                composition = {}
                for s in layer_symbols:
                    composition[s] = composition.get(s, 0) + 1

                layers.append({
                    "z": np.mean([z_positions[i] for i in current_layer_indices]),
                    "indices": current_layer_indices,
                    "symbols": layer_symbols,
                    "composition": composition,
                })
                # Start new layer
                current_layer_indices = [idx]
                current_z = z

        # Don't forget the last layer
        layer_symbols = [symbols[i] for i in current_layer_indices]
        composition = {}
        for s in layer_symbols:
            composition[s] = composition.get(s, 0) + 1
        layers.append({
            "z": np.mean([z_positions[i] for i in current_layer_indices]),
            "indices": current_layer_indices,
            "symbols": layer_symbols,
            "composition": composition,
        })

        return layers

    def get_layer_spacing(self) -> List[float]:
        """
        Get interlayer distances.

        Returns
        -------
        list of float
            Distances between consecutive layers (Angstrom)
        """
        layers = self.identify_layers()
        spacings = []
        for i in range(len(layers) - 1):
            spacing = layers[i + 1]["z"] - layers[i]["z"]
            spacings.append(spacing)
        return spacings

    def get_top_layer(self, tolerance: Union[float, str] = "auto") -> Dict:
        """Get the topmost layer."""
        layers = self.identify_layers(tolerance)
        return layers[-1] if layers else {}

    def get_bottom_layer(self, tolerance: Union[float, str] = "auto") -> Dict:
        """Get the bottommost layer."""
        layers = self.identify_layers(tolerance)
        return layers[0] if layers else {}

    # ========== Polarity Analysis Methods ==========

    def calculate_layer_charges(
        self,
        formal_charges: Optional[Dict[str, float]] = None,
        tolerance: Union[float, str] = "auto",
    ) -> List[Dict]:
        """
        Calculate formal charges per atomic layer.

        Parameters
        ----------
        formal_charges : dict, optional
            Formal charges by element (e.g., {"La": +3, "V": +5, "O": -2}).
            If None, uses DEFAULT_FORMAL_CHARGES.
        tolerance : float or "auto"
            Z-distance tolerance for layer detection (Angstrom)

        Returns
        -------
        list of dict
            List of layer info with charges:
            [{
                "z": z_position,
                "indices": [atom indices],
                "composition": {"La": 2, "O": 2},
                "charge": total_charge,
                "charge_per_area": charge/area,
            }, ...]
        """
        if formal_charges is None:
            formal_charges = DEFAULT_FORMAL_CHARGES

        layers = self.identify_layers(tolerance)
        area = self.get_surface_area()

        for layer in layers:
            total_charge = 0.0
            for symbol, count in layer["composition"].items():
                charge = formal_charges.get(symbol, 0.0)
                total_charge += charge * count
            layer["charge"] = total_charge
            layer["charge_per_area"] = total_charge / area if area > 0 else 0.0

        return layers

    def calculate_dipole_moment(
        self,
        formal_charges: Optional[Dict[str, float]] = None,
    ) -> Tuple[float, np.ndarray]:
        """
        Calculate electric dipole moment of the slab.

        The dipole moment is calculated as sum of q_i * r_i for all atoms.

        Parameters
        ----------
        formal_charges : dict, optional
            Formal charges by element. If None, uses DEFAULT_FORMAL_CHARGES.

        Returns
        -------
        tuple
            (magnitude in e·Å, direction vector [normalized])
        """
        if formal_charges is None:
            formal_charges = DEFAULT_FORMAL_CHARGES

        symbols = self.atoms.get_chemical_symbols()
        positions = self.atoms.get_positions()

        # Calculate center of mass for reference
        masses = self.atoms.get_masses()
        com = np.sum(positions * masses[:, np.newaxis], axis=0) / np.sum(masses)

        # Calculate dipole moment relative to center of mass
        dipole = np.zeros(3)
        for i, (symbol, pos) in enumerate(zip(symbols, positions)):
            charge = formal_charges.get(symbol, 0.0)
            dipole += charge * (pos - com)

        magnitude = np.linalg.norm(dipole)
        direction = dipole / magnitude if magnitude > 1e-10 else np.array([0, 0, 1])

        return magnitude, direction

    def check_polarity(
        self,
        formal_charges: Optional[Dict[str, float]] = None,
        tolerance: Union[float, str] = "auto",
    ) -> Dict:
        """
        Check if surface is polar (has net dipole in z-direction).

        Parameters
        ----------
        formal_charges : dict, optional
            Formal charges by element.
        tolerance : float or "auto"
            Layer detection tolerance.

        Returns
        -------
        dict
            {
                "is_polar": bool,
                "dipole_moment": float (e·Å),
                "dipole_z": float (z-component),
                "top_layer_charge": float,
                "bottom_layer_charge": float,
                "layer_charges": list,
                "recommendation": str,
            }
        """
        layers = self.calculate_layer_charges(formal_charges, tolerance)
        magnitude, direction = self.calculate_dipole_moment(formal_charges)
        dipole_z = magnitude * direction[2]

        top_charge = layers[-1]["charge"] if layers else 0.0
        bottom_charge = layers[0]["charge"] if layers else 0.0

        # Check if alternating charged layers (typical for polar surfaces)
        charges = [layer["charge"] for layer in layers]
        has_alternating = False
        if len(charges) > 1:
            signs = [np.sign(c) for c in charges if abs(c) > 0.1]
            if len(signs) > 1:
                has_alternating = all(signs[i] != signs[i+1] for i in range(len(signs)-1))

        # Determine if polar (significant z-dipole or alternating charges)
        is_polar = abs(dipole_z) > 0.5 or has_alternating

        # Generate recommendation
        if is_polar:
            if abs(top_charge + bottom_charge) < 0.1:
                recommendation = "Surface has symmetric terminations - dipole may cancel."
            else:
                recommendation = (
                    "Polar surface detected. Consider: "
                    "(1) Use symmetric slab with same termination on both sides, "
                    "(2) Apply dipole correction (IDIPOL=3, LDIPOL=.TRUE.), or "
                    "(3) Add compensating charges/vacancies."
                )
        else:
            recommendation = "Non-polar surface - no special treatment needed."

        return {
            "is_polar": is_polar,
            "dipole_moment": magnitude,
            "dipole_z": dipole_z,
            "top_layer_charge": top_charge,
            "bottom_layer_charge": bottom_charge,
            "layer_charges": charges,
            "has_alternating_charges": has_alternating,
            "recommendation": recommendation,
        }

    # ========== Stoichiometry Methods ==========

    def get_stoichiometry(self) -> Dict[str, float]:
        """
        Get normalized stoichiometry.

        Returns
        -------
        dict
            Element ratios normalized to smallest count
            (e.g., {"La": 1.0, "V": 1.0, "O": 3.0})
        """
        counts = self.get_element_count()
        if not counts:
            return {}

        min_count = min(counts.values())
        return {elem: count / min_count for elem, count in counts.items()}

    def check_stoichiometry(
        self,
        expected: Optional[Dict[str, float]] = None,
        tolerance: float = 0.1,
    ) -> Dict:
        """
        Check if stoichiometry is reasonable.

        Parameters
        ----------
        expected : dict, optional
            Expected stoichiometry (e.g., {"La": 1, "V": 1, "O": 3})
        tolerance : float
            Allowed deviation from expected ratios

        Returns
        -------
        dict
            {
                "is_stoichiometric": bool,
                "actual": dict,
                "expected": dict or None,
                "deviation": dict,
                "warnings": list,
            }
        """
        actual = self.get_stoichiometry()
        warnings = []

        if expected is None:
            # Can't check against expected
            return {
                "is_stoichiometric": True,  # Assume OK
                "actual": actual,
                "expected": None,
                "deviation": {},
                "warnings": ["No expected stoichiometry provided for comparison."],
            }

        # Calculate deviations
        deviation = {}
        is_stoichiometric = True

        for elem in set(list(actual.keys()) + list(expected.keys())):
            act = actual.get(elem, 0.0)
            exp = expected.get(elem, 0.0)
            if exp > 0:
                dev = abs(act - exp) / exp
            else:
                dev = float('inf') if act > 0 else 0.0
            deviation[elem] = dev

            if dev > tolerance:
                is_stoichiometric = False
                warnings.append(
                    f"{elem}: expected {exp:.2f}, got {act:.2f} (deviation: {dev*100:.1f}%)"
                )

        return {
            "is_stoichiometric": is_stoichiometric,
            "actual": actual,
            "expected": expected,
            "deviation": deviation,
            "warnings": warnings,
        }

    def get_layer_stoichiometry(self, tolerance: Union[float, str] = "auto") -> List[Dict]:
        """
        Get stoichiometry per layer.

        Returns
        -------
        list of dict
            [{
                "z": z_position,
                "composition": {"La": 2, "O": 2},
                "stoichiometry": {"La": 1.0, "O": 1.0},
            }, ...]
        """
        layers = self.identify_layers(tolerance)

        for layer in layers:
            comp = layer["composition"]
            if comp:
                min_count = min(comp.values())
                layer["stoichiometry"] = {
                    elem: count / min_count for elem, count in comp.items()
                }
            else:
                layer["stoichiometry"] = {}

        return layers

    # ========== Symmetric Slab Methods ==========

    def _composition_matches_ratio(
        self, layer_comp: Dict[str, int], target: Dict[str, int]
    ) -> bool:
        """
        Check if a layer composition matches a target by element ratio.

        Parameters
        ----------
        layer_comp : dict
            Layer composition (e.g., {"V": 8, "O": 16})
        target : dict
            Target composition ratio (e.g., {"V": 1, "O": 2})

        Returns
        -------
        bool
            True if element ratios match
        """
        if not layer_comp or not target:
            return False

        # Must have exactly the same elements
        if set(layer_comp.keys()) != set(target.keys()):
            return False

        # Calculate ratios for both
        layer_min = min(layer_comp.values())
        target_min = min(target.values())

        layer_ratio = {k: v / layer_min for k, v in layer_comp.items()}
        target_ratio = {k: v / target_min for k, v in target.items()}

        # Check if ratios match (within tolerance)
        for elem in layer_ratio:
            if abs(layer_ratio[elem] - target_ratio[elem]) > 0.1:
                return False
        return True

    def _find_matching_layers(
        self, layers: List[Dict], termination: Dict[str, int]
    ) -> List[int]:
        """
        Find layer indices that match the target termination.

        Parameters
        ----------
        layers : list of dict
            Layer info from identify_layers()
        termination : dict
            Target termination composition (e.g., {"La": 1, "O": 1})

        Returns
        -------
        list of int
            Indices of matching layers (sorted by z from bottom to top)
        """
        matching = []
        for i, layer in enumerate(layers):
            if self._composition_matches_ratio(layer["composition"], termination):
                matching.append(i)
        return matching

    def trim_to_symmetric_termination(
        self,
        termination: Optional[Dict[str, int]] = None,
        termination_name: Optional[str] = None,
        min_layers: int = 5,
        tolerance: Union[float, str] = "auto",
    ) -> "SlabStructure":
        """
        Trim slab to ensure matching top/bottom terminations.

        Creates a truly symmetric slab by finding layers that match the
        target termination and trimming to the outermost matching layers.

        Parameters
        ----------
        termination : dict, optional
            Target termination composition as element ratios
            (e.g., {"La": 1, "O": 1} for LaO or {"V": 1, "O": 2} for VO2).
            If None, uses the top layer's composition as target.
        termination_name : str, optional
            Name for the termination (e.g., "LaO", "VO2")
        min_layers : int
            Minimum number of layers to keep (default: 5)
        tolerance : float or "auto"
            Layer detection tolerance

        Returns
        -------
        SlabStructure
            New slab with matching top/bottom terminations

        Examples
        --------
        >>> # Create LaO-terminated symmetric slab
        >>> symmetric = slab.trim_to_symmetric_termination(
        ...     termination={"La": 1, "O": 1}, min_layers=5
        ... )
        >>> layers = symmetric.identify_layers()
        >>> print(layers[0]["composition"])  # Bottom: LaO
        >>> print(layers[-1]["composition"])  # Top: also LaO
        """
        layers = self.identify_layers(tolerance)

        if len(layers) < min_layers:
            raise ValueError(
                f"Slab has only {len(layers)} layers, need at least {min_layers} "
                "for symmetric trimming. Create a larger slab first."
            )

        # If no termination specified, use top layer composition as target
        if termination is None:
            termination = layers[-1]["composition"]

        # Find layers matching the target termination
        matching = self._find_matching_layers(layers, termination)

        if len(matching) < 2:
            raise ValueError(
                f"Found only {len(matching)} layers matching termination {termination}. "
                "Need at least 2 (one for top, one for bottom). "
                "Try a different termination or larger slab."
            )

        # Find the best top/bottom pair with at least min_layers between them
        best_bottom = None
        best_top = None
        best_n_layers = 0

        for bottom_idx in matching:
            for top_idx in reversed(matching):
                if top_idx <= bottom_idx:
                    continue
                n_layers = top_idx - bottom_idx + 1
                if n_layers >= min_layers and n_layers > best_n_layers:
                    best_bottom = bottom_idx
                    best_top = top_idx
                    best_n_layers = n_layers
                    break  # Found best top for this bottom

        if best_bottom is None or best_top is None:
            # Try with relaxed min_layers
            for bottom_idx in matching:
                for top_idx in reversed(matching):
                    if top_idx > bottom_idx:
                        best_bottom = bottom_idx
                        best_top = top_idx
                        best_n_layers = top_idx - bottom_idx + 1
                        break
                if best_bottom is not None:
                    break

        if best_bottom is None or best_top is None:
            raise ValueError(
                f"Could not find valid top/bottom pair with termination {termination}. "
                f"Matching layers: {matching}"
            )

        # Get atom indices to keep (all atoms in layers from bottom to top inclusive)
        keep_indices = []
        for i in range(best_bottom, best_top + 1):
            keep_indices.extend(layers[i]["indices"])

        # Create new slab with only the selected atoms
        new_atoms = self.atoms[keep_indices]
        new_atoms.center(axis=2)

        # Determine termination name
        if termination_name is None:
            # Try to generate from composition
            parts = [f"{elem}{count}" if count > 1 else elem
                     for elem, count in sorted(termination.items())]
            termination_name = "".join(parts)

        new_slab = SlabStructure(
            new_atoms,
            miller_index=self.miller_index,
            termination=f"{termination_name}_symmetric",
            n_layers=best_n_layers,
        )

        return new_slab

    def add_vacuum(self, vacuum: float) -> "SlabStructure":
        """
        Add vacuum to the slab.

        Parameters
        ----------
        vacuum : float
            Vacuum to add (Angstrom)

        Returns
        -------
        SlabStructure
            New slab with added vacuum
        """
        new_atoms = self.atoms.copy()
        add_vacuum(new_atoms, vacuum)
        return SlabStructure(
            new_atoms,
            miller_index=self.miller_index,
            termination=self.termination,
            n_layers=self.n_layers
        )

    def center_in_cell(self, axis: int = 2) -> "SlabStructure":
        """
        Center the slab in the cell.

        Parameters
        ----------
        axis : int
            Axis to center (default: 2 for z)

        Returns
        -------
        SlabStructure
            Centered slab
        """
        new_atoms = self.atoms.copy()
        new_atoms.center(axis=axis)
        new_slab = SlabStructure(
            new_atoms,
            miller_index=self.miller_index,
            termination=self.termination,
            n_layers=self.n_layers
        )
        new_slab._fixed_indices = self._fixed_indices.copy()
        return new_slab

    def repeat_xy(self, nx: int, ny: int) -> "SlabStructure":
        """
        Repeat the slab in x and y directions.

        Parameters
        ----------
        nx, ny : int
            Number of repetitions

        Returns
        -------
        SlabStructure
            Repeated slab (with constraints reapplied if present)
        """
        n_atoms_original = len(self.atoms)
        new_atoms = self.atoms.repeat((nx, ny, 1))
        new_slab = SlabStructure(
            new_atoms,
            miller_index=self.miller_index,
            termination=self.termination,
            n_layers=self.n_layers
        )

        # Reapply constraints if there were fixed atoms
        if self._fixed_indices:
            # Map old fixed indices to new indices in repeated structure
            new_fixed = []
            for rep in range(nx * ny):
                for old_idx in self._fixed_indices:
                    new_fixed.append(rep * n_atoms_original + old_idx)
            new_slab.fix_atoms(new_fixed)

        return new_slab

    def copy(self) -> "SlabStructure":
        """Return a copy of this slab."""
        new_slab = SlabStructure(
            self.atoms.copy(),
            miller_index=self.miller_index,
            termination=self.termination,
            n_layers=self.n_layers
        )
        new_slab._fixed_indices = self._fixed_indices.copy()
        return new_slab

    def __repr__(self) -> str:
        mi = "".join(map(str, self.miller_index)) if self.miller_index else "?"
        return (
            f"SlabStructure({self.formula}, "
            f"({mi}), "
            f"{self.n_atoms} atoms, "
            f"area={self.get_surface_area():.1f} A^2)"
        )


class SurfaceBuilder:
    """
    Builder class for creating surface slabs from bulk structures.

    Examples
    --------
    >>> bulk = BulkStructure.from_cif("LaVO3.cif")
    >>> builder = SurfaceBuilder(bulk)
    >>> slab = builder.create_surface(
    ...     miller_index=(0, 0, 1),
    ...     layers=6,
    ...     vacuum=15.0
    ... )
    >>> slab.fix_bottom_layers(2)
    """

    def __init__(self, bulk: Union[BulkStructure, Atoms]):
        """
        Initialize SurfaceBuilder.

        Parameters
        ----------
        bulk : BulkStructure or Atoms
            Bulk structure to create surface from
        """
        if isinstance(bulk, BulkStructure):
            self.bulk = bulk.atoms.copy()
        else:
            self.bulk = bulk.copy()

    def create_surface(
        self,
        miller_index: Tuple[int, int, int],
        layers: int = 6,
        vacuum: float = 15.0,
        termination: Optional[int] = None,
        periodic: bool = True,
        fix_bottom: int = 0,
    ) -> SlabStructure:
        """
        Create a surface slab.

        Parameters
        ----------
        miller_index : tuple
            Miller indices (h, k, l)
        layers : int
            Number of layers
        vacuum : float
            Vacuum thickness in Angstrom
        termination : int, optional
            Termination index (for multi-termination surfaces)
        periodic : bool
            Whether to make the surface periodic
        fix_bottom : int
            Number of bottom layers to fix (constrain). Default is 0 (no fixing).

        Returns
        -------
        SlabStructure
            Surface slab
        """
        # Create surface using ASE
        slab = surface(
            self.bulk,
            indices=miller_index,
            layers=layers,
            vacuum=vacuum / 2,  # ASE adds vacuum on both sides
            periodic=periodic
        )

        # Add remaining vacuum to top
        slab.center(axis=2)

        # Determine termination label
        term_label = None
        if termination is not None:
            term_label = f"term_{termination}"

        slab_structure = SlabStructure(
            slab,
            miller_index=miller_index,
            termination=term_label,
            n_layers=layers
        )

        # Fix bottom layers if requested
        if fix_bottom > 0:
            slab_structure.fix_bottom_layers(fix_bottom)

        return slab_structure

    def create_surface_with_size(
        self,
        miller_index: Tuple[int, int, int],
        size: Tuple[int, int, int],
        vacuum: float = 15.0,
    ) -> SlabStructure:
        """
        Create a surface slab with specific supercell size.

        Parameters
        ----------
        miller_index : tuple
            Miller indices (h, k, l)
        size : tuple
            Supercell size (nx, ny, nz) where nz is layers
        vacuum : float
            Vacuum thickness in Angstrom

        Returns
        -------
        SlabStructure
            Surface slab
        """
        nx, ny, layers = size

        # Create surface
        slab = surface(
            self.bulk,
            indices=miller_index,
            layers=layers,
            vacuum=vacuum / 2,
            periodic=True
        )

        # Repeat in x and y
        slab = slab.repeat((nx, ny, 1))
        slab.center(axis=2)

        return SlabStructure(
            slab,
            miller_index=miller_index,
            n_layers=layers
        )

    def get_all_terminations(
        self,
        miller_index: Tuple[int, int, int],
        layers: int = 6,
        vacuum: float = 15.0,
    ) -> List[SlabStructure]:
        """
        Generate all possible surface terminations.

        This shifts the bulk cell to expose different atomic planes.

        Parameters
        ----------
        miller_index : tuple
            Miller indices
        layers : int
            Number of layers
        vacuum : float
            Vacuum thickness

        Returns
        -------
        list
            List of SlabStructure objects for each termination
        """
        terminations = []

        # Get unique z-layers in the bulk
        bulk_z = self.bulk.positions[:, 2]
        unique_z = np.sort(np.unique(np.round(bulk_z, decimals=2)))
        n_unique = len(unique_z)

        for i, z_shift in enumerate(unique_z):
            # Create shifted bulk
            shifted_bulk = self.bulk.copy()
            shifted_bulk.positions[:, 2] -= z_shift
            shifted_bulk.wrap()

            # Create surface from shifted bulk
            slab = surface(
                shifted_bulk,
                indices=miller_index,
                layers=layers,
                vacuum=vacuum / 2,
                periodic=True
            )
            slab.center(axis=2)

            term_slab = SlabStructure(
                slab,
                miller_index=miller_index,
                termination=f"term_{i}",
                n_layers=layers
            )
            terminations.append(term_slab)

        return terminations

    def identify_termination(self, slab: SlabStructure) -> dict:
        """
        Identify the termination of a slab based on surface composition.

        Parameters
        ----------
        slab : SlabStructure
            Slab to analyze

        Returns
        -------
        dict
            Termination information including surface composition
        """
        surface_indices = slab.get_surface_atoms(z_threshold=0.15)
        surface_symbols = [slab.atoms[i].symbol for i in surface_indices]

        # Count surface atoms
        composition = {}
        for symbol in surface_symbols:
            composition[symbol] = composition.get(symbol, 0) + 1

        # Determine dominant termination
        dominant = max(composition.items(), key=lambda x: x[1])

        return {
            "surface_indices": surface_indices,
            "composition": composition,
            "dominant_element": dominant[0],
            "dominant_count": dominant[1],
        }

    def create_symmetric_slab(
        self,
        miller_index: Tuple[int, int, int],
        layers: int = 7,
        vacuum: float = 15.0,
        termination: Optional[Dict[str, int]] = None,
        min_layers: int = 5,
        fix_bottom: int = 2,
        tolerance: Union[float, str] = "auto",
    ) -> SlabStructure:
        """
        Create a symmetric slab with same termination on top and bottom.

        Creates an oversized slab and trims it to ensure the top and bottom
        layers have matching compositions. This is essential for polar
        surfaces to cancel the dipole moment.

        Parameters
        ----------
        miller_index : tuple
            Miller indices (h, k, l)
        layers : int
            Target number of layers (actual may differ slightly to achieve symmetry)
        vacuum : float
            Vacuum thickness in Angstrom
        termination : dict, optional
            Target termination composition as element ratios
            (e.g., {"La": 1, "O": 1} for LaO-terminated).
            If None, uses the top layer composition of the initial slab.
        min_layers : int
            Minimum number of layers to keep (default: 5)
        fix_bottom : int
            Number of bottom layers to fix
        tolerance : float or "auto"
            Layer detection tolerance

        Returns
        -------
        SlabStructure
            Symmetric slab with matching top/bottom terminations

        Notes
        -----
        For perovskites (ABO3) along (001), you can request either:
        - AO termination: termination={"La": 1, "O": 1}
        - BO2 termination: termination={"V": 1, "O": 2}

        Examples
        --------
        >>> builder = SurfaceBuilder(bulk)
        >>> # LaO-terminated symmetric slab
        >>> slab = builder.create_symmetric_slab(
        ...     (0, 0, 1), layers=7,
        ...     termination={"La": 1, "O": 1}
        ... )
        """
        # Create oversized slab to have enough layers for trimming
        oversized_layers = max(layers + 4, 10)

        # Use base SurfaceBuilder.create_surface to avoid recursion in subclasses
        oversized = SurfaceBuilder.create_surface(
            self,
            miller_index=miller_index,
            layers=oversized_layers,
            vacuum=vacuum,
            fix_bottom=0,  # Will apply constraints after trimming
        )

        # If no termination specified, find best one automatically
        if termination is None:
            slab = self._find_best_symmetric_termination(
                oversized, min_layers, tolerance
            )
            if slab is None:
                # Try with larger slab
                oversized = SurfaceBuilder.create_surface(
                    self,
                    miller_index=miller_index,
                    layers=oversized_layers + 6,
                    vacuum=vacuum,
                    fix_bottom=0,
                )
                slab = self._find_best_symmetric_termination(
                    oversized, min_layers, tolerance
                )
            if slab is None:
                raise ValueError(
                    "Could not find any termination that creates a symmetric slab. "
                    "Try specifying a termination explicitly, e.g., "
                    "termination={'La': 1, 'O': 1} or termination={'V': 1, 'O': 2}"
                )
        else:
            # Use specified termination
            try:
                slab = oversized.trim_to_symmetric_termination(
                    termination=termination,
                    min_layers=min_layers,
                    tolerance=tolerance,
                )
            except ValueError:
                # Fallback: try with even larger slab
                oversized = SurfaceBuilder.create_surface(
                    self,
                    miller_index=miller_index,
                    layers=oversized_layers + 6,
                    vacuum=vacuum,
                    fix_bottom=0,
                )
                slab = oversized.trim_to_symmetric_termination(
                    termination=termination,
                    min_layers=min_layers,
                    tolerance=tolerance,
                )

        # Add vacuum and center
        slab = slab.add_vacuum(vacuum)
        slab = slab.center_in_cell()

        # Set termination name (use "symmetric" if no specific termination was requested)
        if termination is None:
            slab.termination = "symmetric"

        # Apply constraints
        if fix_bottom > 0:
            slab.fix_bottom_layers(fix_bottom)

        return slab

    def _find_best_symmetric_termination(
        self,
        slab: SlabStructure,
        min_layers: int,
        tolerance: Union[float, str],
    ) -> Optional[SlabStructure]:
        """
        Try all unique layer compositions and find the best symmetric slab.

        Returns the symmetric slab with the most layers, or None if no
        valid symmetric termination is found.
        """
        layers = slab.identify_layers(tolerance)

        # Get unique layer compositions
        unique_compositions = []
        seen_ratios = []
        for layer in layers:
            comp = layer["composition"]
            if not comp:
                continue
            # Normalize to ratio for comparison
            min_val = min(comp.values())
            ratio = tuple(sorted((k, v / min_val) for k, v in comp.items()))
            if ratio not in seen_ratios:
                seen_ratios.append(ratio)
                unique_compositions.append(comp)

        # Try each composition and keep the best result
        best_slab = None
        best_n_layers = 0

        for comp in unique_compositions:
            try:
                trimmed = slab.trim_to_symmetric_termination(
                    termination=comp,
                    min_layers=min_layers,
                    tolerance=tolerance,
                )
                n_layers = len(trimmed.identify_layers(tolerance))
                if n_layers > best_n_layers:
                    best_n_layers = n_layers
                    best_slab = trimmed
            except ValueError:
                continue

        return best_slab


class PerovskiteSurfaceBuilder(SurfaceBuilder):
    """
    Specialized surface builder for ABO3 perovskite structures.

    Handles named terminations (AO, BO2), automatic polarity compensation,
    and stoichiometry validation.

    Examples
    --------
    >>> bulk = BulkStructure.from_cif("LaVO3.cif")
    >>> builder = PerovskiteSurfaceBuilder(bulk, A_site="La", B_site="V")
    >>> slab = builder.create_surface(
    ...     miller_index=(0, 0, 1),
    ...     termination="LaO",
    ...     layers=7,
    ...     symmetric=True
    ... )
    """

    def __init__(
        self,
        bulk: Union[BulkStructure, Atoms],
        A_site: Optional[str] = None,
        B_site: Optional[str] = None,
        anion: str = "O",
    ):
        """
        Initialize PerovskiteSurfaceBuilder.

        Parameters
        ----------
        bulk : BulkStructure or Atoms
            ABO3 perovskite bulk structure
        A_site : str, optional
            A-site cation (e.g., "La", "Sr"). Auto-detected if None.
        B_site : str, optional
            B-site cation (e.g., "V", "Ti"). Auto-detected if None.
        anion : str
            Anion species (default "O")
        """
        super().__init__(bulk)
        self.anion = anion

        # Auto-detect A and B sites if not provided
        symbols = list(set(self.bulk.get_chemical_symbols()))
        symbols = [s for s in symbols if s != anion]

        if A_site is None or B_site is None:
            detected = self._detect_sites(symbols)
            A_site = A_site or detected.get("A_site")
            B_site = B_site or detected.get("B_site")

        self.A_site = A_site
        self.B_site = B_site

    def _detect_sites(self, symbols: List[str]) -> Dict[str, str]:
        """Auto-detect A and B site cations."""
        A_candidates = [s for s in symbols if s in PEROVSKITE_SITES["A_site"]]
        B_candidates = [s for s in symbols if s in PEROVSKITE_SITES["B_site"]]

        result = {}
        if A_candidates:
            result["A_site"] = A_candidates[0]
        if B_candidates:
            result["B_site"] = B_candidates[0]

        return result

    def get_termination_options(
        self,
        miller_index: Tuple[int, int, int] = (0, 0, 1),
    ) -> List[str]:
        """
        Get available termination names for given Miller index.

        Parameters
        ----------
        miller_index : tuple
            Miller indices

        Returns
        -------
        list of str
            Available termination names (e.g., ["LaO", "VO2"])
        """
        if miller_index == (0, 0, 1):
            return [
                f"{self.A_site}{self.anion}",  # e.g., "LaO"
                f"{self.B_site}{self.anion}2",  # e.g., "VO2"
            ]
        elif miller_index == (1, 1, 0):
            return [
                f"{self.A_site}{self.B_site}{self.anion}",
                f"{self.anion}2",
            ]
        else:
            return ["term_0", "term_1"]

    def _get_termination_composition(
        self, termination: Optional[str]
    ) -> Optional[Dict[str, int]]:
        """
        Convert named termination to composition dict.

        Parameters
        ----------
        termination : str or None
            Named termination (e.g., "LaO", "VO2")

        Returns
        -------
        dict or None
            Composition as element ratios (e.g., {"La": 1, "O": 1})
        """
        if termination is None:
            return None

        # AO termination (e.g., LaO)
        if termination == f"{self.A_site}{self.anion}":
            return {self.A_site: 1, self.anion: 1}
        # BO2 termination (e.g., VO2)
        elif termination == f"{self.B_site}{self.anion}2":
            return {self.B_site: 1, self.anion: 2}

        return None

    def create_surface(
        self,
        miller_index: Tuple[int, int, int],
        termination: Optional[str] = None,
        layers: int = 7,
        vacuum: float = 15.0,
        symmetric: bool = True,
        fix_bottom: int = 2,
        supercell: Optional[Tuple[int, int]] = None,
    ) -> SlabStructure:
        """
        Create perovskite surface with specified termination.

        Parameters
        ----------
        miller_index : tuple
            Miller indices (h, k, l)
        termination : str, optional
            Desired termination (e.g., "LaO", "VO2", "AO", "BO2")
            If None, uses first available termination.
        layers : int
            Number of layers
        vacuum : float
            Vacuum thickness in Angstrom
        symmetric : bool
            If True, create symmetric slab (same termination on both sides)
        fix_bottom : int
            Number of bottom layers to fix
        supercell : tuple, optional
            Supercell size (nx, ny)

        Returns
        -------
        SlabStructure
            Perovskite surface slab
        """
        # Normalize termination name
        if termination:
            termination = termination.replace("AO", f"{self.A_site}{self.anion}")
            termination = termination.replace("BO2", f"{self.B_site}{self.anion}2")

        # Get all terminations and find the matching one
        all_terms = self.get_all_terminations(miller_index, layers, vacuum)

        best_slab = None
        best_match = 0

        for slab in all_terms:
            term_info = self.identify_termination(slab)
            comp = term_info["composition"]

            # Check if this termination matches the requested one
            if termination:
                if termination == f"{self.A_site}{self.anion}":
                    # Looking for AO termination
                    if self.A_site in comp and self.anion in comp:
                        match = comp.get(self.A_site, 0) + comp.get(self.anion, 0)
                        if match > best_match:
                            best_match = match
                            best_slab = slab
                elif termination == f"{self.B_site}{self.anion}2":
                    # Looking for BO2 termination
                    if self.B_site in comp:
                        match = comp.get(self.B_site, 0)
                        if match > best_match:
                            best_match = match
                            best_slab = slab
            else:
                # No specific termination requested, use first
                best_slab = slab
                break

        if best_slab is None:
            # Fall back to basic surface creation
            best_slab = super().create_surface(
                miller_index, layers, vacuum, fix_bottom=fix_bottom
            )

        best_slab.termination = termination

        # Create symmetric slab if requested
        if symmetric:
            # Convert named termination to composition dict
            term_dict = self._get_termination_composition(termination)

            # Use the trimming approach for true symmetry
            best_slab = super().create_symmetric_slab(
                miller_index=miller_index,
                layers=layers,
                vacuum=vacuum,
                termination=term_dict,
                min_layers=max(layers - 2, 5),
                fix_bottom=fix_bottom,
            )
            best_slab.termination = f"{termination}_symmetric" if termination else "symmetric"

        # Apply supercell
        if supercell:
            best_slab = best_slab.repeat_xy(supercell[0], supercell[1])

        # Fix bottom layers (may already be applied by create_symmetric_slab)
        if fix_bottom > 0 and not best_slab.fixed_indices:
            best_slab.fix_bottom_layers(fix_bottom)

        return best_slab

    def analyze_surface(self, slab: SlabStructure) -> Dict:
        """
        Analyze perovskite surface for termination, polarity, and stoichiometry.

        Parameters
        ----------
        slab : SlabStructure
            Perovskite slab to analyze

        Returns
        -------
        dict
            Comprehensive analysis including termination, polarity,
            stoichiometry, and recommendations
        """
        term_info = self.identify_termination(slab)
        polarity = slab.check_polarity()
        stoich = slab.check_stoichiometry(
            expected={self.A_site: 1, self.B_site: 1, self.anion: 3}
        )

        # Determine termination type
        comp = term_info["composition"]
        if self.A_site in comp and self.anion in comp and self.B_site not in comp:
            term_type = f"{self.A_site}{self.anion}"
        elif self.B_site in comp:
            term_type = f"{self.B_site}{self.anion}2"
        else:
            term_type = "mixed"

        return {
            "termination_type": term_type,
            "surface_composition": comp,
            "polarity": polarity,
            "stoichiometry": stoich,
            "A_site": self.A_site,
            "B_site": self.B_site,
            "anion": self.anion,
        }


class RocksaltSurfaceBuilder(SurfaceBuilder):
    """
    Specialized surface builder for rocksalt (MX) structures like NiO, MgO.

    Examples
    --------
    >>> bulk = BulkStructure.from_cif("NiO.cif")
    >>> builder = RocksaltSurfaceBuilder(bulk, cation="Ni", anion="O")
    >>> slab = builder.create_surface(
    ...     miller_index=(0, 0, 1),
    ...     layers=6
    ... )
    """

    def __init__(
        self,
        bulk: Union[BulkStructure, Atoms],
        cation: Optional[str] = None,
        anion: str = "O",
    ):
        """
        Initialize RocksaltSurfaceBuilder.

        Parameters
        ----------
        bulk : BulkStructure or Atoms
            Rocksalt bulk structure
        cation : str, optional
            Cation species (e.g., "Ni", "Mg"). Auto-detected if None.
        anion : str
            Anion species (default "O")
        """
        super().__init__(bulk)
        self.anion = anion

        if cation is None:
            symbols = list(set(self.bulk.get_chemical_symbols()))
            cation = [s for s in symbols if s != anion][0] if len(symbols) > 1 else None

        self.cation = cation

    def get_termination_options(
        self,
        miller_index: Tuple[int, int, int] = (0, 0, 1),
    ) -> List[str]:
        """Get available termination names."""
        if miller_index == (0, 0, 1) or miller_index == (1, 1, 0):
            return [f"{self.cation}{self.anion}"]  # Non-polar, mixed
        elif miller_index == (1, 1, 1):
            return [self.cation, self.anion]  # Polar
        else:
            return ["term_0"]

    def create_surface(
        self,
        miller_index: Tuple[int, int, int],
        termination: Optional[str] = None,
        layers: int = 6,
        vacuum: float = 15.0,
        symmetric: bool = False,
        fix_bottom: int = 2,
        supercell: Optional[Tuple[int, int]] = None,
    ) -> SlabStructure:
        """
        Create rocksalt surface.

        Parameters
        ----------
        miller_index : tuple
            Miller indices. (001) and (110) are non-polar, (111) is polar.
        termination : str, optional
            For (111): cation name (e.g., "Ni") or "O"
        layers : int
            Number of layers
        vacuum : float
            Vacuum thickness
        symmetric : bool
            Create symmetric slab (important for polar (111))
        fix_bottom : int
            Number of bottom layers to fix
        supercell : tuple, optional
            Supercell size (nx, ny)

        Returns
        -------
        SlabStructure
            Rocksalt surface slab
        """
        # Check if surface is polar
        is_polar = miller_index == (1, 1, 1)

        if is_polar and symmetric:
            slab = super().create_symmetric_slab(
                miller_index, layers, vacuum, fix_bottom
            )
        else:
            slab = super().create_surface(
                miller_index, layers, vacuum, fix_bottom=fix_bottom
            )

        slab.termination = termination or f"{self.cation}{self.anion}"

        if supercell:
            slab = slab.repeat_xy(supercell[0], supercell[1])

        # Warn about polar surfaces
        if is_polar and not symmetric:
            polarity = slab.check_polarity()
            if polarity["is_polar"]:
                print(f"Warning: {miller_index} surface is polar. "
                      f"Consider using symmetric=True or applying dipole correction.")

        return slab


class FluoriteSurfaceBuilder(SurfaceBuilder):
    """
    Specialized surface builder for fluorite (MX2) structures like CeO2, ZrO2.

    Examples
    --------
    >>> bulk = BulkStructure.from_cif("CeO2.cif")
    >>> builder = FluoriteSurfaceBuilder(bulk, cation="Ce", anion="O")
    >>> slab = builder.create_surface(
    ...     miller_index=(1, 1, 1),  # Most stable
    ...     layers=6
    ... )
    """

    def __init__(
        self,
        bulk: Union[BulkStructure, Atoms],
        cation: Optional[str] = None,
        anion: str = "O",
    ):
        """
        Initialize FluoriteSurfaceBuilder.

        Parameters
        ----------
        bulk : BulkStructure or Atoms
            Fluorite bulk structure
        cation : str, optional
            Cation species (e.g., "Ce", "Zr"). Auto-detected if None.
        anion : str
            Anion species (default "O")
        """
        super().__init__(bulk)
        self.anion = anion

        if cation is None:
            symbols = list(set(self.bulk.get_chemical_symbols()))
            cation = [s for s in symbols if s != anion][0] if len(symbols) > 1 else None

        self.cation = cation

    def get_termination_options(
        self,
        miller_index: Tuple[int, int, int] = (1, 1, 1),
    ) -> List[str]:
        """Get available termination names."""
        if miller_index == (1, 1, 1):
            return [f"{self.anion}"]  # O-terminated most stable
        elif miller_index == (1, 1, 0):
            return [f"{self.cation}{self.anion}", self.anion]
        elif miller_index == (0, 0, 1):
            return [self.cation, f"{self.anion}2"]
        else:
            return ["term_0"]

    def create_surface(
        self,
        miller_index: Tuple[int, int, int],
        termination: Optional[str] = None,
        layers: int = 6,
        vacuum: float = 15.0,
        symmetric: bool = False,
        fix_bottom: int = 2,
        supercell: Optional[Tuple[int, int]] = None,
    ) -> SlabStructure:
        """
        Create fluorite surface.

        Parameters
        ----------
        miller_index : tuple
            Miller indices. (111) is most stable and non-polar.
            (110) and (100) are polar.
        termination : str, optional
            Desired termination
        layers : int
            Number of layers
        vacuum : float
            Vacuum thickness
        symmetric : bool
            Create symmetric slab
        fix_bottom : int
            Number of bottom layers to fix
        supercell : tuple, optional
            Supercell size (nx, ny)

        Returns
        -------
        SlabStructure
            Fluorite surface slab
        """
        # (111) is non-polar for fluorite, others are polar
        is_polar = miller_index != (1, 1, 1)

        if is_polar and symmetric:
            slab = super().create_symmetric_slab(
                miller_index, layers, vacuum, fix_bottom
            )
        else:
            slab = super().create_surface(
                miller_index, layers, vacuum, fix_bottom=fix_bottom
            )

        slab.termination = termination or self.anion

        if supercell:
            slab = slab.repeat_xy(supercell[0], supercell[1])

        return slab


def create_slab_from_cif(
    cif_file: Union[str, Path],
    miller_index: Tuple[int, int, int],
    layers: int = 6,
    vacuum: float = 15.0,
    supercell: Optional[Tuple[int, int]] = None,
    fix_bottom: int = 2,
) -> SlabStructure:
    """
    Convenience function to create a slab directly from a CIF file.

    Parameters
    ----------
    cif_file : str or Path
        Path to CIF file
    miller_index : tuple
        Miller indices (h, k, l)
    layers : int
        Number of layers
    vacuum : float
        Vacuum thickness
    supercell : tuple, optional
        Supercell size (nx, ny)
    fix_bottom : int
        Number of bottom layers to fix

    Returns
    -------
    SlabStructure
        Surface slab with constraints
    """
    # Load bulk
    bulk = BulkStructure.from_cif(cif_file)

    # Create surface
    builder = SurfaceBuilder(bulk)
    slab = builder.create_surface(
        miller_index=miller_index,
        layers=layers,
        vacuum=vacuum
    )

    # Apply supercell if specified
    if supercell is not None:
        slab = slab.repeat_xy(supercell[0], supercell[1])

    # Fix bottom layers
    if fix_bottom > 0:
        slab.fix_bottom_layers(fix_bottom)

    return slab
