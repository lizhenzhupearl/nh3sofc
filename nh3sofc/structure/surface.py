"""Surface generation from bulk structures."""

from pathlib import Path
from typing import Optional, List, Union, Tuple
import numpy as np
from ase import Atoms
from ase.build import surface, add_vacuum
from ase.constraints import FixAtoms
from ase.io import write

from ..core.base import BaseStructure, get_bottom_atoms, get_surface_atoms
from ..core.constants import DEFAULT_SURFACE_PARAMS
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
            Repeated slab
        """
        new_atoms = self.atoms.repeat((nx, ny, 1))
        new_slab = SlabStructure(
            new_atoms,
            miller_index=self.miller_index,
            termination=self.termination,
            n_layers=self.n_layers
        )
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

        return SlabStructure(
            slab,
            miller_index=miller_index,
            termination=term_label,
            n_layers=layers
        )

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
