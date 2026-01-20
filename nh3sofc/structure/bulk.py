"""Bulk structure loading and manipulation."""

from pathlib import Path
from typing import Optional, List, Union, Tuple
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.build import make_supercell

from ..core.base import BaseStructure


class BulkStructure(BaseStructure):
    """
    Class for handling bulk crystal structures.

    Provides methods for loading structures from files, creating supercells,
    and manipulating bulk structures for surface generation.

    Examples
    --------
    >>> bulk = BulkStructure.from_cif("LaVO3.cif")
    >>> supercell = bulk.make_supercell([2, 2, 2])
    >>> print(supercell.formula)
    La8V8O24
    """

    def __init__(self, atoms: Atoms):
        """
        Initialize BulkStructure.

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object representing the bulk structure
        """
        super().__init__(atoms)

        # Ensure periodic boundary conditions
        self.atoms.pbc = True

    @classmethod
    def from_cif(cls, filename: Union[str, Path]) -> "BulkStructure":
        """
        Read bulk structure from CIF file.

        Parameters
        ----------
        filename : str or Path
            Path to CIF file

        Returns
        -------
        BulkStructure
            New BulkStructure object
        """
        atoms = read(str(filename), format="cif")
        return cls(atoms)

    @classmethod
    def from_poscar(cls, filename: Union[str, Path]) -> "BulkStructure":
        """
        Read bulk structure from VASP POSCAR file.

        Parameters
        ----------
        filename : str or Path
            Path to POSCAR file

        Returns
        -------
        BulkStructure
            New BulkStructure object
        """
        atoms = read(str(filename), format="vasp")
        return cls(atoms)

    @classmethod
    def from_atoms(cls, atoms: Atoms) -> "BulkStructure":
        """
        Create BulkStructure from ASE Atoms object.

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object

        Returns
        -------
        BulkStructure
            New BulkStructure object
        """
        return cls(atoms.copy())

    def make_supercell(
        self,
        size: Union[List[int], Tuple[int, int, int], np.ndarray]
    ) -> "BulkStructure":
        """
        Create a supercell of the bulk structure.

        Parameters
        ----------
        size : list or tuple or array
            Supercell size [nx, ny, nz] or transformation matrix

        Returns
        -------
        BulkStructure
            New BulkStructure with supercell

        Examples
        --------
        >>> supercell = bulk.make_supercell([2, 2, 2])
        >>> # Or with transformation matrix
        >>> supercell = bulk.make_supercell([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        """
        size = np.array(size)

        # If 1D array, convert to diagonal matrix
        if size.ndim == 1:
            size = np.diag(size)

        supercell_atoms = make_supercell(self.atoms, size)
        return BulkStructure(supercell_atoms)

    def repeat(self, rep: Union[int, Tuple[int, int, int]]) -> "BulkStructure":
        """
        Repeat the structure along each axis.

        Parameters
        ----------
        rep : int or tuple
            Number of repetitions along each axis

        Returns
        -------
        BulkStructure
            New repeated structure
        """
        new_atoms = self.atoms.repeat(rep)
        return BulkStructure(new_atoms)

    def get_lattice_parameters(self) -> dict:
        """
        Get lattice parameters (a, b, c, alpha, beta, gamma).

        Returns
        -------
        dict
            Lattice parameters
        """
        cell = self.atoms.get_cell()
        lengths = cell.lengths()
        angles = cell.angles()

        return {
            "a": lengths[0],
            "b": lengths[1],
            "c": lengths[2],
            "alpha": angles[0],
            "beta": angles[1],
            "gamma": angles[2],
            "volume": cell.volume,
        }

    def get_spacegroup(self, symprec: float = 1e-3) -> dict:
        """
        Get spacegroup information using spglib.

        Parameters
        ----------
        symprec : float
            Symmetry precision for spglib

        Returns
        -------
        dict
            Spacegroup information (number, symbol, etc.)

        Raises
        ------
        ImportError
            If spglib is not installed
        """
        try:
            import spglib
        except ImportError:
            raise ImportError("spglib is required for spacegroup analysis. "
                            "Install with: pip install spglib")

        cell = (
            self.atoms.get_cell(),
            self.atoms.get_scaled_positions(),
            self.atoms.get_atomic_numbers()
        )

        spg = spglib.get_spacegroup(cell, symprec=symprec)
        spg_dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)

        if spg_dataset is None:
            return {"symbol": "Unknown", "number": 0}

        return {
            "symbol": spg,
            "number": spg_dataset["number"],
            "international": spg_dataset["international"],
            "hall": spg_dataset["hall"],
            "pointgroup": spg_dataset["pointgroup"],
        }

    def standardize(self, symprec: float = 1e-3) -> "BulkStructure":
        """
        Standardize the structure using spglib.

        Parameters
        ----------
        symprec : float
            Symmetry precision

        Returns
        -------
        BulkStructure
            Standardized structure
        """
        try:
            import spglib
        except ImportError:
            raise ImportError("spglib is required. Install with: pip install spglib")

        cell = (
            self.atoms.get_cell(),
            self.atoms.get_scaled_positions(),
            self.atoms.get_atomic_numbers()
        )

        std_cell = spglib.standardize_cell(cell, symprec=symprec)

        if std_cell is None:
            raise ValueError("Could not standardize cell")

        lattice, positions, numbers = std_cell

        from ase.data import chemical_symbols
        symbols = [chemical_symbols[n] for n in numbers]

        new_atoms = Atoms(
            symbols=symbols,
            scaled_positions=positions,
            cell=lattice,
            pbc=True
        )

        return BulkStructure(new_atoms)

    def get_primitive(self, symprec: float = 1e-3) -> "BulkStructure":
        """
        Find primitive cell using spglib.

        Parameters
        ----------
        symprec : float
            Symmetry precision

        Returns
        -------
        BulkStructure
            Primitive cell structure
        """
        try:
            import spglib
        except ImportError:
            raise ImportError("spglib is required. Install with: pip install spglib")

        cell = (
            self.atoms.get_cell(),
            self.atoms.get_scaled_positions(),
            self.atoms.get_atomic_numbers()
        )

        prim_cell = spglib.find_primitive(cell, symprec=symprec)

        if prim_cell is None:
            raise ValueError("Could not find primitive cell")

        lattice, positions, numbers = prim_cell

        from ase.data import chemical_symbols
        symbols = [chemical_symbols[n] for n in numbers]

        new_atoms = Atoms(
            symbols=symbols,
            scaled_positions=positions,
            cell=lattice,
            pbc=True
        )

        return BulkStructure(new_atoms)

    def center(self, axis: Optional[int] = None) -> "BulkStructure":
        """
        Center the structure in the cell.

        Parameters
        ----------
        axis : int, optional
            Axis to center (0, 1, 2). If None, center all axes.

        Returns
        -------
        BulkStructure
            Centered structure
        """
        new_atoms = self.atoms.copy()
        if axis is None:
            new_atoms.center()
        else:
            new_atoms.center(axis=axis)
        return BulkStructure(new_atoms)

    def to_poscar(self, filename: Union[str, Path]) -> None:
        """
        Write structure to VASP POSCAR format.

        Parameters
        ----------
        filename : str or Path
            Output filename
        """
        write(str(filename), self.atoms, format="vasp", direct=True)

    def to_cif(self, filename: Union[str, Path]) -> None:
        """
        Write structure to CIF format.

        Parameters
        ----------
        filename : str or Path
            Output filename
        """
        write(str(filename), self.atoms, format="cif")

    def copy(self) -> "BulkStructure":
        """Return a copy of this structure."""
        return BulkStructure(self.atoms.copy())

    def __repr__(self) -> str:
        params = self.get_lattice_parameters()
        return (
            f"BulkStructure({self.formula}, "
            f"a={params['a']:.3f}, b={params['b']:.3f}, c={params['c']:.3f})"
        )
