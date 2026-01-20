"""Base classes for NH3-SOFC package."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, List, Dict, Any, Union, Tuple
import numpy as np
from ase import Atoms
from ase.io import read, write


class BaseStructure(ABC):
    """Abstract base class for structure objects."""

    def __init__(self, atoms: Optional[Atoms] = None):
        """
        Initialize base structure.

        Parameters
        ----------
        atoms : Atoms, optional
            ASE Atoms object
        """
        self._atoms = atoms

    @property
    def atoms(self) -> Atoms:
        """Get the ASE Atoms object."""
        if self._atoms is None:
            raise ValueError("No atoms object available")
        return self._atoms

    @atoms.setter
    def atoms(self, value: Atoms):
        """Set the ASE Atoms object."""
        self._atoms = value

    @property
    def n_atoms(self) -> int:
        """Number of atoms."""
        return len(self.atoms)

    @property
    def symbols(self) -> List[str]:
        """List of chemical symbols."""
        return list(self.atoms.get_chemical_symbols())

    @property
    def positions(self) -> np.ndarray:
        """Atomic positions."""
        return self.atoms.get_positions()

    @property
    def cell(self) -> np.ndarray:
        """Unit cell."""
        return np.array(self.atoms.get_cell())

    @property
    def formula(self) -> str:
        """Chemical formula."""
        return self.atoms.get_chemical_formula()

    def copy(self) -> "BaseStructure":
        """Return a copy of this structure."""
        new_obj = self.__class__.__new__(self.__class__)
        new_obj._atoms = self.atoms.copy()
        return new_obj

    def write(self, filename: Union[str, Path], format: Optional[str] = None) -> None:
        """
        Write structure to file.

        Parameters
        ----------
        filename : str or Path
            Output filename
        format : str, optional
            File format (auto-detected if not specified)
        """
        write(str(filename), self.atoms, format=format)

    @classmethod
    def from_file(cls, filename: Union[str, Path], format: Optional[str] = None) -> "BaseStructure":
        """
        Read structure from file.

        Parameters
        ----------
        filename : str or Path
            Input filename
        format : str, optional
            File format (auto-detected if not specified)

        Returns
        -------
        BaseStructure
            New structure object
        """
        atoms = read(str(filename), format=format)
        return cls(atoms)

    def get_element_indices(self, element: str) -> List[int]:
        """
        Get indices of atoms of a specific element.

        Parameters
        ----------
        element : str
            Chemical symbol

        Returns
        -------
        list
            Indices of matching atoms
        """
        return [i for i, s in enumerate(self.symbols) if s == element]

    def get_unique_elements(self) -> List[str]:
        """Get list of unique elements in structure."""
        return list(set(self.symbols))

    def get_element_count(self) -> Dict[str, int]:
        """Get count of each element."""
        counts = {}
        for symbol in self.symbols:
            counts[symbol] = counts.get(symbol, 0) + 1
        return counts


class BaseWorkflow(ABC):
    """Abstract base class for workflows."""

    def __init__(
        self,
        atoms: Atoms,
        work_dir: Union[str, Path],
        calculator: str = "vasp",
        **kwargs
    ):
        """
        Initialize workflow.

        Parameters
        ----------
        atoms : Atoms
            Input structure
        work_dir : str or Path
            Working directory
        calculator : str
            Calculator type ('vasp' or 'mace')
        **kwargs : dict
            Additional parameters
        """
        self.atoms = atoms.copy()
        self.work_dir = Path(work_dir)
        self.calculator = calculator
        self.params = kwargs
        self.results = {}

        # Create working directory
        self.work_dir.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def setup(self) -> None:
        """Set up the workflow (generate input files, etc.)."""
        pass

    @abstractmethod
    def run(self) -> Dict[str, Any]:
        """Run the workflow and return results."""
        pass

    def get_results(self) -> Dict[str, Any]:
        """Get workflow results."""
        return self.results


class BaseCalculator(ABC):
    """Abstract base class for calculator interfaces."""

    def __init__(self, atoms: Atoms, work_dir: Union[str, Path], **kwargs):
        """
        Initialize calculator.

        Parameters
        ----------
        atoms : Atoms
            Input structure
        work_dir : str or Path
            Working directory
        **kwargs : dict
            Calculator-specific parameters
        """
        self.atoms = atoms.copy()
        self.work_dir = Path(work_dir)
        self.params = kwargs

        # Create working directory
        self.work_dir.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def generate_inputs(self) -> None:
        """Generate input files for the calculation."""
        pass

    @abstractmethod
    def parse_outputs(self) -> Dict[str, Any]:
        """Parse output files and return results."""
        pass


def get_surface_atoms(
    atoms: Atoms,
    z_threshold: float = 0.2
) -> Tuple[List[int], float]:
    """
    Identify surface atoms (atoms in the top region of the slab).

    Parameters
    ----------
    atoms : Atoms
        Slab structure
    z_threshold : float
        Fraction of slab height to consider as surface (0.2 = top 20%)

    Returns
    -------
    tuple
        (list of surface atom indices, z_cutoff value)
    """
    z_positions = atoms.positions[:, 2]
    z_min, z_max = z_positions.min(), z_positions.max()
    z_range = z_max - z_min
    z_cutoff = z_max - z_threshold * z_range

    surface_indices = [i for i, z in enumerate(z_positions) if z > z_cutoff]

    return surface_indices, z_cutoff


def get_bottom_atoms(
    atoms: Atoms,
    z_threshold: float = 0.2
) -> Tuple[List[int], float]:
    """
    Identify bottom atoms (atoms in the bottom region of the slab).

    Parameters
    ----------
    atoms : Atoms
        Slab structure
    z_threshold : float
        Fraction of slab height to consider as bottom (0.2 = bottom 20%)

    Returns
    -------
    tuple
        (list of bottom atom indices, z_cutoff value)
    """
    z_positions = atoms.positions[:, 2]
    z_min, z_max = z_positions.min(), z_positions.max()
    z_range = z_max - z_min
    z_cutoff = z_min + z_threshold * z_range

    bottom_indices = [i for i, z in enumerate(z_positions) if z < z_cutoff]

    return bottom_indices, z_cutoff


def calculate_rmsd(atoms1: Atoms, atoms2: Atoms) -> float:
    """
    Calculate RMSD between two structures.

    Parameters
    ----------
    atoms1, atoms2 : Atoms
        Structures to compare

    Returns
    -------
    float
        RMSD value in Angstroms
    """
    if len(atoms1) != len(atoms2):
        raise ValueError("Structures must have same number of atoms")

    pos1 = atoms1.get_positions()
    pos2 = atoms2.get_positions()

    diff = pos1 - pos2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    return rmsd
