"""IO utilities for structure file handling.

Provides proper POSCAR and CIF writing with standard formats.
Addresses ASE's default behavior of writing repeated atom types.
"""

from pathlib import Path
from typing import Optional, Union, List, Dict, Any
import numpy as np
from ase import Atoms
from ase.io import write as ase_write


def sort_atoms_by_element(atoms: Atoms) -> Atoms:
    """
    Sort atoms by chemical symbol (alphabetically).

    Parameters
    ----------
    atoms : Atoms
        Input structure

    Returns
    -------
    Atoms
        Sorted structure with atoms grouped by element

    Examples
    --------
    >>> sorted_atoms = sort_atoms_by_element(atoms)
    """
    symbols = atoms.get_chemical_symbols()
    unique_symbols = sorted(set(symbols))
    sort_indices = []
    for sym in unique_symbols:
        for i, s in enumerate(symbols):
            if s == sym:
                sort_indices.append(i)
    return atoms[sort_indices]


def write_poscar(
    atoms: Atoms,
    filename: Union[str, Path],
    comment: Optional[str] = None,
    direct: bool = True,
    sort: bool = True,
    vasp5: bool = True,
) -> Path:
    """
    Write structure to POSCAR with standard VASP format.

    Ensures proper format with unique element labels only (e.g., "La N O")
    instead of repeated labels for each atom (e.g., "La O N La N O O").

    Parameters
    ----------
    atoms : Atoms
        Structure to write
    filename : str or Path
        Output filename
    comment : str, optional
        Comment line (default: chemical formula)
    direct : bool
        Use direct (fractional) coordinates (default: True)
    sort : bool
        Sort atoms by element type (default: True)
    vasp5 : bool
        Use VASP5 format with element symbols (default: True)

    Returns
    -------
    Path
        Path to written file

    Examples
    --------
    >>> from ase.build import bulk
    >>> atoms = bulk("NaCl", "rocksalt", a=5.64)
    >>> write_poscar(atoms, "POSCAR")
    """
    filepath = Path(filename)
    atoms = atoms.copy()

    # Sort atoms by element for proper grouping
    if sort:
        atoms = sort_atoms_by_element(atoms)

    # Get element information
    symbols = atoms.get_chemical_symbols()
    unique_elements = []
    element_counts = []

    current_element = None
    current_count = 0

    for symbol in symbols:
        if symbol != current_element:
            if current_element is not None:
                unique_elements.append(current_element)
                element_counts.append(current_count)
            current_element = symbol
            current_count = 1
        else:
            current_count += 1

    # Don't forget the last element
    if current_element is not None:
        unique_elements.append(current_element)
        element_counts.append(current_count)

    # Build POSCAR content
    lines = []

    # Line 1: Comment
    if comment is None:
        comment = atoms.get_chemical_formula()
    lines.append(comment)

    # Line 2: Scale factor
    lines.append("1.0")

    # Lines 3-5: Cell vectors
    cell = atoms.get_cell()
    for i in range(3):
        lines.append(f"  {cell[i, 0]:21.16f}  {cell[i, 1]:21.16f}  {cell[i, 2]:21.16f}")

    # Line 6: Element symbols (VASP5 format)
    if vasp5:
        lines.append("  " + "  ".join(unique_elements))

    # Line 7: Element counts
    lines.append("  " + "  ".join(str(c) for c in element_counts))

    # Line 8: Coordinate type
    if direct:
        lines.append("Direct")
        positions = atoms.get_scaled_positions()
    else:
        lines.append("Cartesian")
        positions = atoms.get_positions()

    # Remaining lines: Coordinates
    for pos in positions:
        lines.append(f"  {pos[0]:19.16f}  {pos[1]:19.16f}  {pos[2]:19.16f}")

    # Write to file
    with open(filepath, "w") as f:
        f.write("\n".join(lines) + "\n")

    return filepath


def write_cif(
    atoms: Atoms,
    filename: Union[str, Path],
) -> Path:
    """
    Write structure to CIF format.

    Parameters
    ----------
    atoms : Atoms
        Structure to write
    filename : str or Path
        Output filename

    Returns
    -------
    Path
        Path to written file

    Examples
    --------
    >>> write_cif(atoms, "structure.cif")
    """
    filepath = Path(filename)
    ase_write(str(filepath), atoms, format="cif")
    return filepath


def save_structure(
    atoms: Atoms,
    work_dir: Union[str, Path],
    name: str,
    formats: Optional[List[str]] = None,
    comment: Optional[str] = None,
) -> Dict[str, Path]:
    """
    Save structure to work directory in multiple formats.

    Parameters
    ----------
    atoms : Atoms
        Structure to save
    work_dir : str or Path
        Directory to save files
    name : str
        Base name for files (without extension)
    formats : list of str, optional
        Formats to save: "poscar", "cif", or both (default: ["poscar"])
    comment : str, optional
        Comment for POSCAR file

    Returns
    -------
    dict
        Dictionary mapping format to file path

    Examples
    --------
    >>> paths = save_structure(atoms, "./configs", "relaxed", formats=["poscar", "cif"])
    >>> print(paths["poscar"])
    ./configs/relaxed.vasp
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    if formats is None:
        formats = ["poscar"]

    saved_files = {}

    for fmt in formats:
        fmt = fmt.lower()
        if fmt == "poscar":
            filepath = work_dir / f"{name}.vasp"
            write_poscar(atoms, filepath, comment=comment)
            saved_files["poscar"] = filepath
        elif fmt == "cif":
            filepath = work_dir / f"{name}.cif"
            write_cif(atoms, filepath)
            saved_files["cif"] = filepath
        else:
            raise ValueError(f"Unsupported format: {fmt}. Use 'poscar' or 'cif'.")

    return saved_files


def save_configurations(
    configurations: List[Atoms],
    work_dir: Union[str, Path],
    name_prefix: str = "config",
    formats: Optional[List[str]] = None,
    names: Optional[List[str]] = None,
) -> List[Dict[str, Path]]:
    """
    Save multiple configurations to work directory.

    Parameters
    ----------
    configurations : list of Atoms
        List of structures to save
    work_dir : str or Path
        Directory to save files
    name_prefix : str
        Prefix for file names (default: "config")
    formats : list of str, optional
        Formats to save: "poscar", "cif", or both (default: ["poscar"])
    names : list of str, optional
        Custom names for each configuration (overrides prefix)

    Returns
    -------
    list of dict
        List of dictionaries mapping format to file path

    Examples
    --------
    >>> configs = [atoms1, atoms2, atoms3]
    >>> paths = save_configurations(configs, "./decomposition", name_prefix="step")
    >>> # Creates: step_001.vasp, step_002.vasp, step_003.vasp

    >>> # With custom names
    >>> paths = save_configurations(
    ...     configs, "./intermediates",
    ...     names=["initial", "TS", "final"],
    ...     formats=["poscar", "cif"]
    ... )
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    if formats is None:
        formats = ["poscar"]

    all_paths = []

    for i, atoms in enumerate(configurations):
        if names is not None and i < len(names):
            name = names[i]
        else:
            name = f"{name_prefix}_{i + 1:03d}"

        paths = save_structure(atoms, work_dir, name, formats=formats)
        all_paths.append(paths)

    return all_paths
