"""IO utilities for structure file handling.

Provides proper POSCAR and CIF writing with standard formats.
Addresses ASE's default behavior of writing repeated atom types.
"""

from pathlib import Path
from typing import Optional, Union, List, Dict, Any
from datetime import datetime
import numpy as np
from ase import Atoms
from ase.io import write as ase_write


def generate_work_dir(
    atoms: Atoms,
    base_dir: str = "work",
    prefix: Optional[str] = None,
    include_timestamp: bool = False,
) -> Path:
    """
    Generate a meaningful work directory name based on structure properties.

    Parameters
    ----------
    atoms : Atoms
        Structure to analyze for naming
    base_dir : str
        Base directory (default: "work")
    prefix : str, optional
        Optional prefix to add (e.g., "surface", "bulk", "adsorbate")
    include_timestamp : bool
        Include timestamp in directory name (default: False)

    Returns
    -------
    Path
        Generated work directory path

    Examples
    --------
    >>> atoms = bulk("LaVO3", "perovskite", a=3.9)
    >>> generate_work_dir(atoms)
    Path('work/LaVO3_5atoms')

    >>> generate_work_dir(surface_atoms, prefix="surface")
    Path('work/surface_La8V8O24_40atoms')
    """
    # Get formula and atom count
    formula = atoms.get_chemical_formula(mode='metal')
    n_atoms = len(atoms)

    # Build directory name components
    components = []

    if prefix:
        components.append(prefix)

    components.append(formula)
    components.append(f"{n_atoms}atoms")

    # Check if it's a slab (has vacuum in z-direction)
    if atoms.pbc.any():
        cell = atoms.get_cell()
        positions = atoms.get_positions()
        z_range = positions[:, 2].max() - positions[:, 2].min()
        cell_z = cell[2, 2]
        if cell_z > z_range + 5.0:  # More than 5 Å vacuum
            if "surface" not in (prefix or ""):
                components.insert(0, "slab")

    if include_timestamp:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        components.append(timestamp)

    dir_name = "_".join(components)
    work_dir = Path(base_dir) / dir_name

    return work_dir


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
    work_dir: Optional[Union[str, Path]] = None,
    name: str = "structure",
    formats: Optional[List[str]] = None,
    comment: Optional[str] = None,
    prefix: Optional[str] = None,
) -> Dict[str, Path]:
    """
    Save structure to work directory in multiple formats.

    Parameters
    ----------
    atoms : Atoms
        Structure to save
    work_dir : str or Path, optional
        Directory to save files. If not provided, automatically generates
        a meaningful directory name based on the structure (e.g., "work/LaVO3_40atoms")
    name : str
        Base name for files (default: "structure")
    formats : list of str, optional
        Formats to save: "poscar", "cif", or both (default: ["poscar"])
    comment : str, optional
        Comment for POSCAR file
    prefix : str, optional
        Prefix for auto-generated work_dir (e.g., "surface", "bulk")

    Returns
    -------
    dict
        Dictionary mapping format to file path, plus "work_dir" key

    Examples
    --------
    >>> # Auto-generate work directory
    >>> paths = save_structure(atoms, name="relaxed")
    >>> print(paths["work_dir"])
    work/LaVO3_40atoms
    >>> print(paths["poscar"])
    work/LaVO3_40atoms/relaxed.vasp

    >>> # Specify work directory explicitly
    >>> paths = save_structure(atoms, "./my_configs", "relaxed", formats=["poscar", "cif"])
    >>> print(paths["poscar"])
    ./my_configs/relaxed.vasp
    """
    # Auto-generate work_dir if not provided
    if work_dir is None:
        work_dir = generate_work_dir(atoms, prefix=prefix)
    else:
        work_dir = Path(work_dir)

    work_dir.mkdir(parents=True, exist_ok=True)

    if formats is None:
        formats = ["poscar"]

    saved_files = {"work_dir": work_dir}

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
    work_dir: Optional[Union[str, Path]] = None,
    name_prefix: str = "config",
    formats: Optional[List[str]] = None,
    names: Optional[List[str]] = None,
    prefix: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Save multiple configurations to work directory.

    Parameters
    ----------
    configurations : list of Atoms
        List of structures to save
    work_dir : str or Path, optional
        Directory to save files. If not provided, automatically generates
        a meaningful directory name based on the first structure
    name_prefix : str
        Prefix for file names (default: "config")
    formats : list of str, optional
        Formats to save: "poscar", "cif", or both (default: ["poscar"])
    names : list of str, optional
        Custom names for each configuration (overrides prefix)
    prefix : str, optional
        Prefix for auto-generated work_dir (e.g., "decomposition", "screening")

    Returns
    -------
    dict
        Dictionary with "work_dir" and "configs" (list of path dicts)

    Examples
    --------
    >>> # Auto-generate work directory from first structure
    >>> configs = [atoms1, atoms2, atoms3]
    >>> result = save_configurations(configs, name_prefix="step")
    >>> print(result["work_dir"])
    work/LaVO3_40atoms
    >>> # Creates: work/LaVO3_40atoms/step_001.vasp, step_002.vasp, step_003.vasp

    >>> # Specify work directory explicitly
    >>> result = save_configurations(configs, "./decomposition", name_prefix="step")

    >>> # With custom names
    >>> result = save_configurations(
    ...     configs,
    ...     names=["initial", "TS", "final"],
    ...     formats=["poscar", "cif"]
    ... )
    """
    if not configurations:
        raise ValueError("No configurations provided")

    # Auto-generate work_dir if not provided (use first config for naming)
    if work_dir is None:
        work_dir = generate_work_dir(configurations[0], prefix=prefix)
    else:
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
        # Remove work_dir from individual paths since we return it at top level
        paths.pop("work_dir", None)
        all_paths.append(paths)

    return {"work_dir": work_dir, "configs": all_paths}
