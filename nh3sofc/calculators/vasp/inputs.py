"""VASP input file generation.

Generates INCAR, KPOINTS, POTCAR, and POSCAR files for VASP calculations.
Supports different calculation types: relax, static, NEB, MD, frequency.
"""

import os
from pathlib import Path
from typing import Optional, List, Union, Dict, Any
import numpy as np
from ase import Atoms
from ase.io import write as ase_write

from ...core.constants import (
    DEFAULT_VASP_PARAMS,
    HUBBARD_U,
    VDW_METHODS,
)


class VASPInputGenerator:
    """
    Generator for VASP input files.

    Examples
    --------
    >>> from nh3sofc.calculators.vasp import VASPInputGenerator
    >>>
    >>> vasp = VASPInputGenerator(atoms, calc_type="relax", work_dir="./calc")
    >>> vasp.generate_all(encut=520, hubbard_u={"V": 3.25}, vdw="D3BJ")
    """

    # Preset INCAR parameters for different calculation types
    CALC_PRESETS = {
        "static": {
            "IBRION": -1,
            "NSW": 0,
            "NELM": 200,
        },
        "relax": {
            "IBRION": 2,
            "NSW": 300,
            "ISIF": 2,  # Relax ions only (for slabs)
            "EDIFFG": -0.02,
            "NELM": 200,
        },
        "relax_cell": {
            "IBRION": 2,
            "NSW": 300,
            "ISIF": 3,  # Relax ions and cell
            "EDIFFG": -0.02,
            "NELM": 200,
        },
        "neb": {
            "IBRION": 3,
            "POTIM": 0,
            "ICHAIN": 0,
            "IMAGES": 7,
            "SPRING": -5,
            "LCLIMB": True,
            "NSW": 200,
            "EDIFFG": -0.05,
        },
        "frequency": {
            "IBRION": 5,
            "NSW": 1,
            "NFREE": 2,
            "POTIM": 0.015,
            "EDIFF": 1e-7,
        },
        "md_nvt": {
            "IBRION": 0,
            "NSW": 10000,
            "POTIM": 1.0,  # 1 fs timestep
            "SMASS": 0,    # Nose-Hoover thermostat
            "TEBEG": 673,  # Starting temperature
            "TEEND": 673,  # Ending temperature
            "ISIF": 2,
            "NELM": 200,
        },
        "md_nve": {
            "IBRION": 0,
            "NSW": 10000,
            "POTIM": 1.0,
            "SMASS": -1,   # NVE ensemble
            "TEBEG": 673,
            "ISIF": 2,
        },
    }

    def __init__(
        self,
        atoms: Atoms,
        calc_type: str = "relax",
        work_dir: Union[str, Path] = "./",
        **kwargs,
    ):
        """
        Initialize VASPInputGenerator.

        Parameters
        ----------
        atoms : Atoms
            Input structure
        calc_type : str
            Calculation type: "static", "relax", "relax_cell", "neb",
            "frequency", "md_nvt", "md_nve"
        work_dir : str or Path
            Working directory for output files
        **kwargs : dict
            Additional parameters to override defaults
        """
        self.atoms = atoms.copy()
        self.calc_type = calc_type
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # Initialize INCAR parameters
        self.incar_params = DEFAULT_VASP_PARAMS.copy()

        # Apply preset for calculation type
        if calc_type in self.CALC_PRESETS:
            self.incar_params.update(self.CALC_PRESETS[calc_type])

        # Apply user overrides
        self.incar_params.update(kwargs)

    def sort_atoms_by_element(self) -> None:
        """
        Sort atoms by chemical symbol (alphabetically).

        This ensures proper POSCAR format where atoms are grouped by element
        (e.g., "La La La V V V O O O") instead of scattered ordering
        (e.g., "La V O La V O") which some visualization tools can't read.

        Also ensures consistency between POSCAR, POTCAR, and INCAR.
        """
        symbols = self.atoms.get_chemical_symbols()
        unique_symbols = sorted(set(symbols))
        sort_indices = []
        for sym in unique_symbols:
            for i, s in enumerate(symbols):
                if s == sym:
                    sort_indices.append(i)
        self.atoms = self.atoms[sort_indices]

    def set_encut(self, encut: float) -> None:
        """Set plane-wave cutoff energy."""
        self.incar_params["ENCUT"] = encut

    def set_kpoints(
        self,
        kspacing: Optional[float] = None,
        kpoints: Optional[List[int]] = None,
    ) -> None:
        """
        Set k-point parameters.

        Parameters
        ----------
        kspacing : float, optional
            K-point spacing (generates automatic mesh)
        kpoints : list, optional
            Explicit k-point mesh [kx, ky, kz]
        """
        if kspacing is not None:
            self.incar_params["KSPACING"] = kspacing
            self._kpoints_mesh = None
        elif kpoints is not None:
            self._kpoints_mesh = kpoints
            if "KSPACING" in self.incar_params:
                del self.incar_params["KSPACING"]

    def set_hubbard_u(self, u_values: Dict[str, float]) -> None:
        """
        Set Hubbard U parameters.

        Parameters
        ----------
        u_values : dict
            Dictionary of {element: U_value}
        """
        # Get unique elements in order
        symbols = self.atoms.get_chemical_symbols()
        unique_elements = []
        for s in symbols:
            if s not in unique_elements:
                unique_elements.append(s)

        # Build LDAU parameters
        ldaul = []
        ldauu = []
        ldauj = []

        for elem in unique_elements:
            if elem in u_values and u_values[elem] > 0:
                ldaul.append(2)  # d-orbitals
                ldauu.append(u_values[elem])
                ldauj.append(0)
            else:
                ldaul.append(-1)
                ldauu.append(0)
                ldauj.append(0)

        self.incar_params["LDAU"] = True
        self.incar_params["LDAUTYPE"] = 2
        self.incar_params["LDAUL"] = " ".join(map(str, ldaul))
        self.incar_params["LDAUU"] = " ".join(map(str, ldauu))
        self.incar_params["LDAUJ"] = " ".join(map(str, ldauj))
        self.incar_params["LDAUPRINT"] = 1
        self.incar_params["LMAXMIX"] = 4

    def set_vdw(self, method: str) -> None:
        """
        Set van der Waals correction.

        Parameters
        ----------
        method : str
            VdW method: "D3", "D3BJ", "D4", "TS", "TSHIRSH"
        """
        if method.upper() in VDW_METHODS:
            self.incar_params.update(VDW_METHODS[method.upper()])
        else:
            raise ValueError(f"Unknown vdW method: {method}")

    def set_spin_polarized(self, spin: bool = True, magmom: Optional[List[float]] = None) -> None:
        """
        Set spin polarization.

        Parameters
        ----------
        spin : bool
            Enable spin polarization
        magmom : list, optional
            Initial magnetic moments per atom
        """
        self.incar_params["ISPIN"] = 2 if spin else 1

        if magmom is not None:
            self.incar_params["MAGMOM"] = " ".join(map(str, magmom))
        elif spin:
            # Default: set small moment on transition metals
            tm = ["V", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", "Cr"]
            magmom = []
            for symbol in self.atoms.get_chemical_symbols():
                if symbol in tm:
                    magmom.append(2.0)
                else:
                    magmom.append(0.0)
            if any(m > 0 for m in magmom):
                self.incar_params["MAGMOM"] = " ".join(map(str, magmom))

    def set_surface_settings(self) -> None:
        """Apply settings optimized for surface calculations."""
        self.incar_params["ISIF"] = 2  # Don't relax cell
        self.incar_params["IDIPOL"] = 3  # Dipole correction in z
        self.incar_params["LDIPOL"] = True

    def generate_incar(self) -> str:
        """
        Generate INCAR file content.

        Returns
        -------
        str
            INCAR file content
        """
        lines = [f"# VASP INCAR - {self.calc_type} calculation"]
        lines.append(f"# Generated by nh3sofc")
        lines.append("")

        # Group parameters by category
        categories = {
            "Electronic": ["ENCUT", "EDIFF", "NELM", "ISMEAR", "SIGMA", "PREC",
                          "ALGO", "LREAL", "ISPIN", "MAGMOM"],
            "Ionic": ["IBRION", "NSW", "ISIF", "EDIFFG", "POTIM", "NFREE"],
            "DFT+U": ["LDAU", "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ",
                     "LDAUPRINT", "LMAXMIX"],
            "vdW": ["IVDW"],
            "Output": ["LWAVE", "LCHARG", "LORBIT", "NCORE", "KPAR"],
            "Dipole": ["IDIPOL", "LDIPOL"],
            "MD": ["SMASS", "TEBEG", "TEEND"],
            "NEB": ["ICHAIN", "IMAGES", "SPRING", "LCLIMB"],
        }

        written = set()

        for category, params in categories.items():
            category_lines = []
            for param in params:
                if param in self.incar_params:
                    value = self.incar_params[param]
                    if isinstance(value, bool):
                        value = ".TRUE." if value else ".FALSE."
                    category_lines.append(f"  {param} = {value}")
                    written.add(param)

            if category_lines:
                lines.append(f"# {category}")
                lines.extend(category_lines)
                lines.append("")

        # Write any remaining parameters
        remaining = []
        for param, value in self.incar_params.items():
            if param not in written:
                if isinstance(value, bool):
                    value = ".TRUE." if value else ".FALSE."
                remaining.append(f"  {param} = {value}")

        if remaining:
            lines.append("# Other")
            lines.extend(remaining)

        return "\n".join(lines)

    def generate_kpoints(self, kspacing: Optional[float] = None) -> str:
        """
        Generate KPOINTS file content.

        Parameters
        ----------
        kspacing : float, optional
            K-point spacing for automatic mesh

        Returns
        -------
        str
            KPOINTS file content
        """
        if hasattr(self, "_kpoints_mesh") and self._kpoints_mesh is not None:
            kpts = self._kpoints_mesh
        elif kspacing is not None or "KSPACING" in self.incar_params:
            # Auto-generate from spacing
            spacing = kspacing or self.incar_params.get("KSPACING", 0.03)
            cell = self.atoms.get_cell()
            reciprocal = cell.reciprocal() * 2 * np.pi
            lengths = np.linalg.norm(reciprocal, axis=1)
            kpts = [max(1, int(np.ceil(l / spacing))) for l in lengths]
        else:
            # Default: Gamma only
            kpts = [1, 1, 1]

        lines = [
            "Automatic mesh",
            "0",
            "Gamma",  # Gamma-centered
            f"  {kpts[0]}  {kpts[1]}  {kpts[2]}",
            "  0  0  0",
        ]

        return "\n".join(lines)

    def generate_poscar(self) -> str:
        """
        Generate POSCAR file content with standard VASP format.

        Uses proper format with unique element labels only (e.g., "La N O")
        instead of repeated labels for each atom (e.g., "La O N La N O O").

        Note: If using generate_all() with sort_atoms=True (default), atoms
        are already sorted by element type before this method is called.

        Returns
        -------
        str
            POSCAR file content
        """
        from ...core.io import write_poscar
        from io import StringIO
        import tempfile
        import os

        # Write to temporary file and read back content
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vasp", delete=False) as f:
            temp_path = f.name

        try:
            write_poscar(
                self.atoms,
                temp_path,
                comment=self.atoms.get_chemical_formula(),
                direct=True,
                sort=False,  # Already sorted if sort_atoms=True in generate_all()
            )
            with open(temp_path, "r") as f:
                content = f.read()
        finally:
            os.unlink(temp_path)

        return content

    def generate_potcar(
        self,
        potcar_dir: Optional[Union[str, Path]] = None,
        potcar_map: Optional[Dict[str, str]] = None,
    ) -> Optional[str]:
        """
        Generate POTCAR by concatenating pseudopotential files.

        Parameters
        ----------
        potcar_dir : str or Path, optional
            Directory containing POTCAR files. Defaults to VASP_PP_PATH env var.
        potcar_map : dict, optional
            Mapping of element to POTCAR variant (e.g., {"O": "O_s"})

        Returns
        -------
        str or None
            POTCAR content, or None if POTCARs not found
        """
        if potcar_dir is None:
            potcar_dir = os.environ.get("VASP_PP_PATH")
            if potcar_dir is None:
                print("Warning: VASP_PP_PATH not set, skipping POTCAR generation")
                return None

        potcar_dir = Path(potcar_dir)

        if potcar_map is None:
            potcar_map = {}

        # Get unique elements in order of appearance
        symbols = self.atoms.get_chemical_symbols()
        unique_elements = []
        for s in symbols:
            if s not in unique_elements:
                unique_elements.append(s)

        # Concatenate POTCARs
        potcar_content = []

        for elem in unique_elements:
            # Use mapped name or default
            potcar_name = potcar_map.get(elem, elem)

            # Try common locations
            possible_paths = [
                potcar_dir / "potpaw_PBE" / potcar_name / "POTCAR",
                potcar_dir / "potpaw_PBE.54" / potcar_name / "POTCAR",
                potcar_dir / "PBE" / potcar_name / "POTCAR",
                potcar_dir / potcar_name / "POTCAR",
            ]

            found = False
            for path in possible_paths:
                if path.exists():
                    with open(path, "r") as f:
                        potcar_content.append(f.read())
                    found = True
                    break

            if not found:
                print(f"Warning: POTCAR not found for {elem}")
                return None

        return "".join(potcar_content)

    def generate_all(
        self,
        encut: Optional[float] = None,
        kspacing: Optional[float] = None,
        hubbard_u: Optional[Dict[str, float]] = None,
        vdw: Optional[str] = None,
        spin_polarized: bool = True,
        is_surface: bool = False,
        potcar_dir: Optional[str] = None,
        potcar_map: Optional[Dict[str, str]] = None,
        sort_atoms: bool = True,
    ) -> Dict[str, Path]:
        """
        Generate all VASP input files.

        Parameters
        ----------
        encut : float, optional
            Plane-wave cutoff
        kspacing : float, optional
            K-point spacing
        hubbard_u : dict, optional
            Hubbard U values by element
        vdw : str, optional
            VdW correction method
        spin_polarized : bool
            Enable spin polarization
        is_surface : bool
            Apply surface-specific settings
        potcar_dir : str, optional
            POTCAR directory
        potcar_map : dict, optional
            Element to POTCAR variant mapping
        sort_atoms : bool
            If True (default), sort atoms by element type before writing files.
            This ensures proper POSCAR format and consistency across all files.

        Returns
        -------
        dict
            Dictionary mapping file names to paths
        """
        # Sort atoms by element for proper POSCAR format and consistency
        if sort_atoms:
            self.sort_atoms_by_element()

        # Apply settings
        if encut is not None:
            self.set_encut(encut)
        if kspacing is not None:
            self.set_kpoints(kspacing=kspacing)
        if hubbard_u is not None:
            self.set_hubbard_u(hubbard_u)
        if vdw is not None:
            self.set_vdw(vdw)
        if spin_polarized:
            self.set_spin_polarized(True)
        if is_surface:
            self.set_surface_settings()

        # Generate files
        files = {}

        # INCAR
        incar_path = self.work_dir / "INCAR"
        with open(incar_path, "w") as f:
            f.write(self.generate_incar())
        files["INCAR"] = incar_path

        # KPOINTS
        kpoints_path = self.work_dir / "KPOINTS"
        with open(kpoints_path, "w") as f:
            f.write(self.generate_kpoints())
        files["KPOINTS"] = kpoints_path

        # POSCAR
        poscar_path = self.work_dir / "POSCAR"
        with open(poscar_path, "w") as f:
            f.write(self.generate_poscar())
        files["POSCAR"] = poscar_path

        # POTCAR (optional)
        potcar_content = self.generate_potcar(potcar_dir, potcar_map)
        if potcar_content is not None:
            potcar_path = self.work_dir / "POTCAR"
            with open(potcar_path, "w") as f:
                f.write(potcar_content)
            files["POTCAR"] = potcar_path

        return files


def create_vasp_inputs(
    atoms: Atoms,
    work_dir: Union[str, Path],
    calc_type: str = "relax",
    encut: float = 520,
    kspacing: float = 0.03,
    hubbard_u: Optional[Dict[str, float]] = None,
    vdw: Optional[str] = None,
    is_surface: bool = True,
    sort_atoms: bool = True,
) -> Dict[str, Path]:
    """
    Convenience function to create VASP inputs.

    Parameters
    ----------
    atoms : Atoms
        Structure
    work_dir : str or Path
        Output directory
    calc_type : str
        Calculation type
    encut : float
        Cutoff energy
    kspacing : float
        K-point spacing
    hubbard_u : dict, optional
        Hubbard U values
    vdw : str, optional
        VdW method
    is_surface : bool
        Surface calculation
    sort_atoms : bool
        If True (default), sort atoms by element type for proper POSCAR format

    Returns
    -------
    dict
        Generated file paths
    """
    gen = VASPInputGenerator(atoms, calc_type=calc_type, work_dir=work_dir)
    return gen.generate_all(
        encut=encut,
        kspacing=kspacing,
        hubbard_u=hubbard_u,
        vdw=vdw,
        is_surface=is_surface,
        sort_atoms=sort_atoms,
    )
