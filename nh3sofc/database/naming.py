"""Standardized naming conventions for calculations.

Provides consistent directory and file naming for organizing
DFT calculations and results.
"""

from pathlib import Path
from typing import Optional, Dict, Any, List
import re


class NamingConvention:
    """
    Standardized naming for NH3-SOFC calculations.

    Format: {material}_{miller}_{termination}_{calc_type}_{suffix}

    Examples
    --------
    >>> namer = NamingConvention("LaVO3", miller=(0,0,1), termination="LaO")
    >>> namer.get_slab_name()
    'LaVO3_001_LaO'
    >>> namer.get_calc_dir("relax", suffix="NH3")
    'LaVO3_001_LaO_relax_NH3'
    """

    def __init__(
        self,
        material: str,
        miller: tuple = (0, 0, 1),
        termination: Optional[str] = None,
        vacancy_concentration: float = 0.0,
        nitrogen_fraction: float = 0.0,
    ):
        """
        Initialize NamingConvention.

        Parameters
        ----------
        material : str
            Base material name (e.g., "LaVO3")
        miller : tuple
            Miller indices (h, k, l)
        termination : str, optional
            Surface termination (e.g., "LaO", "VO2")
        vacancy_concentration : float
            Oxygen vacancy concentration (0-1)
        nitrogen_fraction : float
            Nitrogen substitution fraction (0-1)
        """
        self.material = material
        self.miller = miller
        self.termination = termination
        self.vacancy_concentration = vacancy_concentration
        self.nitrogen_fraction = nitrogen_fraction

    @property
    def miller_str(self) -> str:
        """Get Miller indices as string (e.g., '001')."""
        return "".join(str(abs(i)) for i in self.miller)

    @property
    def material_str(self) -> str:
        """Get material string with modifications."""
        name = self.material

        # Add nitrogen fraction if present
        if self.nitrogen_fraction > 0:
            n_frac = f"N{self.nitrogen_fraction:.2f}".replace(".", "")
            name = f"{name}_{n_frac}"

        # Add vacancy concentration if present
        if self.vacancy_concentration > 0:
            vac = f"vac{self.vacancy_concentration:.2f}".replace(".", "")
            name = f"{name}_{vac}"

        return name

    def get_slab_name(self) -> str:
        """
        Get standardized slab name.

        Returns
        -------
        str
            Slab name (e.g., 'LaVO3_001_LaO')
        """
        parts = [self.material_str, self.miller_str]

        if self.termination:
            parts.append(self.termination)

        return "_".join(parts)

    def get_calc_dir(
        self,
        calc_type: str,
        suffix: Optional[str] = None,
        config_id: Optional[int] = None,
    ) -> str:
        """
        Get calculation directory name.

        Parameters
        ----------
        calc_type : str
            Calculation type (relax, neb, freq, md)
        suffix : str, optional
            Additional suffix (e.g., adsorbate name)
        config_id : int, optional
            Configuration number

        Returns
        -------
        str
            Directory name
        """
        parts = [self.get_slab_name(), calc_type]

        if suffix:
            parts.append(suffix)

        if config_id is not None:
            parts.append(f"config{config_id:03d}")

        return "_".join(parts)

    def get_decomposition_dir(
        self,
        step: str,
        config_id: int = 0,
    ) -> str:
        """
        Get decomposition step directory name.

        Parameters
        ----------
        step : str
            Decomposition step (NH3, NH2_H, NH_2H, N_3H)
        config_id : int
            Configuration number

        Returns
        -------
        str
            Directory name
        """
        return self.get_calc_dir("decomp", suffix=step, config_id=config_id)

    def get_neb_dir(
        self,
        initial_step: str,
        final_step: str,
    ) -> str:
        """
        Get NEB calculation directory name.

        Parameters
        ----------
        initial_step : str
            Initial state
        final_step : str
            Final state

        Returns
        -------
        str
            Directory name
        """
        suffix = f"{initial_step}_to_{final_step}"
        return self.get_calc_dir("neb", suffix=suffix)

    def parse_name(self, name: str) -> Dict[str, Any]:
        """
        Parse a directory/file name to extract components.

        Parameters
        ----------
        name : str
            Name to parse

        Returns
        -------
        dict
            Parsed components
        """
        parts = name.split("_")

        result = {
            "material": parts[0] if parts else "",
            "miller": parts[1] if len(parts) > 1 else "",
            "termination": None,
            "calc_type": None,
            "suffix": None,
            "config_id": None,
        }

        # Try to identify other parts
        for i, part in enumerate(parts[2:], 2):
            if part in ["relax", "neb", "freq", "md", "decomp", "screen"]:
                result["calc_type"] = part
            elif part.startswith("config"):
                try:
                    result["config_id"] = int(part[6:])
                except ValueError:
                    pass
            elif i == 2 and part not in ["relax", "neb", "freq", "md"]:
                result["termination"] = part

        return result


class DirectoryStructure:
    """
    Manage calculation directory structure.

    Creates and navigates standardized directory hierarchy.
    """

    def __init__(self, base_dir: str = "./calculations"):
        """
        Initialize DirectoryStructure.

        Parameters
        ----------
        base_dir : str
            Base directory for all calculations
        """
        self.base_dir = Path(base_dir)

    def create_structure(
        self,
        namer: NamingConvention,
        calc_types: Optional[List[str]] = None,
    ) -> Dict[str, Path]:
        """
        Create directory structure for a material.

        Parameters
        ----------
        namer : NamingConvention
            Naming convention object
        calc_types : list, optional
            Calculation types to create directories for

        Returns
        -------
        dict
            Created directory paths
        """
        if calc_types is None:
            calc_types = ["bulk", "slab", "relax", "neb", "freq"]

        paths = {}
        material_dir = self.base_dir / namer.material_str

        for calc_type in calc_types:
            calc_dir = material_dir / calc_type
            calc_dir.mkdir(parents=True, exist_ok=True)
            paths[calc_type] = calc_dir

        return paths

    def get_calc_path(
        self,
        namer: NamingConvention,
        calc_type: str,
        suffix: Optional[str] = None,
        config_id: Optional[int] = None,
    ) -> Path:
        """
        Get path for a specific calculation.

        Parameters
        ----------
        namer : NamingConvention
            Naming convention
        calc_type : str
            Calculation type
        suffix : str, optional
            Additional suffix
        config_id : int, optional
            Configuration ID

        Returns
        -------
        Path
            Calculation directory path
        """
        dir_name = namer.get_calc_dir(calc_type, suffix, config_id)
        return self.base_dir / namer.material_str / calc_type / dir_name

    def find_calculations(
        self,
        material: Optional[str] = None,
        calc_type: Optional[str] = None,
        pattern: Optional[str] = None,
    ) -> List[Path]:
        """
        Find calculation directories matching criteria.

        Parameters
        ----------
        material : str, optional
            Material to filter by
        calc_type : str, optional
            Calculation type to filter by
        pattern : str, optional
            Glob pattern

        Returns
        -------
        list
            Matching directory paths
        """
        if pattern:
            return list(self.base_dir.glob(pattern))

        search_path = self.base_dir
        if material:
            search_path = search_path / material
        if calc_type:
            search_path = search_path / calc_type if material else search_path

        if search_path.exists():
            return [p for p in search_path.iterdir() if p.is_dir()]
        return []

    def get_latest_calculation(
        self,
        namer: NamingConvention,
        calc_type: str,
    ) -> Optional[Path]:
        """
        Get the most recent calculation directory.

        Parameters
        ----------
        namer : NamingConvention
            Naming convention
        calc_type : str
            Calculation type

        Returns
        -------
        Path or None
            Latest calculation path
        """
        base_name = namer.get_calc_dir(calc_type)
        pattern = f"{base_name}*"

        calc_dir = self.base_dir / namer.material_str / calc_type
        if not calc_dir.exists():
            return None

        matches = sorted(calc_dir.glob(pattern), key=lambda p: p.stat().st_mtime)
        return matches[-1] if matches else None


def generate_calc_id(
    material: str,
    miller: tuple,
    calc_type: str,
    adsorbate: Optional[str] = None,
) -> str:
    """
    Generate a unique calculation ID.

    Parameters
    ----------
    material : str
        Material name
    miller : tuple
        Miller indices
    calc_type : str
        Calculation type
    adsorbate : str, optional
        Adsorbate name

    Returns
    -------
    str
        Unique calculation ID
    """
    miller_str = "".join(str(abs(i)) for i in miller)
    parts = [material, miller_str, calc_type]

    if adsorbate:
        parts.append(adsorbate)

    return "_".join(parts)


def parse_calc_id(calc_id: str) -> Dict[str, str]:
    """
    Parse a calculation ID.

    Parameters
    ----------
    calc_id : str
        Calculation ID

    Returns
    -------
    dict
        Parsed components
    """
    parts = calc_id.split("_")

    return {
        "material": parts[0] if parts else "",
        "miller": parts[1] if len(parts) > 1 else "",
        "calc_type": parts[2] if len(parts) > 2 else "",
        "adsorbate": parts[3] if len(parts) > 3 else None,
    }
