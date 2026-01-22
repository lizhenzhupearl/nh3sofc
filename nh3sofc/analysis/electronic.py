"""Electronic structure analysis for catalyst descriptors.

Provides tools for extracting electronic structure descriptors from
DFT calculations, including d-band center analysis from VASP DOSCAR files.
"""

from pathlib import Path
from typing import Optional, List, Union, Dict, Any, Tuple
import numpy as np


class DBandAnalyzer:
    """
    Analyze d-band electronic structure from VASP DOSCAR files.

    The d-band center is a key descriptor for predicting adsorption
    energies on transition metal surfaces (Hammer-Norskov d-band model).

    The d-band center is defined as:
        ε_d = ∫ ε * ρ_d(ε) dε / ∫ ρ_d(ε) dε

    where ρ_d(ε) is the d-band density of states.

    Examples
    --------
    >>> analyzer = DBandAnalyzer.from_doscar("DOSCAR")
    >>> d_center = analyzer.get_d_band_center(atom_index=0)
    >>> print(f"d-band center: {d_center:.3f} eV")

    >>> # Get d-band properties for surface atoms
    >>> props = analyzer.get_d_band_properties([0, 1, 2])
    >>> for atom, data in props.items():
    ...     print(f"Atom {atom}: center={data['center']:.3f}, width={data['width']:.3f}")
    """

    def __init__(
        self,
        energy: np.ndarray,
        dos_total: np.ndarray,
        dos_projected: Optional[Dict[int, Dict[str, np.ndarray]]] = None,
        e_fermi: float = 0.0,
        is_spin_polarized: bool = False,
    ):
        """
        Initialize DBandAnalyzer.

        Parameters
        ----------
        energy : np.ndarray
            Energy grid (eV)
        dos_total : np.ndarray
            Total density of states
        dos_projected : dict, optional
            Projected DOS per atom: {atom_index: {'s': dos_s, 'p': dos_p, 'd': dos_d, ...}}
        e_fermi : float
            Fermi energy (eV)
        is_spin_polarized : bool
            Whether calculation is spin-polarized
        """
        self.energy = np.array(energy)
        self.dos_total = np.array(dos_total)
        self.dos_projected = dos_projected or {}
        self.e_fermi = e_fermi
        self.is_spin_polarized = is_spin_polarized

        # Shift energy to Fermi level
        self.energy_shifted = self.energy - self.e_fermi

    @classmethod
    def from_doscar(
        cls,
        filepath: Union[str, Path],
        atoms_to_parse: Optional[List[int]] = None,
    ) -> "DBandAnalyzer":
        """
        Create DBandAnalyzer from VASP DOSCAR file.

        Parameters
        ----------
        filepath : str or Path
            Path to DOSCAR file
        atoms_to_parse : list, optional
            Specific atom indices to parse (0-indexed). If None, parse all.

        Returns
        -------
        DBandAnalyzer
            Analyzer instance with parsed DOS data
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"DOSCAR not found: {filepath}")

        with open(filepath, "r") as f:
            lines = f.readlines()

        # Parse header
        # Line 1: natom, natom, ?, ?
        # Line 6: E_max, E_min, NEDOS, E_fermi, ?
        natom = int(lines[0].split()[0])
        header = lines[5].split()
        e_max = float(header[0])
        e_min = float(header[1])
        nedos = int(header[2])
        e_fermi = float(header[3])

        # Parse total DOS (starts at line 6)
        dos_start = 6
        energy = []
        dos_total = []
        dos_integrated = []

        for i in range(nedos):
            parts = lines[dos_start + i].split()
            energy.append(float(parts[0]))
            dos_total.append(float(parts[1]))
            if len(parts) > 2:
                # Could be spin-polarized or integrated DOS
                if len(parts) >= 3:
                    dos_integrated.append(float(parts[2]))

        energy = np.array(energy)
        dos_total = np.array(dos_total)

        # Detect spin polarization from total DOS columns
        is_spin_polarized = len(lines[dos_start].split()) > 3

        # Parse projected DOS for each atom
        dos_projected = {}
        atom_dos_start = dos_start + nedos

        # Check if projected DOS is present
        if len(lines) > atom_dos_start:
            if atoms_to_parse is None:
                atoms_to_parse = list(range(natom))

            for atom_idx in range(natom):
                # Each atom block starts with a header line, then NEDOS lines
                block_start = atom_dos_start + atom_idx * (nedos + 1) + 1

                if block_start + nedos > len(lines):
                    break

                if atom_idx not in atoms_to_parse:
                    continue

                # Parse orbital-projected DOS
                # VASP format depends on LORBIT setting
                # LORBIT=10: s, p, d (for each spin if spin-polarized)
                # LORBIT=11: s, py, pz, px, dxy, dyz, dz2, dxz, dx2-y2

                dos_s = []
                dos_p = []
                dos_d = []

                for i in range(nedos):
                    line_idx = block_start + i
                    if line_idx >= len(lines):
                        break
                    parts = lines[line_idx].split()

                    if len(parts) >= 4:
                        # At minimum: energy, s, p, d
                        dos_s.append(float(parts[1]))
                        dos_p.append(float(parts[2]))
                        dos_d.append(float(parts[3]))
                    elif len(parts) >= 10:
                        # LORBIT=11 format: energy, s, py, pz, px, dxy, dyz, dz2, dxz, dx2-y2
                        dos_s.append(float(parts[1]))
                        dos_p.append(sum(float(parts[j]) for j in [2, 3, 4]))
                        dos_d.append(sum(float(parts[j]) for j in [5, 6, 7, 8, 9]))

                if dos_s:
                    dos_projected[atom_idx] = {
                        "s": np.array(dos_s),
                        "p": np.array(dos_p),
                        "d": np.array(dos_d),
                    }

        return cls(
            energy=energy,
            dos_total=dos_total,
            dos_projected=dos_projected,
            e_fermi=e_fermi,
            is_spin_polarized=is_spin_polarized,
        )

    def get_d_band_center(
        self,
        atom_index: int,
        energy_range: Optional[Tuple[float, float]] = None,
        use_fermi_shifted: bool = True,
    ) -> float:
        """
        Calculate d-band center for a specific atom.

        ε_d = ∫ ε * ρ_d(ε) dε / ∫ ρ_d(ε) dε

        Parameters
        ----------
        atom_index : int
            Atom index (0-indexed)
        energy_range : tuple, optional
            (E_min, E_max) to integrate over. Default: (-10, 2) eV relative to Fermi
        use_fermi_shifted : bool
            If True, return center relative to Fermi level

        Returns
        -------
        float
            D-band center in eV
        """
        if atom_index not in self.dos_projected:
            raise ValueError(
                f"No projected DOS for atom {atom_index}. "
                f"Available atoms: {list(self.dos_projected.keys())}"
            )

        dos_d = self.dos_projected[atom_index]["d"]
        energy = self.energy_shifted if use_fermi_shifted else self.energy

        # Default energy range: typical d-band region
        if energy_range is None:
            energy_range = (-10.0, 2.0)

        # Apply energy mask
        mask = (energy >= energy_range[0]) & (energy <= energy_range[1])
        e_masked = energy[mask]
        d_masked = dos_d[mask]

        # Numerical integration using trapezoidal rule
        numerator = np.trapz(e_masked * d_masked, e_masked)
        denominator = np.trapz(d_masked, e_masked)

        if abs(denominator) < 1e-10:
            raise ValueError(f"No d-band DOS found for atom {atom_index}")

        return numerator / denominator

    def get_d_band_width(
        self,
        atom_index: int,
        energy_range: Optional[Tuple[float, float]] = None,
    ) -> float:
        """
        Calculate d-band width (second moment) for a specific atom.

        W_d = sqrt(∫ (ε - ε_d)² * ρ_d(ε) dε / ∫ ρ_d(ε) dε)

        Parameters
        ----------
        atom_index : int
            Atom index (0-indexed)
        energy_range : tuple, optional
            Energy range for integration

        Returns
        -------
        float
            D-band width in eV
        """
        if atom_index not in self.dos_projected:
            raise ValueError(f"No projected DOS for atom {atom_index}")

        center = self.get_d_band_center(atom_index, energy_range)
        dos_d = self.dos_projected[atom_index]["d"]
        energy = self.energy_shifted

        if energy_range is None:
            energy_range = (-10.0, 2.0)

        mask = (energy >= energy_range[0]) & (energy <= energy_range[1])
        e_masked = energy[mask]
        d_masked = dos_d[mask]

        # Second moment
        numerator = np.trapz((e_masked - center) ** 2 * d_masked, e_masked)
        denominator = np.trapz(d_masked, e_masked)

        if abs(denominator) < 1e-10:
            return 0.0

        return np.sqrt(numerator / denominator)

    def get_d_band_filling(
        self,
        atom_index: int,
        energy_range: Optional[Tuple[float, float]] = None,
    ) -> float:
        """
        Calculate d-band filling (fraction of occupied d-states).

        f_d = ∫_{-∞}^{E_F} ρ_d(ε) dε / ∫ ρ_d(ε) dε

        Parameters
        ----------
        atom_index : int
            Atom index (0-indexed)
        energy_range : tuple, optional
            Total energy range for d-band

        Returns
        -------
        float
            D-band filling (0 to 1)
        """
        if atom_index not in self.dos_projected:
            raise ValueError(f"No projected DOS for atom {atom_index}")

        dos_d = self.dos_projected[atom_index]["d"]
        energy = self.energy_shifted

        if energy_range is None:
            energy_range = (-10.0, 5.0)

        # Full d-band
        mask_full = (energy >= energy_range[0]) & (energy <= energy_range[1])
        total = np.trapz(dos_d[mask_full], energy[mask_full])

        # Occupied (below Fermi level)
        mask_occ = (energy >= energy_range[0]) & (energy <= 0.0)
        occupied = np.trapz(dos_d[mask_occ], energy[mask_occ])

        if abs(total) < 1e-10:
            return 0.0

        return occupied / total

    def get_d_band_properties(
        self,
        atom_indices: Optional[List[int]] = None,
        energy_range: Optional[Tuple[float, float]] = None,
    ) -> Dict[int, Dict[str, float]]:
        """
        Get comprehensive d-band properties for multiple atoms.

        Parameters
        ----------
        atom_indices : list, optional
            Atom indices to analyze. Default: all available
        energy_range : tuple, optional
            Energy range for integration

        Returns
        -------
        dict
            {atom_index: {'center': float, 'width': float, 'filling': float}}
        """
        if atom_indices is None:
            atom_indices = list(self.dos_projected.keys())

        results = {}
        for idx in atom_indices:
            if idx in self.dos_projected:
                results[idx] = {
                    "center": self.get_d_band_center(idx, energy_range),
                    "width": self.get_d_band_width(idx, energy_range),
                    "filling": self.get_d_band_filling(idx, energy_range),
                }

        return results

    def get_surface_average_d_band_center(
        self,
        surface_atom_indices: List[int],
        weights: Optional[List[float]] = None,
        energy_range: Optional[Tuple[float, float]] = None,
    ) -> float:
        """
        Calculate weighted average d-band center for surface atoms.

        Parameters
        ----------
        surface_atom_indices : list
            Indices of surface atoms
        weights : list, optional
            Weights for averaging (e.g., coordination numbers)
        energy_range : tuple, optional
            Energy range for integration

        Returns
        -------
        float
            Average d-band center in eV
        """
        centers = []
        for idx in surface_atom_indices:
            if idx in self.dos_projected:
                centers.append(self.get_d_band_center(idx, energy_range))

        if not centers:
            raise ValueError("No d-band data for specified surface atoms")

        if weights is None:
            return np.mean(centers)
        else:
            weights = np.array(weights[: len(centers)])
            return np.average(centers, weights=weights)


class DOSAnalyzer:
    """
    General density of states analysis utilities.

    Provides tools for analyzing total and projected DOS beyond
    d-band specific analysis.

    Examples
    --------
    >>> analyzer = DOSAnalyzer.from_doscar("DOSCAR")
    >>> gap = analyzer.get_band_gap()
    >>> print(f"Band gap: {gap:.3f} eV")
    """

    def __init__(
        self,
        energy: np.ndarray,
        dos: np.ndarray,
        e_fermi: float = 0.0,
    ):
        """
        Initialize DOSAnalyzer.

        Parameters
        ----------
        energy : np.ndarray
            Energy grid (eV)
        dos : np.ndarray
            Density of states
        e_fermi : float
            Fermi energy
        """
        self.energy = np.array(energy)
        self.dos = np.array(dos)
        self.e_fermi = e_fermi
        self.energy_shifted = self.energy - self.e_fermi

    @classmethod
    def from_doscar(cls, filepath: Union[str, Path]) -> "DOSAnalyzer":
        """Create DOSAnalyzer from VASP DOSCAR file."""
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"DOSCAR not found: {filepath}")

        with open(filepath, "r") as f:
            lines = f.readlines()

        # Parse header
        header = lines[5].split()
        nedos = int(header[2])
        e_fermi = float(header[3])

        # Parse total DOS
        energy = []
        dos = []

        for i in range(nedos):
            parts = lines[6 + i].split()
            energy.append(float(parts[0]))
            dos.append(float(parts[1]))

        return cls(
            energy=np.array(energy),
            dos=np.array(dos),
            e_fermi=e_fermi,
        )

    def get_band_gap(self, threshold: float = 0.01) -> float:
        """
        Calculate band gap from DOS.

        Parameters
        ----------
        threshold : float
            DOS threshold for identifying gap region

        Returns
        -------
        float
            Band gap in eV (0 if metallic)
        """
        # Find where DOS is below threshold
        low_dos = self.dos < threshold

        # Find gap around Fermi level
        fermi_idx = np.argmin(np.abs(self.energy_shifted))

        # Search for gap below Fermi level (VBM)
        vbm_idx = fermi_idx
        for i in range(fermi_idx, -1, -1):
            if not low_dos[i]:
                vbm_idx = i
                break

        # Search for gap above Fermi level (CBM)
        cbm_idx = fermi_idx
        for i in range(fermi_idx, len(low_dos)):
            if not low_dos[i]:
                cbm_idx = i
                break

        if vbm_idx == cbm_idx:
            return 0.0  # Metallic

        gap = self.energy_shifted[cbm_idx] - self.energy_shifted[vbm_idx]
        return max(0.0, gap)

    def get_dos_at_fermi(self) -> float:
        """Get DOS at Fermi level (indicator of metallicity)."""
        fermi_idx = np.argmin(np.abs(self.energy_shifted))
        return self.dos[fermi_idx]


def calculate_d_band_center(
    doscar_path: Union[str, Path],
    atom_indices: List[int],
    energy_range: Optional[Tuple[float, float]] = None,
) -> Dict[int, float]:
    """
    Convenience function to calculate d-band centers from DOSCAR.

    Parameters
    ----------
    doscar_path : str or Path
        Path to VASP DOSCAR file
    atom_indices : list
        Atom indices to analyze (0-indexed)
    energy_range : tuple, optional
        Energy range for integration (default: -10 to 2 eV)

    Returns
    -------
    dict
        {atom_index: d_band_center}

    Examples
    --------
    >>> centers = calculate_d_band_center("DOSCAR", [0, 1, 2])
    >>> print(f"Average d-band center: {np.mean(list(centers.values())):.3f} eV")
    """
    analyzer = DBandAnalyzer.from_doscar(doscar_path, atoms_to_parse=atom_indices)

    results = {}
    for idx in atom_indices:
        try:
            results[idx] = analyzer.get_d_band_center(idx, energy_range)
        except ValueError:
            pass

    return results


def get_surface_d_band_center(
    doscar_path: Union[str, Path],
    surface_atom_indices: List[int],
    energy_range: Optional[Tuple[float, float]] = None,
) -> Tuple[float, Dict[int, float]]:
    """
    Calculate average d-band center for surface atoms.

    Parameters
    ----------
    doscar_path : str or Path
        Path to VASP DOSCAR file
    surface_atom_indices : list
        Indices of surface atoms
    energy_range : tuple, optional
        Energy range for integration

    Returns
    -------
    tuple
        (average_center, {atom_index: center})

    Examples
    --------
    >>> avg, per_atom = get_surface_d_band_center("DOSCAR", [0, 1, 2, 3])
    >>> print(f"Surface average d-band center: {avg:.3f} eV")
    """
    analyzer = DBandAnalyzer.from_doscar(doscar_path, atoms_to_parse=surface_atom_indices)

    centers = {}
    for idx in surface_atom_indices:
        try:
            centers[idx] = analyzer.get_d_band_center(idx, energy_range)
        except ValueError:
            pass

    if not centers:
        raise ValueError("No d-band data found for surface atoms")

    average = np.mean(list(centers.values()))
    return average, centers
