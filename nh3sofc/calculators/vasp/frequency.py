"""VASP frequency/vibrational analysis setup and parsing.

Supports:
- Frequency calculations (IBRION=5,6,7)
- Partial Hessian for adsorbates only
- Thermochemistry from frequencies
"""

from pathlib import Path
from typing import Optional, List, Union, Dict, Any
import numpy as np
from ase import Atoms

from .inputs import VASPInputGenerator
from .outputs import VASPOutputParser
from ...core.constants import KB_EV, H_EV_S


class FrequencyCalculation:
    """
    Setup and analysis of VASP frequency calculations.

    Examples
    --------
    >>> freq = FrequencyCalculation(atoms, adsorbate_indices=[36, 37, 38, 39])
    >>> freq.generate_inputs(work_dir="./freq_calc")
    >>>
    >>> # After calculation completes:
    >>> freq.parse_outputs(work_dir="./freq_calc")
    >>> zpe = freq.get_zpe()
    """

    def __init__(
        self,
        atoms: Atoms,
        adsorbate_indices: Optional[List[int]] = None,
        method: str = "finite_difference",
    ):
        """
        Initialize FrequencyCalculation.

        Parameters
        ----------
        atoms : Atoms
            Optimized structure
        adsorbate_indices : list, optional
            Indices of adsorbate atoms. If provided, only these atoms
            are displaced (partial Hessian).
        method : str
            Method: "finite_difference" (IBRION=5/6) or "dfpt" (IBRION=7/8)
        """
        self.atoms = atoms.copy()
        self.adsorbate_indices = adsorbate_indices
        self.method = method

        # Parsed data
        self.frequencies = None  # in cm^-1
        self.eigenvectors = None
        self.intensities = None

    def generate_inputs(
        self,
        work_dir: Union[str, Path],
        nfree: int = 2,
        potim: float = 0.015,
        encut: Optional[float] = None,
        kspacing: Optional[float] = None,
        hubbard_u: Optional[Dict[str, float]] = None,
        vdw: Optional[str] = None,
    ) -> Dict[str, Path]:
        """
        Generate VASP input files for frequency calculation.

        Parameters
        ----------
        work_dir : str or Path
            Output directory
        nfree : int
            Number of displacements (1 or 2)
        potim : float
            Displacement step size in Angstrom
        encut : float, optional
            Cutoff energy
        kspacing : float, optional
            K-point spacing
        hubbard_u : dict, optional
            Hubbard U values
        vdw : str, optional
            VdW correction

        Returns
        -------
        dict
            Generated file paths
        """
        work_dir = Path(work_dir)

        # Use partial Hessian if adsorbate indices specified
        atoms = self.atoms.copy()

        if self.adsorbate_indices is not None:
            # Set selective dynamics - only displace adsorbate atoms
            from ase.constraints import FixAtoms

            # Fix all atoms except adsorbates
            fix_indices = [i for i in range(len(atoms))
                          if i not in self.adsorbate_indices]
            atoms.set_constraint(FixAtoms(indices=fix_indices))

        # Create generator with frequency preset
        if self.method == "dfpt":
            ibrion = 7  # DFPT for insulators
        else:
            ibrion = 5  # Finite differences

        gen = VASPInputGenerator(
            atoms,
            calc_type="frequency",
            work_dir=work_dir,
            IBRION=ibrion,
            NFREE=nfree,
            POTIM=potim,
            EDIFF=1e-7,  # Tighter convergence for frequencies
            NSW=1,
        )

        return gen.generate_all(
            encut=encut,
            kspacing=kspacing,
            hubbard_u=hubbard_u,
            vdw=vdw,
            is_surface=True,
        )

    def parse_outputs(self, work_dir: Union[str, Path]) -> None:
        """
        Parse frequency calculation outputs.

        Parameters
        ----------
        work_dir : str or Path
            Calculation directory
        """
        work_dir = Path(work_dir)
        outcar_path = work_dir / "OUTCAR"

        if not outcar_path.exists():
            raise FileNotFoundError(f"OUTCAR not found in {work_dir}")

        frequencies = []
        imaginary = []

        with open(outcar_path, "r") as f:
            in_freq_block = False

            for line in f:
                # Look for frequency block
                if "Eigenvectors and eigenvalues of the dynamical matrix" in line:
                    in_freq_block = True
                    frequencies = []
                    imaginary = []
                    continue

                if in_freq_block:
                    # Parse frequency lines
                    if "f  =" in line or "f/i=" in line:
                        parts = line.split()
                        # Find the cm-1 value
                        for i, p in enumerate(parts):
                            if p == "cm-1":
                                freq = float(parts[i - 1])
                                frequencies.append(freq)
                                imaginary.append("f/i" in line)
                                break

                    # End of frequency block
                    if "---" in line and frequencies:
                        break

        self.frequencies = np.array(frequencies)
        self.imaginary = np.array(imaginary)

    def get_frequencies(self, include_imaginary: bool = False) -> np.ndarray:
        """
        Get vibrational frequencies.

        Parameters
        ----------
        include_imaginary : bool
            Include imaginary frequencies (negative values)

        Returns
        -------
        ndarray
            Frequencies in cm^-1
        """
        if self.frequencies is None:
            raise ValueError("No frequencies parsed. Run parse_outputs first.")

        if include_imaginary:
            freqs = self.frequencies.copy()
            freqs[self.imaginary] *= -1  # Make imaginary negative
            return freqs
        else:
            return self.frequencies[~self.imaginary]

    def get_frequencies_ev(self, include_imaginary: bool = False) -> np.ndarray:
        """
        Get frequencies in eV.

        Parameters
        ----------
        include_imaginary : bool
            Include imaginary frequencies

        Returns
        -------
        ndarray
            Frequencies in eV
        """
        cm1_to_ev = 1.23981e-4  # 1 cm^-1 in eV
        return self.get_frequencies(include_imaginary) * cm1_to_ev

    def get_zpe(self) -> float:
        """
        Calculate zero-point energy.

        ZPE = (1/2) * sum(h * nu_i)

        Returns
        -------
        float
            Zero-point energy in eV
        """
        freqs_ev = self.get_frequencies_ev(include_imaginary=False)
        return 0.5 * np.sum(freqs_ev)

    def get_vibrational_entropy(self, temperature: float) -> float:
        """
        Calculate vibrational entropy.

        S_vib = sum_i [ (h*nu_i/kT) / (exp(h*nu_i/kT) - 1) - ln(1 - exp(-h*nu_i/kT)) ]

        Parameters
        ----------
        temperature : float
            Temperature in K

        Returns
        -------
        float
            Vibrational entropy in eV/K
        """
        freqs_ev = self.get_frequencies_ev(include_imaginary=False)
        kT = KB_EV * temperature

        S = 0.0
        for nu in freqs_ev:
            if nu > 0:
                x = nu / kT
                if x < 100:  # Avoid overflow
                    S += x / (np.exp(x) - 1) - np.log(1 - np.exp(-x))

        return KB_EV * S

    def get_vibrational_energy(self, temperature: float) -> float:
        """
        Calculate vibrational internal energy.

        E_vib = sum_i [ h*nu_i * (1/2 + 1/(exp(h*nu_i/kT) - 1)) ]

        Parameters
        ----------
        temperature : float
            Temperature in K

        Returns
        -------
        float
            Vibrational energy in eV
        """
        freqs_ev = self.get_frequencies_ev(include_imaginary=False)
        kT = KB_EV * temperature

        E = 0.0
        for nu in freqs_ev:
            if nu > 0:
                x = nu / kT
                if x < 100:
                    E += nu * (0.5 + 1 / (np.exp(x) - 1))
                else:
                    E += nu * 0.5  # Ground state only

        return E

    def get_helmholtz_energy(self, temperature: float) -> float:
        """
        Calculate vibrational Helmholtz free energy.

        A_vib = E_vib - T * S_vib = ZPE + kT * sum_i ln(1 - exp(-h*nu_i/kT))

        Parameters
        ----------
        temperature : float
            Temperature in K

        Returns
        -------
        float
            Helmholtz free energy in eV
        """
        freqs_ev = self.get_frequencies_ev(include_imaginary=False)
        kT = KB_EV * temperature

        A = 0.0
        for nu in freqs_ev:
            if nu > 0:
                x = nu / kT
                A += 0.5 * nu  # ZPE contribution
                if x < 100:
                    A += kT * np.log(1 - np.exp(-x))

        return A

    def get_thermal_correction(self, temperature: float) -> Dict[str, float]:
        """
        Get all thermal corrections.

        Parameters
        ----------
        temperature : float
            Temperature in K

        Returns
        -------
        dict
            Dictionary with ZPE, E_vib, S_vib, A_vib, G_correction
        """
        zpe = self.get_zpe()
        e_vib = self.get_vibrational_energy(temperature)
        s_vib = self.get_vibrational_entropy(temperature)
        a_vib = self.get_helmholtz_energy(temperature)

        # For surface species, A ≈ G (no PV term)
        return {
            "ZPE": zpe,
            "E_vib": e_vib,
            "S_vib": s_vib,
            "TS_vib": temperature * s_vib,
            "A_vib": a_vib,
            "temperature": temperature,
        }

    def get_summary(self) -> Dict[str, Any]:
        """
        Get summary of frequency calculation.

        Returns
        -------
        dict
            Summary including frequencies, ZPE, etc.
        """
        if self.frequencies is None:
            return {"parsed": False}

        n_real = np.sum(~self.imaginary)
        n_imag = np.sum(self.imaginary)

        return {
            "parsed": True,
            "n_frequencies": len(self.frequencies),
            "n_real": n_real,
            "n_imaginary": n_imag,
            "frequencies_cm1": self.get_frequencies().tolist(),
            "imaginary_cm1": self.frequencies[self.imaginary].tolist() if n_imag > 0 else [],
            "ZPE_eV": self.get_zpe(),
            "ZPE_kJ_mol": self.get_zpe() * 96.485,
        }


def calculate_zpe_from_outcar(outcar_path: Union[str, Path]) -> float:
    """
    Quick function to calculate ZPE from OUTCAR.

    Parameters
    ----------
    outcar_path : str or Path
        Path to OUTCAR file

    Returns
    -------
    float
        Zero-point energy in eV
    """
    freq_calc = FrequencyCalculation(Atoms())  # Dummy atoms
    freq_calc.parse_outputs(Path(outcar_path).parent)
    return freq_calc.get_zpe()


def get_thermal_corrections(
    outcar_path: Union[str, Path],
    temperature: float,
) -> Dict[str, float]:
    """
    Get thermal corrections from OUTCAR.

    Parameters
    ----------
    outcar_path : str or Path
        Path to OUTCAR
    temperature : float
        Temperature in K

    Returns
    -------
    dict
        Thermal corrections
    """
    freq_calc = FrequencyCalculation(Atoms())
    freq_calc.parse_outputs(Path(outcar_path).parent)
    return freq_calc.get_thermal_correction(temperature)
