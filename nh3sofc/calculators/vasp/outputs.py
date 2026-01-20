"""VASP output file parsing.

Parses OUTCAR, vasprun.xml, CONTCAR, and other VASP output files.
"""

import re
from pathlib import Path
from typing import Optional, List, Union, Dict, Any, Tuple
import numpy as np
from ase import Atoms
from ase.io import read as ase_read


class VASPOutputParser:
    """
    Parser for VASP output files.

    Extracts energies, forces, stresses, and other properties from
    VASP calculations.

    Examples
    --------
    >>> parser = VASPOutputParser("./calculation")
    >>> energy = parser.get_energy()
    >>> forces = parser.get_forces()
    >>> atoms = parser.get_final_structure()
    """

    def __init__(self, work_dir: Union[str, Path]):
        """
        Initialize VASPOutputParser.

        Parameters
        ----------
        work_dir : str or Path
            Directory containing VASP output files
        """
        self.work_dir = Path(work_dir)

        # Check for output files
        self.outcar_path = self.work_dir / "OUTCAR"
        self.vasprun_path = self.work_dir / "vasprun.xml"
        self.contcar_path = self.work_dir / "CONTCAR"
        self.oszicar_path = self.work_dir / "OSZICAR"

    def is_converged(self) -> bool:
        """
        Check if the calculation converged.

        Returns
        -------
        bool
            True if calculation converged
        """
        if not self.outcar_path.exists():
            return False

        with open(self.outcar_path, "r") as f:
            content = f.read()

        # Check for ionic convergence
        if "reached required accuracy" in content:
            return True

        # Check for electronic convergence on last ionic step
        if "EDIFF is reached" in content:
            return True

        return False

    def get_energy(self, per_atom: bool = False) -> Optional[float]:
        """
        Get the final total energy.

        Parameters
        ----------
        per_atom : bool
            Return energy per atom

        Returns
        -------
        float or None
            Total energy in eV
        """
        if not self.outcar_path.exists():
            return None

        energy = None

        with open(self.outcar_path, "r") as f:
            for line in f:
                # Look for "energy without entropy" (most accurate)
                if "energy  without entropy" in line:
                    parts = line.split()
                    energy = float(parts[-1])
                # Fallback to TOTEN
                elif "TOTEN" in line and energy is None:
                    match = re.search(r"TOTEN\s*=\s*([-\d.]+)", line)
                    if match:
                        energy = float(match.group(1))

        if energy is not None and per_atom:
            n_atoms = self.get_n_atoms()
            if n_atoms:
                energy /= n_atoms

        return energy

    def get_energy_series(self) -> List[float]:
        """
        Get energies from all ionic steps.

        Returns
        -------
        list
            List of energies for each ionic step
        """
        if not self.oszicar_path.exists():
            return []

        energies = []
        with open(self.oszicar_path, "r") as f:
            for line in f:
                if "E0=" in line:
                    match = re.search(r"E0=\s*([-\d.E+]+)", line)
                    if match:
                        energies.append(float(match.group(1)))

        return energies

    def get_forces(self) -> Optional[np.ndarray]:
        """
        Get forces on atoms.

        Returns
        -------
        ndarray or None
            Forces array of shape (n_atoms, 3) in eV/A
        """
        if not self.outcar_path.exists():
            return None

        forces = []
        in_forces_block = False

        with open(self.outcar_path, "r") as f:
            for line in f:
                if "TOTAL-FORCE" in line:
                    in_forces_block = True
                    forces = []
                    next(f)  # Skip header line
                    continue

                if in_forces_block:
                    if "---" in line:
                        in_forces_block = False
                        continue
                    parts = line.split()
                    if len(parts) >= 6:
                        forces.append([float(parts[3]), float(parts[4]), float(parts[5])])

        if forces:
            return np.array(forces)
        return None

    def get_stress(self) -> Optional[np.ndarray]:
        """
        Get stress tensor.

        Returns
        -------
        ndarray or None
            Stress tensor (Voigt notation) in kBar
        """
        if not self.outcar_path.exists():
            return None

        stress = None

        with open(self.outcar_path, "r") as f:
            for line in f:
                if "in kB" in line:
                    parts = line.split()
                    if len(parts) >= 8:
                        stress = np.array([float(x) for x in parts[2:8]])

        return stress

    def get_max_force(self) -> Optional[float]:
        """
        Get maximum force magnitude.

        Returns
        -------
        float or None
            Maximum force in eV/A
        """
        forces = self.get_forces()
        if forces is None:
            return None

        force_magnitudes = np.linalg.norm(forces, axis=1)
        return float(np.max(force_magnitudes))

    def get_n_atoms(self) -> Optional[int]:
        """
        Get number of atoms.

        Returns
        -------
        int or None
            Number of atoms
        """
        if self.contcar_path.exists():
            try:
                atoms = ase_read(str(self.contcar_path), format="vasp")
                return len(atoms)
            except Exception:
                pass

        poscar_path = self.work_dir / "POSCAR"
        if poscar_path.exists():
            try:
                atoms = ase_read(str(poscar_path), format="vasp")
                return len(atoms)
            except Exception:
                pass

        return None

    def get_final_structure(self) -> Optional[Atoms]:
        """
        Get the final optimized structure.

        Returns
        -------
        Atoms or None
            Final structure
        """
        if self.contcar_path.exists():
            try:
                return ase_read(str(self.contcar_path), format="vasp")
            except Exception:
                pass

        # Try vasprun.xml
        if self.vasprun_path.exists():
            try:
                return ase_read(str(self.vasprun_path), format="vasp-xml", index=-1)
            except Exception:
                pass

        return None

    def get_trajectory(self) -> Optional[List[Atoms]]:
        """
        Get trajectory from all ionic steps.

        Returns
        -------
        list or None
            List of Atoms objects for each step
        """
        if self.vasprun_path.exists():
            try:
                return ase_read(str(self.vasprun_path), format="vasp-xml", index=":")
            except Exception:
                pass

        return None

    def get_fermi_energy(self) -> Optional[float]:
        """
        Get Fermi energy.

        Returns
        -------
        float or None
            Fermi energy in eV
        """
        if not self.outcar_path.exists():
            return None

        with open(self.outcar_path, "r") as f:
            for line in f:
                if "E-fermi" in line:
                    parts = line.split()
                    return float(parts[2])

        return None

    def get_band_gap(self) -> Optional[Dict[str, float]]:
        """
        Get band gap information.

        Returns
        -------
        dict or None
            Band gap information
        """
        if not self.outcar_path.exists():
            return None

        vbm = None
        cbm = None

        # This is a simplified approach; proper band gap calculation
        # requires parsing eigenvalues from vasprun.xml
        with open(self.outcar_path, "r") as f:
            content = f.read()

        # Look for direct/indirect gap info if available
        gap_match = re.search(r"band gap\s*=\s*([\d.]+)", content, re.IGNORECASE)
        if gap_match:
            return {"gap": float(gap_match.group(1))}

        return None

    def get_magnetic_moment(self) -> Optional[float]:
        """
        Get total magnetic moment.

        Returns
        -------
        float or None
            Total magnetic moment in Bohr magnetons
        """
        if not self.outcar_path.exists():
            return None

        mag = None
        with open(self.outcar_path, "r") as f:
            for line in f:
                if "number of electron" in line and "magnetization" in line:
                    parts = line.split()
                    mag = float(parts[-1])

        return mag

    def get_elapsed_time(self) -> Optional[float]:
        """
        Get total elapsed time.

        Returns
        -------
        float or None
            Elapsed time in seconds
        """
        if not self.outcar_path.exists():
            return None

        with open(self.outcar_path, "r") as f:
            content = f.read()

        match = re.search(r"Elapsed time \(sec\):\s*([\d.]+)", content)
        if match:
            return float(match.group(1))

        return None

    def get_n_ionic_steps(self) -> int:
        """
        Get number of ionic steps performed.

        Returns
        -------
        int
            Number of ionic steps
        """
        return len(self.get_energy_series())

    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of the calculation results.

        Returns
        -------
        dict
            Summary of calculation results
        """
        return {
            "converged": self.is_converged(),
            "energy": self.get_energy(),
            "energy_per_atom": self.get_energy(per_atom=True),
            "max_force": self.get_max_force(),
            "n_atoms": self.get_n_atoms(),
            "n_ionic_steps": self.get_n_ionic_steps(),
            "magnetic_moment": self.get_magnetic_moment(),
            "fermi_energy": self.get_fermi_energy(),
            "elapsed_time": self.get_elapsed_time(),
        }

    def __repr__(self) -> str:
        status = "converged" if self.is_converged() else "not converged"
        energy = self.get_energy()
        energy_str = f"{energy:.4f} eV" if energy else "N/A"
        return f"VASPOutputParser({self.work_dir}, {status}, E={energy_str})"


class NEBOutputParser:
    """
    Parser for NEB calculation outputs.
    """

    def __init__(self, work_dir: Union[str, Path]):
        """
        Initialize NEBOutputParser.

        Parameters
        ----------
        work_dir : str or Path
            NEB calculation directory
        """
        self.work_dir = Path(work_dir)

    def get_image_dirs(self) -> List[Path]:
        """Get directories for each NEB image."""
        dirs = []
        for d in sorted(self.work_dir.iterdir()):
            if d.is_dir() and d.name.isdigit():
                dirs.append(d)
        return dirs

    def get_energies(self) -> List[float]:
        """
        Get energies for all images.

        Returns
        -------
        list
            Energies for each image
        """
        energies = []
        for image_dir in self.get_image_dirs():
            parser = VASPOutputParser(image_dir)
            energy = parser.get_energy()
            if energy is not None:
                energies.append(energy)
        return energies

    def get_barrier(self) -> Tuple[float, float]:
        """
        Get forward and reverse barriers.

        Returns
        -------
        tuple
            (forward_barrier, reverse_barrier)
        """
        energies = self.get_energies()
        if len(energies) < 3:
            return (0.0, 0.0)

        e_initial = energies[0]
        e_final = energies[-1]
        e_max = max(energies)

        forward = e_max - e_initial
        reverse = e_max - e_final

        return (forward, reverse)

    def get_reaction_energy(self) -> float:
        """
        Get reaction energy (E_final - E_initial).

        Returns
        -------
        float
            Reaction energy
        """
        energies = self.get_energies()
        if len(energies) < 2:
            return 0.0
        return energies[-1] - energies[0]


def parse_vasp_calculation(work_dir: Union[str, Path]) -> Dict[str, Any]:
    """
    Convenience function to parse a VASP calculation.

    Parameters
    ----------
    work_dir : str or Path
        Calculation directory

    Returns
    -------
    dict
        Calculation results
    """
    parser = VASPOutputParser(work_dir)
    results = parser.get_summary()
    results["final_structure"] = parser.get_final_structure()
    return results
