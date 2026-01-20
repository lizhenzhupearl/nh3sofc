"""Frequency and thermochemistry workflow.

Provides automated vibrational analysis and thermodynamic property
calculations including ZPE, entropy, and Gibbs free energy.
"""

from pathlib import Path
from typing import Optional, List, Union, Dict, Any
import numpy as np
from ase import Atoms
from ase.io import write as ase_write, read as ase_read

from ..core.base import BaseWorkflow
from ..core.constants import KB_EV, H_EV_S, C_CMS, EV_TO_J
from ..calculators.vasp.inputs import VASPInputGenerator
from ..calculators.vasp.outputs import VASPOutputParser
from ..calculators.vasp.frequency import (
    FrequencyCalculation,
    calculate_zpe_from_outcar,
    get_thermal_corrections,
)
from ..jobs.pbs import create_vasp_job


class FrequencyWorkflow(BaseWorkflow):
    """
    Workflow for vibrational frequency calculations.

    Calculates frequencies, ZPE, and thermodynamic properties.

    Examples
    --------
    >>> # VASP frequency calculation
    >>> wf = FrequencyWorkflow(
    ...     relaxed_atoms,
    ...     work_dir="./freq",
    ...     selective_atoms=[10, 11, 12],  # Only NH3 atoms
    ... )
    >>> wf.setup()
    >>> # After VASP job completes:
    >>> results = wf.parse_results()
    >>> zpe = results["zpe"]
    >>> gibbs = wf.get_gibbs_energy(T=673, p=1.0)
    """

    def __init__(
        self,
        atoms: Atoms,
        work_dir: Union[str, Path] = "./",
        calculator: str = "vasp",
        selective_atoms: Optional[List[int]] = None,
        nfree: int = 2,
        potim: float = 0.015,
        # VASP parameters
        encut: float = 520,
        kspacing: float = 0.03,
        hubbard_u: Optional[Dict[str, float]] = None,
        vdw: Optional[str] = None,
        is_surface: bool = True,
        # PBS parameters
        nodes: int = 1,
        ppn: int = 24,
        walltime: str = "24:00:00",
        queue: Optional[str] = None,
        # Temperature for thermodynamics
        temperature: float = 298.15,
        pressure: float = 1.0,
        **kwargs,
    ):
        """
        Initialize FrequencyWorkflow.

        Parameters
        ----------
        atoms : Atoms
            Relaxed structure
        work_dir : str or Path
            Working directory
        calculator : str
            Calculator type: "vasp" or "mace"
        selective_atoms : list, optional
            Indices of atoms to include in frequency calculation.
            If None, uses atoms with constraints (or all if no constraints).
        nfree : int
            Number of displacements per atom direction (VASP NFREE)
        potim : float
            Displacement step size (VASP POTIM)
        encut : float
            VASP cutoff energy
        kspacing : float
            K-point spacing
        hubbard_u : dict, optional
            Hubbard U values
        vdw : str, optional
            VdW correction
        is_surface : bool
            Surface calculation settings
        nodes : int
            PBS nodes
        ppn : int
            Processors per node
        walltime : str
            PBS walltime
        queue : str, optional
            PBS queue
        temperature : float
            Temperature for thermodynamics (K)
        pressure : float
            Pressure for thermodynamics (bar)
        **kwargs : dict
            Additional parameters
        """
        super().__init__(atoms, work_dir, calculator, **kwargs)

        self.selective_atoms = selective_atoms
        self.nfree = nfree
        self.potim = potim

        # VASP settings
        self.encut = encut
        self.kspacing = kspacing
        self.hubbard_u = hubbard_u
        self.vdw = vdw
        self.is_surface = is_surface

        # PBS settings
        self.nodes = nodes
        self.ppn = ppn
        self.walltime = walltime
        self.queue = queue

        # Thermodynamics
        self.temperature = temperature
        self.pressure = pressure

        # Results storage
        self.frequencies = None
        self.freq_calc = None

    def _determine_selective_atoms(self) -> Optional[List[int]]:
        """Determine which atoms to include in frequency calculation."""
        if self.selective_atoms is not None:
            return self.selective_atoms

        # Check for constraints
        constraints = self.atoms.constraints
        if constraints:
            from ase.constraints import FixAtoms
            fixed_indices = set()
            for c in constraints:
                if isinstance(c, FixAtoms):
                    fixed_indices.update(c.index)

            # Free atoms are those not fixed
            all_indices = set(range(len(self.atoms)))
            free_indices = all_indices - fixed_indices

            if free_indices:
                return sorted(list(free_indices))

        return None

    def setup(self) -> Dict[str, Path]:
        """
        Set up the frequency calculation.

        Returns
        -------
        dict
            Dictionary of generated file paths
        """
        files = {}

        if self.calculator == "vasp":
            return self._setup_vasp()
        else:
            # For MACE, just save structure
            init_path = self.work_dir / "initial.traj"
            ase_write(str(init_path), self.atoms)
            files["initial"] = init_path
            return files

    def _setup_vasp(self) -> Dict[str, Path]:
        """Set up VASP frequency calculation."""
        selective = self._determine_selective_atoms()

        # Generate VASP inputs for frequency calculation
        gen = VASPInputGenerator(
            self.atoms,
            calc_type="frequency",
            work_dir=self.work_dir,
            NFREE=self.nfree,
            POTIM=self.potim,
        )

        # Set selective dynamics if needed
        if selective:
            gen.set_selective_dynamics(selective)

        vasp_files = gen.generate_all(
            encut=self.encut,
            kspacing=self.kspacing,
            hubbard_u=self.hubbard_u,
            vdw=self.vdw,
            is_surface=self.is_surface,
        )

        # Generate PBS script
        pbs_path = create_vasp_job(
            work_dir=self.work_dir,
            nodes=self.nodes,
            ppn=self.ppn,
            walltime=self.walltime,
            queue=self.queue,
        )

        vasp_files["PBS"] = pbs_path

        # Save initial structure
        init_path = self.work_dir / "initial.traj"
        ase_write(str(init_path), self.atoms)
        vasp_files["initial"] = init_path

        n_atoms = len(selective) if selective else len(self.atoms)
        print(f"Set up frequency calculation for {n_atoms} atoms")
        print(f"Submit with: qsub {pbs_path}")

        return vasp_files

    def run(self) -> Dict[str, Any]:
        """
        Run the frequency calculation.

        For VASP: returns instructions (manual submission required).
        For MACE: runs frequency calculation directly.

        Returns
        -------
        dict
            Results dictionary
        """
        if self.calculator == "vasp":
            pbs_path = self.work_dir / "run.pbs"
            return {
                "status": "setup_complete",
                "message": f"Submit job with: qsub {pbs_path}",
                "work_dir": str(self.work_dir),
            }

        elif self.calculator == "mace":
            return self._run_mace()

        else:
            raise ValueError(f"Unknown calculator: {self.calculator}")

    def _run_mace(self) -> Dict[str, Any]:
        """Run frequency calculation with MACE."""
        try:
            from mace.calculators import MACECalculator
        except ImportError:
            raise ImportError("MACE not installed")

        from ase.vibrations import Vibrations

        # Load model
        if hasattr(self, "mace_model") and self.mace_model:
            calc = MACECalculator(model_path=self.mace_model, device="cpu")
        else:
            try:
                from mace.calculators import mace_mp
                calc = mace_mp(model="medium", device="cpu")
            except Exception:
                raise ValueError("No MACE model specified")

        atoms = self.atoms.copy()
        atoms.calc = calc

        # Determine which atoms to vibrate
        selective = self._determine_selective_atoms()

        # Run vibrations
        vib = Vibrations(atoms, name=str(self.work_dir / "vib"), indices=selective)
        vib.run()

        # Get frequencies
        self.frequencies = vib.get_frequencies()

        # Calculate ZPE
        real_freqs = [f for f in self.frequencies if not np.isnan(f) and f > 0]
        zpe = 0.5 * sum(real_freqs) * H_EV_S * C_CMS * 100  # Convert to eV

        self.results = {
            "frequencies": self.frequencies,
            "frequencies_real": real_freqs,
            "n_imaginary": sum(1 for f in self.frequencies
                              if not np.isnan(f) and f < 0),
            "zpe": zpe,
        }

        # Save summary
        vib.summary(log=str(self.work_dir / "vibrations.log"))

        return self.results

    def parse_results(self) -> Dict[str, Any]:
        """
        Parse results from completed VASP frequency calculation.

        Returns
        -------
        dict
            Results dictionary
        """
        # Use FrequencyCalculation parser
        self.freq_calc = FrequencyCalculation(self.work_dir)

        self.frequencies = self.freq_calc.get_frequencies()
        real_freqs = [f for f in self.frequencies
                     if not np.isnan(f) and f.real > 0]

        self.results = {
            "frequencies": self.frequencies,
            "frequencies_real": [f.real for f in real_freqs],
            "n_imaginary": self.freq_calc.count_imaginary(),
            "zpe": self.freq_calc.get_zpe(),
            "electronic_energy": self.freq_calc.get_electronic_energy(),
        }

        return self.results

    def get_thermal_corrections(
        self,
        temperature: Optional[float] = None,
    ) -> Dict[str, float]:
        """
        Get thermal corrections at given temperature.

        Parameters
        ----------
        temperature : float, optional
            Temperature in K. Defaults to self.temperature.

        Returns
        -------
        dict
            Thermal corrections (ZPE, H, TS, etc.)
        """
        if self.freq_calc is None:
            self.parse_results()

        T = temperature if temperature is not None else self.temperature

        return {
            "zpe": self.freq_calc.get_zpe(),
            "thermal_correction": self.freq_calc.get_thermal_correction(T),
            "entropy": self.freq_calc.get_vibrational_entropy(T),
            "helmholtz": self.freq_calc.get_helmholtz_energy(T),
        }

    def get_gibbs_energy(
        self,
        temperature: Optional[float] = None,
        pressure: Optional[float] = None,
        electronic_energy: Optional[float] = None,
    ) -> float:
        """
        Calculate Gibbs free energy.

        G = E_elec + ZPE + H(T) - TS

        For adsorbates, uses harmonic approximation.

        Parameters
        ----------
        temperature : float, optional
            Temperature in K
        pressure : float, optional
            Pressure in bar (for gas molecules)
        electronic_energy : float, optional
            Electronic energy. If None, reads from calculation.

        Returns
        -------
        float
            Gibbs free energy in eV
        """
        if self.freq_calc is None:
            self.parse_results()

        T = temperature if temperature is not None else self.temperature

        # Electronic energy
        if electronic_energy is not None:
            E_elec = electronic_energy
        else:
            E_elec = self.freq_calc.get_electronic_energy()

        # Thermal corrections
        zpe = self.freq_calc.get_zpe()
        H_corr = self.freq_calc.get_thermal_correction(T)
        S = self.freq_calc.get_vibrational_entropy(T)

        # G = E + ZPE + H(T) - TS
        G = E_elec + zpe + H_corr - T * S

        return G

    def print_summary(self) -> None:
        """Print summary of frequency calculation results."""
        if self.results is None:
            self.parse_results()

        print("\n" + "=" * 60)
        print("Frequency Calculation Summary")
        print("=" * 60)

        print(f"\nNumber of frequencies: {len(self.results['frequencies_real'])}")
        print(f"Imaginary frequencies: {self.results['n_imaginary']}")

        print(f"\nZero-Point Energy: {self.results['zpe']:.4f} eV")

        if "electronic_energy" in self.results:
            print(f"Electronic Energy: {self.results['electronic_energy']:.4f} eV")

        # Thermal corrections at T
        T = self.temperature
        try:
            thermal = self.get_thermal_corrections(T)
            print(f"\nAt T = {T:.1f} K:")
            print(f"  Thermal correction: {thermal['thermal_correction']:.4f} eV")
            print(f"  Entropy (TS): {T * thermal['entropy']:.4f} eV")
            print(f"  Helmholtz energy: {thermal['helmholtz']:.4f} eV")
        except Exception as e:
            print(f"\nCould not calculate thermal corrections: {e}")

        print("=" * 60)


class GasPhaseThermo:
    """
    Thermodynamic properties for gas phase molecules.

    Uses ideal gas approximation with translational, rotational,
    and vibrational contributions.
    """

    # Reference energies for common gas molecules (DFT/PBE)
    REFERENCE_ENERGIES = {
        "H2": -6.77,      # eV (adjust based on your calculations)
        "N2": -16.63,
        "NH3": -19.54,
        "H2O": -14.22,
        "O2": -9.86,
    }

    # Molecular properties
    MOLECULAR_DATA = {
        "H2": {"mass": 2.016, "symmetry": 2, "linear": True, "spin": 0},
        "N2": {"mass": 28.01, "symmetry": 2, "linear": True, "spin": 0},
        "NH3": {"mass": 17.03, "symmetry": 3, "linear": False, "spin": 0},
        "H2O": {"mass": 18.02, "symmetry": 2, "linear": False, "spin": 0},
        "O2": {"mass": 32.00, "symmetry": 2, "linear": True, "spin": 1},
    }

    def __init__(
        self,
        molecule: str,
        electronic_energy: Optional[float] = None,
        frequencies: Optional[List[float]] = None,
    ):
        """
        Initialize GasPhaseThermo.

        Parameters
        ----------
        molecule : str
            Molecule name (H2, N2, NH3, etc.)
        electronic_energy : float, optional
            Electronic energy in eV. Uses reference if not provided.
        frequencies : list, optional
            Vibrational frequencies in cm^-1
        """
        self.molecule = molecule

        if electronic_energy is not None:
            self.E_elec = electronic_energy
        elif molecule in self.REFERENCE_ENERGIES:
            self.E_elec = self.REFERENCE_ENERGIES[molecule]
        else:
            raise ValueError(f"Unknown molecule: {molecule}. Provide electronic_energy.")

        self.frequencies = frequencies or []
        self.mol_data = self.MOLECULAR_DATA.get(molecule, {})

    def get_translational_entropy(self, T: float, p: float = 1.0) -> float:
        """
        Calculate translational entropy.

        Parameters
        ----------
        T : float
            Temperature in K
        p : float
            Pressure in bar

        Returns
        -------
        float
            Translational entropy in eV/K
        """
        from scipy.constants import k, h, N_A

        mass = self.mol_data.get("mass", 28.0) / 1000 / N_A  # kg per molecule
        V = k * T / (p * 1e5)  # Volume per molecule (m^3)

        # Sackur-Tetrode equation
        lambda_th = h / np.sqrt(2 * np.pi * mass * k * T)
        q_trans = V / lambda_th**3

        S_trans = k * (np.log(q_trans) + 5/2)

        return S_trans / EV_TO_J

    def get_rotational_entropy(self, T: float) -> float:
        """
        Calculate rotational entropy.

        Parameters
        ----------
        T : float
            Temperature in K

        Returns
        -------
        float
            Rotational entropy in eV/K
        """
        from scipy.constants import k

        sigma = self.mol_data.get("symmetry", 1)
        linear = self.mol_data.get("linear", False)

        # Approximate rotational contribution
        if linear:
            # Linear molecule: S = k * (ln(T/sigma*theta_rot) + 1)
            S_rot = k * (np.log(T / sigma) + 1)
        else:
            # Non-linear: S = k * (3/2 * ln(T/sigma) + 3/2)
            S_rot = k * (1.5 * np.log(T / sigma) + 1.5)

        return S_rot / EV_TO_J

    def get_vibrational_entropy(self, T: float) -> float:
        """
        Calculate vibrational entropy.

        Parameters
        ----------
        T : float
            Temperature in K

        Returns
        -------
        float
            Vibrational entropy in eV/K
        """
        if not self.frequencies:
            return 0.0

        from scipy.constants import k

        S_vib = 0.0
        for freq in self.frequencies:
            if freq <= 0:
                continue

            # Convert frequency to energy
            hv = freq * H_EV_S * C_CMS * 100  # eV
            x = hv / (KB_EV * T)

            if x > 100:
                continue  # Negligible contribution

            S_vib += k * (x / (np.exp(x) - 1) - np.log(1 - np.exp(-x)))

        return S_vib / EV_TO_J

    def get_zpe(self) -> float:
        """Get zero-point energy in eV."""
        if not self.frequencies:
            return 0.0

        zpe = 0.5 * sum(f for f in self.frequencies if f > 0)
        return zpe * H_EV_S * C_CMS * 100  # eV

    def get_gibbs_energy(
        self,
        T: float,
        p: float = 1.0,
    ) -> float:
        """
        Calculate Gibbs free energy.

        G = E_elec + ZPE + H(T) - T*S(T,p)

        Parameters
        ----------
        T : float
            Temperature in K
        p : float
            Pressure in bar

        Returns
        -------
        float
            Gibbs free energy in eV
        """
        # Zero-point energy
        zpe = self.get_zpe()

        # Enthalpy correction (approximately kT per degree of freedom)
        n_dof = 6 if not self.mol_data.get("linear", False) else 5
        H_corr = n_dof / 2 * KB_EV * T

        # Total entropy
        S_trans = self.get_translational_entropy(T, p)
        S_rot = self.get_rotational_entropy(T)
        S_vib = self.get_vibrational_entropy(T)
        S_total = S_trans + S_rot + S_vib

        # Gibbs energy
        G = self.E_elec + zpe + H_corr - T * S_total

        return G


class ThermochemistryWorkflow:
    """
    High-level workflow for thermochemistry calculations.

    Calculates Gibbs free energies for surfaces and adsorbates,
    and computes reaction thermodynamics.
    """

    def __init__(
        self,
        temperature: float = 298.15,
        pressure: float = 1.0,
    ):
        """
        Initialize ThermochemistryWorkflow.

        Parameters
        ----------
        temperature : float
            Temperature in K
        pressure : float
            Pressure in bar
        """
        self.temperature = temperature
        self.pressure = pressure

        self.energies = {}
        self.corrections = {}

    def add_surface(
        self,
        name: str,
        electronic_energy: float,
        zpe: float = 0.0,
    ) -> None:
        """
        Add surface energy.

        Parameters
        ----------
        name : str
            Surface identifier
        electronic_energy : float
            Electronic energy in eV
        zpe : float
            Zero-point energy correction
        """
        self.energies[name] = {
            "E_elec": electronic_energy,
            "zpe": zpe,
        }

    def add_adsorbate(
        self,
        name: str,
        electronic_energy: float,
        frequencies: Optional[List[float]] = None,
        zpe: Optional[float] = None,
    ) -> None:
        """
        Add adsorbate system energy.

        Parameters
        ----------
        name : str
            Adsorbate identifier (e.g., "NH3*")
        electronic_energy : float
            Electronic energy in eV
        frequencies : list, optional
            Vibrational frequencies in cm^-1
        zpe : float, optional
            Pre-calculated ZPE
        """
        if zpe is None and frequencies:
            zpe = 0.5 * sum(f for f in frequencies if f > 0) * H_EV_S * C_CMS * 100

        self.energies[name] = {
            "E_elec": electronic_energy,
            "frequencies": frequencies or [],
            "zpe": zpe or 0.0,
        }

    def get_gibbs_energy(
        self,
        name: str,
        T: Optional[float] = None,
    ) -> float:
        """
        Get Gibbs free energy for a system.

        Parameters
        ----------
        name : str
            System identifier
        T : float, optional
            Temperature

        Returns
        -------
        float
            Gibbs free energy in eV
        """
        if name not in self.energies:
            raise ValueError(f"Unknown system: {name}")

        T = T if T is not None else self.temperature
        data = self.energies[name]

        E_elec = data["E_elec"]
        zpe = data.get("zpe", 0.0)

        # Thermal correction (simplified harmonic)
        freqs = data.get("frequencies", [])
        H_corr = 0.0
        TS_corr = 0.0

        for freq in freqs:
            if freq <= 0:
                continue
            hv = freq * H_EV_S * C_CMS * 100
            x = hv / (KB_EV * T)
            if x < 100:
                H_corr += hv / (np.exp(x) - 1)
                TS_corr += KB_EV * T * (x / (np.exp(x) - 1) - np.log(1 - np.exp(-x)))

        G = E_elec + zpe + H_corr - TS_corr
        return G

    def calculate_reaction_energy(
        self,
        products: Dict[str, int],
        reactants: Dict[str, int],
        T: Optional[float] = None,
    ) -> float:
        """
        Calculate reaction Gibbs energy.

        ΔG = Σ(n_i * G_i)_products - Σ(n_i * G_i)_reactants

        Parameters
        ----------
        products : dict
            Product species and stoichiometry
        reactants : dict
            Reactant species and stoichiometry
        T : float, optional
            Temperature

        Returns
        -------
        float
            Reaction Gibbs energy in eV
        """
        T = T if T is not None else self.temperature

        G_products = sum(n * self.get_gibbs_energy(name, T)
                        for name, n in products.items())
        G_reactants = sum(n * self.get_gibbs_energy(name, T)
                         for name, n in reactants.items())

        return G_products - G_reactants

    def get_decomposition_energetics(
        self,
        pathway: Dict[str, float],
        T: Optional[float] = None,
    ) -> Dict[str, float]:
        """
        Calculate decomposition pathway energetics.

        Parameters
        ----------
        pathway : dict
            Dictionary of step energies
        T : float, optional
            Temperature

        Returns
        -------
        dict
            Step-by-step energetics
        """
        T = T if T is not None else self.temperature

        results = {}
        ref_energy = pathway.get("NH3*", 0.0)

        for step, E in pathway.items():
            results[step] = E - ref_energy

        return results


def run_frequency_calculation(
    atoms: Atoms,
    work_dir: Union[str, Path],
    calculator: str = "vasp",
    selective_atoms: Optional[List[int]] = None,
    **kwargs,
) -> FrequencyWorkflow:
    """
    Convenience function to set up a frequency calculation.

    Parameters
    ----------
    atoms : Atoms
        Relaxed structure
    work_dir : str or Path
        Working directory
    calculator : str
        Calculator type
    selective_atoms : list, optional
        Atoms to include in frequency calculation
    **kwargs : dict
        Additional parameters

    Returns
    -------
    FrequencyWorkflow
        Configured workflow
    """
    wf = FrequencyWorkflow(
        atoms,
        work_dir=work_dir,
        calculator=calculator,
        selective_atoms=selective_atoms,
        **kwargs,
    )
    wf.setup()
    return wf
