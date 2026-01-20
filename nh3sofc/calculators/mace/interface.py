"""MACE calculator interface for ML force field calculations.

Provides wrapper around MACE for ASE-compatible calculations
with support for foundation models and custom trained models.
"""

from pathlib import Path
from typing import Optional, Union, Dict, Any, List
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator


class MACECalculatorWrapper:
    """
    Wrapper for MACE calculator with NH3-SOFC specific features.

    Provides easy interface to MACE models including:
    - Foundation models (MACE-MP)
    - Custom trained models
    - Uncertainty estimation

    Examples
    --------
    >>> # Use foundation model
    >>> calc = MACECalculatorWrapper.from_foundation_model()
    >>> atoms.calc = calc.get_calculator()
    >>> energy = atoms.get_potential_energy()
    >>>
    >>> # Use custom model
    >>> calc = MACECalculatorWrapper("my_model.model")
    >>> atoms.calc = calc.get_calculator()
    """

    def __init__(
        self,
        model_path: Optional[str] = None,
        device: str = "cpu",
        default_dtype: str = "float64",
        dispersion: bool = False,
        dispersion_xc: str = "pbe",
    ):
        """
        Initialize MACECalculatorWrapper.

        Parameters
        ----------
        model_path : str, optional
            Path to MACE model file. If None, uses foundation model.
        device : str
            Device to run on ("cpu" or "cuda")
        default_dtype : str
            Default data type ("float32" or "float64")
        dispersion : bool
            Whether to add D3 dispersion correction
        dispersion_xc : str
            XC functional for dispersion (if dispersion=True)
        """
        self.model_path = model_path
        self.device = device
        self.default_dtype = default_dtype
        self.dispersion = dispersion
        self.dispersion_xc = dispersion_xc

        self._calculator = None
        self._mace_available = self._check_mace()

    def _check_mace(self) -> bool:
        """Check if MACE is available."""
        try:
            import mace
            return True
        except ImportError:
            return False

    @classmethod
    def from_foundation_model(
        cls,
        model: str = "medium",
        device: str = "cpu",
        dispersion: bool = True,
        **kwargs,
    ) -> "MACECalculatorWrapper":
        """
        Create calculator from MACE foundation model (MACE-MP).

        Parameters
        ----------
        model : str
            Model size: "small", "medium", or "large"
        device : str
            Device to run on
        dispersion : bool
            Whether to add D3 dispersion
        **kwargs : dict
            Additional parameters

        Returns
        -------
        MACECalculatorWrapper
            Calculator wrapper
        """
        wrapper = cls(
            model_path=None,
            device=device,
            dispersion=dispersion,
            **kwargs,
        )
        wrapper._foundation_model = model
        return wrapper

    @classmethod
    def from_custom_model(
        cls,
        model_path: str,
        device: str = "cpu",
        **kwargs,
    ) -> "MACECalculatorWrapper":
        """
        Create calculator from custom trained model.

        Parameters
        ----------
        model_path : str
            Path to model file
        device : str
            Device to run on
        **kwargs : dict
            Additional parameters

        Returns
        -------
        MACECalculatorWrapper
            Calculator wrapper
        """
        return cls(model_path=model_path, device=device, **kwargs)

    def get_calculator(self) -> Calculator:
        """
        Get ASE-compatible calculator.

        Returns
        -------
        Calculator
            MACE calculator for ASE
        """
        if not self._mace_available:
            raise ImportError(
                "MACE not installed. Install with: pip install mace-torch"
            )

        if self._calculator is not None:
            return self._calculator

        if self.model_path:
            # Custom model
            from mace.calculators import MACECalculator

            self._calculator = MACECalculator(
                model_paths=self.model_path,
                device=self.device,
                default_dtype=self.default_dtype,
            )
        else:
            # Foundation model
            from mace.calculators import mace_mp

            model_size = getattr(self, "_foundation_model", "medium")
            self._calculator = mace_mp(
                model=model_size,
                device=self.device,
                default_dtype=self.default_dtype,
                dispersion=self.dispersion,
            )

        return self._calculator

    def calculate(
        self,
        atoms: Atoms,
        properties: List[str] = ["energy", "forces"],
    ) -> Dict[str, Any]:
        """
        Calculate properties for a structure.

        Parameters
        ----------
        atoms : Atoms
            Structure to calculate
        properties : list
            Properties to calculate

        Returns
        -------
        dict
            Calculated properties
        """
        calc = self.get_calculator()
        atoms_copy = atoms.copy()
        atoms_copy.calc = calc

        results = {}

        if "energy" in properties:
            results["energy"] = atoms_copy.get_potential_energy()

        if "forces" in properties:
            results["forces"] = atoms_copy.get_forces()

        if "stress" in properties:
            try:
                results["stress"] = atoms_copy.get_stress()
            except Exception:
                pass

        return results

    def get_energy(self, atoms: Atoms) -> float:
        """Get potential energy."""
        return self.calculate(atoms, ["energy"])["energy"]

    def get_forces(self, atoms: Atoms) -> np.ndarray:
        """Get forces on atoms."""
        return self.calculate(atoms, ["forces"])["forces"]

    def relax(
        self,
        atoms: Atoms,
        fmax: float = 0.05,
        steps: int = 500,
        optimizer: str = "BFGS",
        trajectory: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Relax a structure using MACE.

        Parameters
        ----------
        atoms : Atoms
            Structure to relax
        fmax : float
            Force convergence criterion
        steps : int
            Maximum optimization steps
        optimizer : str
            Optimizer to use ("BFGS", "FIRE", "LBFGS")
        trajectory : str, optional
            Trajectory file path

        Returns
        -------
        dict
            Relaxation results
        """
        from ase.optimize import BFGS, FIRE, LBFGS

        optimizers = {
            "BFGS": BFGS,
            "FIRE": FIRE,
            "LBFGS": LBFGS,
        }

        if optimizer not in optimizers:
            raise ValueError(f"Unknown optimizer: {optimizer}")

        calc = self.get_calculator()
        atoms_copy = atoms.copy()
        atoms_copy.calc = calc

        opt_class = optimizers[optimizer]
        opt = opt_class(atoms_copy, trajectory=trajectory)
        opt.run(fmax=fmax, steps=steps)

        return {
            "converged": opt.converged(),
            "energy": atoms_copy.get_potential_energy(),
            "forces": atoms_copy.get_forces(),
            "max_force": np.max(np.linalg.norm(atoms_copy.get_forces(), axis=1)),
            "n_steps": opt.nsteps,
            "final_structure": atoms_copy,
        }


class MACEEnsemble:
    """
    Ensemble of MACE models for uncertainty estimation.

    Uses multiple models to estimate prediction uncertainty,
    useful for active learning.
    """

    def __init__(
        self,
        model_paths: List[str],
        device: str = "cpu",
    ):
        """
        Initialize MACEEnsemble.

        Parameters
        ----------
        model_paths : list
            Paths to ensemble model files
        device : str
            Device to run on
        """
        self.model_paths = model_paths
        self.device = device
        self.calculators = []

        self._load_models()

    def _load_models(self) -> None:
        """Load all models in ensemble."""
        try:
            from mace.calculators import MACECalculator
        except ImportError:
            raise ImportError("MACE not installed")

        for path in self.model_paths:
            calc = MACECalculator(model_paths=path, device=self.device)
            self.calculators.append(calc)

    def calculate_with_uncertainty(
        self,
        atoms: Atoms,
    ) -> Dict[str, Any]:
        """
        Calculate properties with uncertainty estimation.

        Parameters
        ----------
        atoms : Atoms
            Structure to calculate

        Returns
        -------
        dict
            Properties with uncertainties
        """
        energies = []
        forces_list = []

        for calc in self.calculators:
            atoms_copy = atoms.copy()
            atoms_copy.calc = calc

            energies.append(atoms_copy.get_potential_energy())
            forces_list.append(atoms_copy.get_forces())

        energies = np.array(energies)
        forces = np.array(forces_list)

        return {
            "energy_mean": np.mean(energies),
            "energy_std": np.std(energies),
            "forces_mean": np.mean(forces, axis=0),
            "forces_std": np.std(forces, axis=0),
            "max_force_std": np.max(np.std(np.linalg.norm(forces, axis=2), axis=0)),
        }

    def get_uncertainty(self, atoms: Atoms) -> float:
        """
        Get overall uncertainty for a structure.

        Parameters
        ----------
        atoms : Atoms
            Structure

        Returns
        -------
        float
            Uncertainty metric (energy std + max force std)
        """
        results = self.calculate_with_uncertainty(atoms)
        return results["energy_std"] + results["max_force_std"]

    def needs_retraining(
        self,
        atoms: Atoms,
        energy_threshold: float = 0.1,
        force_threshold: float = 0.5,
    ) -> bool:
        """
        Check if structure has high uncertainty (needs DFT).

        Parameters
        ----------
        atoms : Atoms
            Structure to check
        energy_threshold : float
            Energy uncertainty threshold (eV)
        force_threshold : float
            Force uncertainty threshold (eV/A)

        Returns
        -------
        bool
            True if structure needs DFT calculation
        """
        results = self.calculate_with_uncertainty(atoms)

        if results["energy_std"] > energy_threshold:
            return True
        if results["max_force_std"] > force_threshold:
            return True

        return False


def get_mace_calculator(
    model: str = "foundation",
    model_path: Optional[str] = None,
    device: str = "cpu",
    **kwargs,
) -> Calculator:
    """
    Get a MACE calculator.

    Parameters
    ----------
    model : str
        Model type: "foundation", "small", "medium", "large", or "custom"
    model_path : str, optional
        Path to custom model
    device : str
        Device to run on
    **kwargs : dict
        Additional parameters

    Returns
    -------
    Calculator
        MACE calculator
    """
    if model == "custom":
        if model_path is None:
            raise ValueError("model_path required for custom model")
        wrapper = MACECalculatorWrapper.from_custom_model(model_path, device, **kwargs)
    elif model in ["foundation", "small", "medium", "large"]:
        size = "medium" if model == "foundation" else model
        wrapper = MACECalculatorWrapper.from_foundation_model(size, device, **kwargs)
    else:
        raise ValueError(f"Unknown model type: {model}")

    return wrapper.get_calculator()
