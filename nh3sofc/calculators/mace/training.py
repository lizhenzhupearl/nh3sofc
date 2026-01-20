"""MACE training data generation and management.

Provides tools for extracting training data from VASP calculations
and preparing datasets for MACE model training.
"""

from pathlib import Path
from typing import Optional, List, Dict, Any, Union
import json
import numpy as np
from ase import Atoms
from ase.io import read as ase_read, write as ase_write

from ...calculators.vasp.outputs import VASPOutputParser


class TrainingDataExtractor:
    """
    Extract training data from VASP calculations.

    Reads VASP output files and converts to MACE training format.

    Examples
    --------
    >>> extractor = TrainingDataExtractor()
    >>> extractor.add_from_vasp("./relax_calc")
    >>> extractor.add_from_vasp("./md_calc", max_frames=100)
    >>> extractor.save_xyz("training_data.xyz")
    """

    def __init__(self):
        """Initialize TrainingDataExtractor."""
        self.structures = []
        self.energies = []
        self.forces = []
        self.stresses = []
        self.metadata = []

    def add_structure(
        self,
        atoms: Atoms,
        energy: Optional[float] = None,
        forces: Optional[np.ndarray] = None,
        stress: Optional[np.ndarray] = None,
        source: str = "unknown",
    ) -> None:
        """
        Add a single structure to training data.

        Parameters
        ----------
        atoms : Atoms
            Structure
        energy : float, optional
            Total energy (eV)
        forces : array, optional
            Forces (eV/A)
        stress : array, optional
            Stress tensor
        source : str
            Data source identifier
        """
        self.structures.append(atoms.copy())
        self.energies.append(energy)
        self.forces.append(forces)
        self.stresses.append(stress)
        self.metadata.append({"source": source})

    def add_from_vasp(
        self,
        calc_dir: Union[str, Path],
        include_intermediate: bool = True,
        max_frames: Optional[int] = None,
    ) -> int:
        """
        Extract training data from VASP calculation.

        Parameters
        ----------
        calc_dir : str or Path
            VASP calculation directory
        include_intermediate : bool
            Include intermediate structures from optimization
        max_frames : int, optional
            Maximum number of frames to extract

        Returns
        -------
        int
            Number of structures added
        """
        calc_dir = Path(calc_dir)
        parser = VASPOutputParser(calc_dir)

        count = 0

        # Try to read trajectory from vasprun.xml
        vasprun = calc_dir / "vasprun.xml"
        if vasprun.exists() and include_intermediate:
            try:
                trajectory = ase_read(str(vasprun), index=":")

                if max_frames and len(trajectory) > max_frames:
                    # Sample evenly
                    indices = np.linspace(0, len(trajectory) - 1, max_frames, dtype=int)
                    trajectory = [trajectory[i] for i in indices]

                for atoms in trajectory:
                    energy = atoms.get_potential_energy() if atoms.calc else None
                    forces = atoms.get_forces() if atoms.calc else None

                    self.add_structure(
                        atoms,
                        energy=energy,
                        forces=forces,
                        source=str(calc_dir),
                    )
                    count += 1

            except Exception as e:
                print(f"Could not read vasprun.xml: {e}")

        # If no trajectory, try final structure
        if count == 0:
            try:
                atoms = parser.get_final_structure()
                if atoms:
                    energy = parser.get_energy()
                    forces = parser.get_forces()

                    self.add_structure(
                        atoms,
                        energy=energy,
                        forces=forces,
                        source=str(calc_dir),
                    )
                    count = 1
            except Exception as e:
                print(f"Could not read from {calc_dir}: {e}")

        return count

    def add_from_directory(
        self,
        base_dir: Union[str, Path],
        pattern: str = "**/OUTCAR",
        **kwargs,
    ) -> int:
        """
        Extract training data from multiple VASP calculations.

        Parameters
        ----------
        base_dir : str or Path
            Base directory to search
        pattern : str
            Glob pattern to find calculations
        **kwargs : dict
            Arguments passed to add_from_vasp

        Returns
        -------
        int
            Total number of structures added
        """
        base_dir = Path(base_dir)
        total = 0

        for outcar in base_dir.glob(pattern):
            calc_dir = outcar.parent
            count = self.add_from_vasp(calc_dir, **kwargs)
            total += count

        print(f"Added {total} structures from {base_dir}")
        return total

    def filter_by_energy(
        self,
        max_energy_per_atom: float = 0.0,
        reference_energy: Optional[float] = None,
    ) -> int:
        """
        Filter structures by energy.

        Parameters
        ----------
        max_energy_per_atom : float
            Maximum energy per atom relative to reference
        reference_energy : float, optional
            Reference energy per atom

        Returns
        -------
        int
            Number of structures removed
        """
        if reference_energy is None:
            # Use minimum energy as reference
            valid_energies = [
                e / len(s) for e, s in zip(self.energies, self.structures)
                if e is not None
            ]
            reference_energy = min(valid_energies) if valid_energies else 0.0

        keep_indices = []
        for i, (e, s) in enumerate(zip(self.energies, self.structures)):
            if e is None:
                continue
            e_per_atom = e / len(s)
            if e_per_atom - reference_energy <= max_energy_per_atom:
                keep_indices.append(i)

        removed = len(self.structures) - len(keep_indices)

        self.structures = [self.structures[i] for i in keep_indices]
        self.energies = [self.energies[i] for i in keep_indices]
        self.forces = [self.forces[i] for i in keep_indices]
        self.stresses = [self.stresses[i] for i in keep_indices]
        self.metadata = [self.metadata[i] for i in keep_indices]

        return removed

    def filter_by_forces(self, max_force: float = 10.0) -> int:
        """
        Filter structures by maximum force.

        Parameters
        ----------
        max_force : float
            Maximum allowed force (eV/A)

        Returns
        -------
        int
            Number of structures removed
        """
        keep_indices = []
        for i, forces in enumerate(self.forces):
            if forces is None:
                continue
            max_f = np.max(np.linalg.norm(forces, axis=1))
            if max_f <= max_force:
                keep_indices.append(i)

        removed = len(self.structures) - len(keep_indices)

        self.structures = [self.structures[i] for i in keep_indices]
        self.energies = [self.energies[i] for i in keep_indices]
        self.forces = [self.forces[i] for i in keep_indices]
        self.stresses = [self.stresses[i] for i in keep_indices]
        self.metadata = [self.metadata[i] for i in keep_indices]

        return removed

    def save_xyz(
        self,
        filename: str,
        include_forces: bool = True,
        include_stress: bool = False,
    ) -> None:
        """
        Save training data to extended XYZ format.

        Parameters
        ----------
        filename : str
            Output file path
        include_forces : bool
            Include forces in output
        include_stress : bool
            Include stress in output
        """
        atoms_list = []

        for i, atoms in enumerate(self.structures):
            atoms_copy = atoms.copy()

            # Add energy
            if self.energies[i] is not None:
                atoms_copy.info["energy"] = self.energies[i]
                atoms_copy.info["REF_energy"] = self.energies[i]

            # Add forces
            if include_forces and self.forces[i] is not None:
                atoms_copy.arrays["forces"] = self.forces[i]
                atoms_copy.arrays["REF_forces"] = self.forces[i]

            # Add stress
            if include_stress and self.stresses[i] is not None:
                atoms_copy.info["stress"] = self.stresses[i]
                atoms_copy.info["REF_stress"] = self.stresses[i]

            atoms_list.append(atoms_copy)

        ase_write(filename, atoms_list, format="extxyz")
        print(f"Saved {len(atoms_list)} structures to {filename}")

    def split_data(
        self,
        train_fraction: float = 0.8,
        val_fraction: float = 0.1,
        seed: int = 42,
    ) -> Dict[str, List[int]]:
        """
        Split data into train/val/test sets.

        Parameters
        ----------
        train_fraction : float
            Fraction for training
        val_fraction : float
            Fraction for validation
        seed : int
            Random seed

        Returns
        -------
        dict
            Indices for each split
        """
        np.random.seed(seed)
        n = len(self.structures)
        indices = np.random.permutation(n)

        n_train = int(n * train_fraction)
        n_val = int(n * val_fraction)

        return {
            "train": indices[:n_train].tolist(),
            "val": indices[n_train:n_train + n_val].tolist(),
            "test": indices[n_train + n_val:].tolist(),
        }

    def save_split(
        self,
        output_dir: Union[str, Path],
        train_fraction: float = 0.8,
        val_fraction: float = 0.1,
        seed: int = 42,
    ) -> Dict[str, Path]:
        """
        Save split datasets.

        Parameters
        ----------
        output_dir : str or Path
            Output directory
        train_fraction : float
            Training fraction
        val_fraction : float
            Validation fraction
        seed : int
            Random seed

        Returns
        -------
        dict
            Paths to saved files
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        splits = self.split_data(train_fraction, val_fraction, seed)
        paths = {}

        for split_name, indices in splits.items():
            if not indices:
                continue

            atoms_list = []
            for i in indices:
                atoms = self.structures[i].copy()
                if self.energies[i] is not None:
                    atoms.info["REF_energy"] = self.energies[i]
                if self.forces[i] is not None:
                    atoms.arrays["REF_forces"] = self.forces[i]
                atoms_list.append(atoms)

            filename = output_dir / f"{split_name}.xyz"
            ase_write(str(filename), atoms_list, format="extxyz")
            paths[split_name] = filename

            print(f"Saved {len(atoms_list)} structures to {filename}")

        return paths

    def get_statistics(self) -> Dict[str, Any]:
        """
        Get dataset statistics.

        Returns
        -------
        dict
            Statistics summary
        """
        valid_energies = [e for e in self.energies if e is not None]
        valid_forces = [f for f in self.forces if f is not None]

        stats = {
            "n_structures": len(self.structures),
            "n_with_energy": len(valid_energies),
            "n_with_forces": len(valid_forces),
        }

        if valid_energies:
            energies_per_atom = [
                e / len(s) for e, s in zip(self.energies, self.structures)
                if e is not None
            ]
            stats["energy_per_atom_mean"] = np.mean(energies_per_atom)
            stats["energy_per_atom_std"] = np.std(energies_per_atom)
            stats["energy_per_atom_range"] = (
                min(energies_per_atom),
                max(energies_per_atom),
            )

        if valid_forces:
            max_forces = [
                np.max(np.linalg.norm(f, axis=1)) for f in valid_forces
            ]
            stats["max_force_mean"] = np.mean(max_forces)
            stats["max_force_max"] = np.max(max_forces)

        # Composition statistics
        elements = set()
        for atoms in self.structures:
            elements.update(atoms.get_chemical_symbols())
        stats["elements"] = sorted(list(elements))

        return stats

    def print_summary(self) -> None:
        """Print dataset summary."""
        stats = self.get_statistics()

        print("\n" + "=" * 50)
        print("Training Data Summary")
        print("=" * 50)

        print(f"\nTotal structures: {stats['n_structures']}")
        print(f"With energy: {stats['n_with_energy']}")
        print(f"With forces: {stats['n_with_forces']}")

        if "energy_per_atom_mean" in stats:
            print(f"\nEnergy per atom:")
            print(f"  Mean: {stats['energy_per_atom_mean']:.4f} eV")
            print(f"  Std: {stats['energy_per_atom_std']:.4f} eV")
            print(f"  Range: [{stats['energy_per_atom_range'][0]:.4f}, "
                  f"{stats['energy_per_atom_range'][1]:.4f}] eV")

        if "max_force_mean" in stats:
            print(f"\nMaximum forces:")
            print(f"  Mean: {stats['max_force_mean']:.4f} eV/A")
            print(f"  Max: {stats['max_force_max']:.4f} eV/A")

        print(f"\nElements: {', '.join(stats['elements'])}")
        print("=" * 50)


class MACETrainingConfig:
    """
    Generate MACE training configuration.

    Creates configuration files for MACE training scripts.
    """

    DEFAULT_CONFIG = {
        "model": "MACE",
        "hidden_irreps": "128x0e + 128x1o",
        "r_max": 5.0,
        "batch_size": 10,
        "max_num_epochs": 500,
        "patience": 50,
        "lr": 0.01,
        "weight_decay": 5e-7,
        "energy_weight": 1.0,
        "forces_weight": 100.0,
        "stress_weight": 0.0,
    }

    def __init__(self, **kwargs):
        """
        Initialize MACETrainingConfig.

        Parameters
        ----------
        **kwargs : dict
            Configuration overrides
        """
        self.config = self.DEFAULT_CONFIG.copy()
        self.config.update(kwargs)

    def set_data_paths(
        self,
        train_file: str,
        valid_file: str,
        test_file: Optional[str] = None,
    ) -> None:
        """Set training data file paths."""
        self.config["train_file"] = train_file
        self.config["valid_file"] = valid_file
        if test_file:
            self.config["test_file"] = test_file

    def set_model_params(
        self,
        hidden_irreps: str = "128x0e + 128x1o",
        r_max: float = 5.0,
        num_interactions: int = 2,
    ) -> None:
        """Set model architecture parameters."""
        self.config["hidden_irreps"] = hidden_irreps
        self.config["r_max"] = r_max
        self.config["num_interactions"] = num_interactions

    def set_training_params(
        self,
        batch_size: int = 10,
        max_epochs: int = 500,
        lr: float = 0.01,
    ) -> None:
        """Set training parameters."""
        self.config["batch_size"] = batch_size
        self.config["max_num_epochs"] = max_epochs
        self.config["lr"] = lr

    def save(self, filename: str) -> None:
        """Save configuration to file."""
        with open(filename, "w") as f:
            json.dump(self.config, f, indent=2)
        print(f"Saved config to {filename}")

    def get_command(
        self,
        output_dir: str = "./mace_model",
        name: str = "nh3sofc_model",
    ) -> str:
        """
        Generate MACE training command.

        Parameters
        ----------
        output_dir : str
            Output directory for model
        name : str
            Model name

        Returns
        -------
        str
            Training command
        """
        cmd_parts = [
            "mace_run_train",
            f"--name={name}",
            f"--train_file={self.config.get('train_file', 'train.xyz')}",
            f"--valid_file={self.config.get('valid_file', 'valid.xyz')}",
            f"--model={self.config['model']}",
            f"--hidden_irreps='{self.config['hidden_irreps']}'",
            f"--r_max={self.config['r_max']}",
            f"--batch_size={self.config['batch_size']}",
            f"--max_num_epochs={self.config['max_num_epochs']}",
            f"--lr={self.config['lr']}",
            f"--energy_weight={self.config['energy_weight']}",
            f"--forces_weight={self.config['forces_weight']}",
            f"--results_dir={output_dir}",
        ]

        return " \\\n  ".join(cmd_parts)


def extract_training_data(
    directories: List[str],
    output_file: str = "training_data.xyz",
    **kwargs,
) -> TrainingDataExtractor:
    """
    Extract training data from multiple directories.

    Parameters
    ----------
    directories : list
        List of calculation directories
    output_file : str
        Output XYZ file
    **kwargs : dict
        Arguments for add_from_vasp

    Returns
    -------
    TrainingDataExtractor
        Extractor with data
    """
    extractor = TrainingDataExtractor()

    for directory in directories:
        extractor.add_from_vasp(directory, **kwargs)

    extractor.save_xyz(output_file)
    return extractor
