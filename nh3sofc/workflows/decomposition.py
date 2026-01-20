"""NH3 decomposition pathway workflow.

Automates the generation and calculation of NH3 decomposition intermediates:
NH3* → NH2* + H* → NH* + 2H* → N* + 3H*
"""

from pathlib import Path
from typing import Optional, List, Union, Dict, Any
import numpy as np
from ase import Atoms
from ase.io import write as ase_write, read as ase_read

from ..structure.decomposition import (
    DecompositionBuilder,
    generate_decomposition_pathway,
    save_decomposition_configs,
)
from .relaxation import RelaxationWorkflow, BatchRelaxation


class DecompositionWorkflow:
    """
    Workflow for NH3 decomposition pathway calculations.

    Generates intermediate configurations and sets up relaxation
    calculations for each step.

    Examples
    --------
    >>> # Generate and set up calculations
    >>> wf = DecompositionWorkflow(
    ...     nh3_on_slab,
    ...     work_dir="./decomposition",
    ...     n_configs_per_step=5,
    ... )
    >>> wf.setup()
    >>>
    >>> # After calculations complete:
    >>> results = wf.parse_results()
    >>> energies = wf.get_energy_profile()
    """

    STEPS = ["NH3", "NH2_H", "NH_2H", "N_3H"]

    def __init__(
        self,
        nh3_on_slab: Atoms,
        work_dir: Union[str, Path],
        n_configs_per_step: int = 5,
        calculator: str = "vasp",
        # Calculation parameters
        encut: float = 520,
        kspacing: float = 0.03,
        hubbard_u: Optional[Dict[str, float]] = None,
        vdw: Optional[str] = None,
        # PBS parameters
        nodes: int = 1,
        ppn: int = 24,
        walltime: str = "24:00:00",
        queue: Optional[str] = None,
        random_seed: Optional[int] = None,
        **kwargs,
    ):
        """
        Initialize DecompositionWorkflow.

        Parameters
        ----------
        nh3_on_slab : Atoms
            Optimized NH3 on surface structure
        work_dir : str or Path
            Base working directory
        n_configs_per_step : int
            Number of configurations per decomposition step
        calculator : str
            Calculator type ("vasp" or "mace")
        encut : float
            VASP cutoff energy
        kspacing : float
            K-point spacing
        hubbard_u : dict, optional
            Hubbard U values
        vdw : str, optional
            VdW correction
        nodes : int
            PBS nodes
        ppn : int
            Processors per node
        walltime : str
            PBS walltime
        queue : str, optional
            PBS queue
        random_seed : int, optional
            Random seed for configuration generation
        **kwargs : dict
            Additional parameters
        """
        self.nh3_on_slab = nh3_on_slab.copy()
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        self.n_configs_per_step = n_configs_per_step
        self.calculator = calculator
        self.random_seed = random_seed

        # Calculation parameters
        self.calc_params = {
            "encut": encut,
            "kspacing": kspacing,
            "hubbard_u": hubbard_u,
            "vdw": vdw,
            "nodes": nodes,
            "ppn": ppn,
            "walltime": walltime,
            "queue": queue,
            **kwargs,
        }

        # Storage for configurations and results
        self.pathway = {}
        self.batch_relaxations = {}
        self.results = {}

    def generate_configurations(self) -> Dict[str, List[Atoms]]:
        """
        Generate all decomposition intermediate configurations.

        Returns
        -------
        dict
            Dictionary mapping step names to lists of configurations
        """
        self.pathway = generate_decomposition_pathway(
            self.nh3_on_slab,
            n_configs_per_step=self.n_configs_per_step,
            random_seed=self.random_seed,
        )

        # Save configurations
        save_decomposition_configs(
            self.pathway,
            self.work_dir / "initial_configs",
        )

        # Print summary
        print("Generated decomposition configurations:")
        for step, configs in self.pathway.items():
            print(f"  {step}: {len(configs)} configurations")

        return self.pathway

    def setup(self) -> Dict[str, List[Dict[str, Path]]]:
        """
        Set up relaxation calculations for all configurations.

        Returns
        -------
        dict
            Dictionary of file paths for each step
        """
        if not self.pathway:
            self.generate_configurations()

        all_files = {}

        for step in self.STEPS:
            if step not in self.pathway:
                continue

            configs = self.pathway[step]
            step_dir = self.work_dir / step

            # Create batch relaxation
            batch = BatchRelaxation(
                structures=configs,
                base_dir=step_dir,
                calculator=self.calculator,
                naming_func=lambda i: f"config_{i:03d}",
                **self.calc_params,
            )

            files = batch.setup_all()
            all_files[step] = files
            self.batch_relaxations[step] = batch

        return all_files

    def parse_results(self) -> Dict[str, List[Dict[str, Any]]]:
        """
        Parse results from all completed calculations.

        Returns
        -------
        dict
            Results for each step
        """
        for step, batch in self.batch_relaxations.items():
            self.results[step] = batch.parse_all()

        return self.results

    def get_status(self) -> Dict[str, Dict[str, int]]:
        """
        Get status of all calculations.

        Returns
        -------
        dict
            Status for each step
        """
        status = {}
        for step, batch in self.batch_relaxations.items():
            status[step] = batch.get_status()
        return status

    def get_lowest_energies(self) -> Dict[str, float]:
        """
        Get lowest energy for each decomposition step.

        Returns
        -------
        dict
            Lowest energy for each step
        """
        if not self.results:
            self.parse_results()

        lowest = {}
        for step, results in self.results.items():
            energies = [r["energy"] for r in results
                       if "energy" in r and r["energy"] is not None]
            if energies:
                lowest[step] = min(energies)

        return lowest

    def get_energy_profile(self, reference: str = "NH3") -> Dict[str, float]:
        """
        Get relative energy profile.

        Parameters
        ----------
        reference : str
            Reference state (default: NH3)

        Returns
        -------
        dict
            Relative energies for each step
        """
        lowest = self.get_lowest_energies()

        if reference not in lowest:
            raise ValueError(f"Reference state {reference} not found")

        e_ref = lowest[reference]
        profile = {step: e - e_ref for step, e in lowest.items()}

        return profile

    def get_reaction_energies(self) -> Dict[str, float]:
        """
        Get step-by-step reaction energies.

        Returns
        -------
        dict
            Reaction energy for each step
        """
        profile = self.get_energy_profile()

        reactions = {}
        step_pairs = [
            ("NH3→NH2+H", "NH3", "NH2_H"),
            ("NH2+H→NH+2H", "NH2_H", "NH_2H"),
            ("NH+2H→N+3H", "NH_2H", "N_3H"),
        ]

        for name, initial, final in step_pairs:
            if initial in profile and final in profile:
                reactions[name] = profile[final] - profile[initial]

        return reactions

    def get_best_structures(self) -> Dict[str, Atoms]:
        """
        Get lowest energy structure for each step.

        Returns
        -------
        dict
            Best structure for each step
        """
        if not self.results:
            self.parse_results()

        best = {}
        for step, results in self.results.items():
            completed = [r for r in results
                        if "energy" in r and r["energy"] is not None]
            if completed:
                best_result = min(completed, key=lambda x: x["energy"])
                if "final_structure" in best_result:
                    best[step] = best_result["final_structure"]

        return best

    def print_summary(self) -> None:
        """Print summary of decomposition pathway results."""
        print("\n" + "=" * 60)
        print("NH3 Decomposition Pathway Summary")
        print("=" * 60)

        # Status
        status = self.get_status()
        print("\nCalculation Status:")
        for step, s in status.items():
            total = sum(s.values())
            print(f"  {step}: {s['completed']}/{total} completed")

        # Energy profile
        try:
            profile = self.get_energy_profile()
            print("\nEnergy Profile (relative to NH3*):")
            for step in self.STEPS:
                if step in profile:
                    print(f"  {step:10s}: {profile[step]:+.3f} eV")

            # Reaction energies
            reactions = self.get_reaction_energies()
            print("\nReaction Energies:")
            for name, dE in reactions.items():
                sign = "exo" if dE < 0 else "endo"
                print(f"  {name:20s}: {dE:+.3f} eV ({sign}thermic)")

            # Identify RDS (approximate)
            if reactions:
                rds = max(reactions.items(), key=lambda x: x[1])
                print(f"\nApproximate RDS: {rds[0]} (ΔE = {rds[1]:.3f} eV)")

        except Exception as e:
            print(f"\nCould not calculate energy profile: {e}")

        print("=" * 60)


class H2FormationWorkflow:
    """
    Workflow for H2 formation and desorption studies.

    Studies 2H* → H2* → H2(g) step.
    """

    def __init__(
        self,
        slab: Atoms,
        work_dir: Union[str, Path],
        h_h_distances: Optional[List[float]] = None,
        n_configs: int = 5,
        **kwargs,
    ):
        """
        Initialize H2FormationWorkflow.

        Parameters
        ----------
        slab : Atoms
            Clean slab structure
        work_dir : str or Path
            Working directory
        h_h_distances : list, optional
            H-H distances to sample
        n_configs : int
            Configurations per distance
        **kwargs : dict
            Calculation parameters
        """
        self.slab = slab.copy()
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        if h_h_distances is None:
            h_h_distances = [0.74, 1.0, 1.5, 2.0, 2.5, 3.0]
        self.h_h_distances = h_h_distances
        self.n_configs = n_configs
        self.kwargs = kwargs

        self.configurations = []
        self.batch = None

    def generate_configurations(self) -> List[Atoms]:
        """Generate 2H configurations at various distances."""
        from ..structure.decomposition import DecompositionBuilder

        # Create a dummy NH3 structure to use DecompositionBuilder
        # We'll create H2 configs directly instead
        z_top = self.slab.positions[:, 2].max()
        cell = self.slab.get_cell()
        a, b = cell[0, 0], cell[1, 1]

        configs = []

        for d in self.h_h_distances:
            for i in range(self.n_configs):
                slab_copy = self.slab.copy()

                # Random center
                cx = np.random.uniform(0, a)
                cy = np.random.uniform(0, b)
                angle = np.random.uniform(0, 2 * np.pi)

                # H positions
                h1_pos = [cx + d/2 * np.cos(angle),
                         cy + d/2 * np.sin(angle),
                         z_top + 1.0]
                h2_pos = [cx - d/2 * np.cos(angle),
                         cy - d/2 * np.sin(angle),
                         z_top + 1.0]

                # Add H atoms
                from ase import Atoms as AseAtoms
                h_atoms = AseAtoms("HH", positions=[h1_pos, h2_pos])
                combined = slab_copy + h_atoms

                configs.append(combined)

        self.configurations = configs
        return configs

    def setup(self) -> List[Dict[str, Path]]:
        """Set up relaxation calculations."""
        if not self.configurations:
            self.generate_configurations()

        self.batch = BatchRelaxation(
            structures=self.configurations,
            base_dir=self.work_dir,
            naming_func=lambda i: f"h2_config_{i:03d}",
            **self.kwargs,
        )

        return self.batch.setup_all()

    def parse_results(self) -> List[Dict[str, Any]]:
        """Parse calculation results."""
        return self.batch.parse_all()

    def get_h2_formation_energy(
        self,
        e_slab: float,
        e_h2_gas: float = -6.77,
    ) -> float:
        """
        Calculate H2 formation/desorption energy.

        ΔE = E(slab) + E(H2_gas) - E(2H*/slab)

        Parameters
        ----------
        e_slab : float
            Clean slab energy
        e_h2_gas : float
            H2 gas phase energy (default: -6.77 eV for PBE)

        Returns
        -------
        float
            H2 desorption energy
        """
        results = self.parse_results()
        energies = [r["energy"] for r in results
                   if "energy" in r and r["energy"] is not None]

        if not energies:
            raise ValueError("No completed calculations")

        e_2h_slab = min(energies)
        return e_slab + e_h2_gas - e_2h_slab


def run_decomposition_study(
    nh3_on_slab: Atoms,
    work_dir: Union[str, Path],
    n_configs: int = 5,
    calculator: str = "vasp",
    **kwargs,
) -> DecompositionWorkflow:
    """
    Convenience function to set up a complete decomposition study.

    Parameters
    ----------
    nh3_on_slab : Atoms
        Optimized NH3 on surface
    work_dir : str or Path
        Working directory
    n_configs : int
        Configurations per step
    calculator : str
        Calculator type
    **kwargs : dict
        Additional parameters

    Returns
    -------
    DecompositionWorkflow
        Configured workflow
    """
    wf = DecompositionWorkflow(
        nh3_on_slab,
        work_dir=work_dir,
        n_configs_per_step=n_configs,
        calculator=calculator,
        **kwargs,
    )
    wf.setup()
    return wf
