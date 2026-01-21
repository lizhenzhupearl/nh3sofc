"""Exsolution simulation workflow.

Provides workflows for studying exsolution processes in perovskite materials:
1. Defective perovskite stability
2. Surface segregation energetics
3. Nanoparticle formation energies
4. Integration with NH3 decomposition studies
"""

from pathlib import Path
from typing import Optional, List, Union, Dict, Any
import numpy as np
from ase import Atoms
from ase.io import write as ase_write, read as ase_read

from ..core.base import BaseWorkflow
from ..core.constants import EXSOLUTION_METALS, DEFAULT_EXSOLUTION_PARAMS
from ..structure.exsolution import ExsolutionBuilder, ExsolutionStructure
from ..calculators.vasp.inputs import VASPInputGenerator
from ..calculators.vasp.outputs import VASPOutputParser
from ..jobs.pbs import create_vasp_job


class ExsolutionWorkflow(BaseWorkflow):
    """
    Workflow for studying exsolution processes.

    Generates and manages calculations for the complete exsolution pathway:
    1. Pristine perovskite surface
    2. Defective perovskite (with vacancies)
    3. Surface-segregated metal
    4. Exsolved nanoparticle on surface

    Examples
    --------
    >>> from nh3sofc.structure import BulkStructure, SurfaceBuilder
    >>>
    >>> # Prepare perovskite surface
    >>> bulk = BulkStructure.from_cif("LaSrTiNiO3.cif")
    >>> surface = SurfaceBuilder(bulk).create_surface(
    ...     miller_index=(0, 0, 1),
    ...     layers=6,
    ...     vacuum=15.0
    ... )
    >>>
    >>> # Set up exsolution workflow
    >>> wf = ExsolutionWorkflow(
    ...     atoms=surface.atoms,
    ...     work_dir="./exsolution_study",
    ...     metal="Ni",
    ...     particle_size=13,
    ...     calculator="vasp",
    ... )
    >>> wf.generate_pathway_structures()
    >>> wf.setup()
    >>>
    >>> # After VASP calculations complete:
    >>> results = wf.parse_results()
    >>> print(f"Exsolution energy: {results['exsolution_energy']:.2f} eV")
    """

    def __init__(
        self,
        atoms: Atoms,
        work_dir: Union[str, Path],
        metal: str = "Ni",
        particle_size: int = 13,
        vacancy_fraction: float = 0.1,
        calculator: str = "vasp",
        # VASP parameters
        encut: float = 520,
        kspacing: float = 0.03,
        hubbard_u: Optional[Dict[str, float]] = None,
        vdw: Optional[str] = "D3BJ",
        # PBS parameters
        nodes: int = 2,
        ppn: int = 24,
        walltime: str = "48:00:00",
        queue: Optional[str] = None,
        **kwargs,
    ):
        """
        Initialize ExsolutionWorkflow.

        Parameters
        ----------
        atoms : Atoms
            Base perovskite surface structure
        work_dir : str or Path
            Working directory for calculations
        metal : str
            Exsolution metal (Ni, Co, Fe)
        particle_size : int
            Number of atoms in nanoparticle (1, 4, 7, 13, 19)
        vacancy_fraction : float
            Oxygen vacancy concentration (0.0 to 1.0)
        calculator : str
            Calculator type ("vasp" or "mace")
        encut : float
            VASP cutoff energy
        kspacing : float
            K-point spacing
        hubbard_u : dict, optional
            Hubbard U values for transition metals
        vdw : str, optional
            VdW correction method
        nodes : int
            PBS nodes per calculation
        ppn : int
            Processors per node
        walltime : str
            PBS walltime
        queue : str, optional
            PBS queue name
        **kwargs : dict
            Additional parameters
        """
        super().__init__(atoms, work_dir, calculator, **kwargs)

        self.metal = metal
        self.particle_size = particle_size
        self.vacancy_fraction = vacancy_fraction

        # VASP settings
        self.encut = encut
        self.kspacing = kspacing
        self.hubbard_u = hubbard_u or {}
        self.vdw = vdw

        # Add Hubbard U for exsolution metal if not specified
        if metal not in self.hubbard_u:
            from ..core.constants import HUBBARD_U
            if metal in HUBBARD_U:
                self.hubbard_u[metal] = HUBBARD_U[metal]

        # PBS settings
        self.nodes = nodes
        self.ppn = ppn
        self.walltime = walltime
        self.queue = queue

        # Pathway structures
        self.pathway = []
        self.stage_dirs = {}

    def generate_pathway_structures(self) -> List[Dict[str, Any]]:
        """
        Generate structures for all exsolution pathway stages.

        Creates:
        1. Pristine perovskite surface
        2. Defective perovskite (with vacancies)
        3. Surface-segregated metal
        4. Exsolved nanoparticle

        Returns
        -------
        list
            List of pathway step dictionaries
        """
        builder = ExsolutionBuilder(self.atoms)

        self.pathway = builder.create_exsolution_pathway(
            metal=self.metal,
            particle_size=self.particle_size,
            vacancy_fraction=self.vacancy_fraction,
        )

        print(f"Generated {len(self.pathway)} exsolution pathway structures:")
        for step in self.pathway:
            print(f"  - {step['stage']}: {step['description']}")

        return self.pathway

    def setup(self) -> Dict[str, Dict[str, Path]]:
        """
        Set up calculations for all pathway stages.

        Generates VASP input files and PBS scripts for each stage.

        Returns
        -------
        dict
            Dictionary of stage names to file paths
        """
        if not self.pathway:
            self.generate_pathway_structures()

        all_files = {}

        for step in self.pathway:
            stage = step["stage"]
            atoms = step["atoms"]

            # Create stage directory
            stage_dir = self.work_dir / stage
            stage_dir.mkdir(parents=True, exist_ok=True)
            self.stage_dirs[stage] = stage_dir

            if self.calculator == "vasp":
                # Generate VASP inputs
                gen = VASPInputGenerator(
                    atoms,
                    calc_type="relax",
                    work_dir=stage_dir,
                    NSW=300,
                    EDIFFG=-0.02,
                )

                files = gen.generate_all(
                    encut=self.encut,
                    kspacing=self.kspacing,
                    hubbard_u=self.hubbard_u,
                    vdw=self.vdw,
                    is_surface=True,
                )

                # Generate PBS script
                pbs_path = create_vasp_job(
                    work_dir=stage_dir,
                    nodes=self.nodes,
                    ppn=self.ppn,
                    walltime=self.walltime,
                    queue=self.queue,
                )
                files["PBS"] = pbs_path

                # Save initial structure
                ase_write(str(stage_dir / "initial.vasp"), atoms, format="vasp")

                all_files[stage] = files

        # Create submission script for all stages
        self._create_batch_script()

        print(f"\nSetup complete. To run all calculations:")
        print(f"  cd {self.work_dir}")
        print(f"  bash submit_all.sh")

        return all_files

    def _create_batch_script(self) -> Path:
        """Create batch submission script for all stages."""
        script_path = self.work_dir / "submit_all.sh"

        lines = ["#!/bin/bash", "", "# Submit all exsolution pathway calculations", ""]

        for stage in self.stage_dirs:
            lines.append(f"echo 'Submitting {stage}...'")
            lines.append(f"cd {stage}")
            lines.append("qsub run.pbs")
            lines.append("cd ..")
            lines.append("")

        lines.append("echo 'All jobs submitted.'")

        with open(script_path, "w") as f:
            f.write("\n".join(lines))

        script_path.chmod(0o755)
        return script_path

    def run(self) -> Dict[str, Any]:
        """
        Run or provide instructions for running calculations.

        For VASP: Returns submission instructions.
        For MACE: Runs calculations directly.

        Returns
        -------
        dict
            Results or instructions
        """
        if self.calculator == "vasp":
            return {
                "status": "setup_complete",
                "message": f"Submit jobs with: cd {self.work_dir} && bash submit_all.sh",
                "stages": list(self.stage_dirs.keys()),
            }
        elif self.calculator == "mace":
            return self._run_mace()
        else:
            raise ValueError(f"Unknown calculator: {self.calculator}")

    def _run_mace(self) -> Dict[str, Any]:
        """Run optimizations with MACE calculator."""
        try:
            from mace.calculators import mace_mp
        except ImportError:
            raise ImportError("MACE not installed. Install with: pip install mace-torch")

        calc = mace_mp(model="medium", device="cpu")
        results = {}

        for step in self.pathway:
            stage = step["stage"]
            atoms = step["atoms"].copy()
            atoms.calc = calc

            stage_dir = self.work_dir / stage
            stage_dir.mkdir(parents=True, exist_ok=True)

            from ase.optimize import BFGS
            opt = BFGS(atoms, trajectory=str(stage_dir / "opt.traj"))
            opt.run(fmax=0.02, steps=200)

            results[stage] = {
                "energy": atoms.get_potential_energy(),
                "converged": opt.converged(),
                "n_steps": opt.nsteps,
            }

            ase_write(str(stage_dir / "final.vasp"), atoms, format="vasp")

        self.results = results
        self._calculate_exsolution_energetics(results)
        return self.results

    def parse_results(self) -> Dict[str, Any]:
        """
        Parse results from completed VASP calculations.

        Returns
        -------
        dict
            Results dictionary with energies and exsolution analysis
        """
        results = {}

        for stage, stage_dir in self.stage_dirs.items():
            try:
                parser = VASPOutputParser(stage_dir)
                results[stage] = {
                    "energy": parser.get_energy(),
                    "converged": parser.is_converged(),
                    "max_force": parser.get_max_force(),
                    "final_structure": parser.get_final_structure(),
                }
            except Exception as e:
                results[stage] = {"error": str(e)}

        self.results = results
        self._calculate_exsolution_energetics(results)

        return self.results

    def _calculate_exsolution_energetics(self, results: Dict[str, Any]) -> None:
        """Calculate exsolution energetics from stage energies."""
        if "pristine" not in results or "exsolved" not in results:
            return

        e_pristine = results.get("pristine", {}).get("energy")
        e_defective = results.get("defective", {}).get("energy")
        e_segregated = results.get("segregated", {}).get("energy")
        e_exsolved = results.get("exsolved", {}).get("energy")

        if e_pristine is not None and e_exsolved is not None:
            # Get bulk metal energy
            metal_info = EXSOLUTION_METALS.get(self.metal, {})
            e_bulk_metal = metal_info.get("bulk_energy", -5.5) * self.particle_size

            # Exsolution energy (simplified, without O chemical potential correction)
            e_exsolution = e_exsolved - e_pristine - e_bulk_metal
            self.results["exsolution_energy"] = e_exsolution

            # Segregation energy (if available)
            if e_defective is not None and e_segregated is not None:
                e_segregation = e_segregated - e_defective
                self.results["segregation_energy"] = e_segregation

            # Summary
            self.results["summary"] = {
                "metal": self.metal,
                "particle_size": self.particle_size,
                "vacancy_fraction": self.vacancy_fraction,
                "exsolution_energy_eV": e_exsolution,
                "favorable": e_exsolution < 0,
            }

    def couple_with_decomposition(
        self,
        decomposition_kwargs: Optional[Dict[str, Any]] = None,
    ) -> "DecompositionWorkflow":
        """
        Set up NH3 decomposition study on exsolved particle surface.

        Parameters
        ----------
        decomposition_kwargs : dict, optional
            Parameters for DecompositionWorkflow

        Returns
        -------
        DecompositionWorkflow
            Configured decomposition workflow
        """
        from .decomposition import DecompositionWorkflow

        # Get exsolved structure
        if not self.pathway:
            self.generate_pathway_structures()

        exsolved_step = next(
            (s for s in self.pathway if s["stage"] == "exsolved"),
            None
        )

        if exsolved_step is None:
            raise ValueError("No exsolved structure in pathway")

        exsolved_atoms = exsolved_step["atoms"]

        # Create decomposition workflow
        kwargs = decomposition_kwargs or {}
        decomp_dir = self.work_dir / "nh3_decomposition"

        decomp_wf = DecompositionWorkflow(
            atoms=exsolved_atoms,
            work_dir=decomp_dir,
            calculator=self.calculator,
            **kwargs,
        )

        print(f"Created NH3 decomposition workflow at {decomp_dir}")
        return decomp_wf


class ExsolutionScreeningWorkflow:
    """
    High-throughput screening workflow for exsolution parameters.

    Screens multiple:
    - Exsolution metals (Ni, Co, Fe)
    - Particle sizes
    - Vacancy concentrations

    Examples
    --------
    >>> screening = ExsolutionScreeningWorkflow(
    ...     base_structure=perovskite_surface,
    ...     parameter_space={
    ...         "metal": ["Ni", "Co", "Fe"],
    ...         "particle_size": [1, 4, 13],
    ...         "vacancy_fraction": [0.0, 0.05, 0.1],
    ...     },
    ...     work_dir="./exsolution_screening",
    ... )
    >>> screening.generate_all()
    >>> screening.setup_all()
    >>>
    >>> # After calculations:
    >>> results = screening.parse_all()
    >>> best = screening.get_best_result()
    """

    def __init__(
        self,
        base_structure: Atoms,
        parameter_space: Dict[str, List[Any]],
        work_dir: Union[str, Path],
        calculator: str = "vasp",
        n_configs_per_combo: int = 1,
        **kwargs,
    ):
        """
        Initialize ExsolutionScreeningWorkflow.

        Parameters
        ----------
        base_structure : Atoms
            Base perovskite surface structure
        parameter_space : dict
            Dictionary of parameter lists to screen:
            - "metal": list of metals (e.g., ["Ni", "Co"])
            - "particle_size": list of sizes (e.g., [1, 13])
            - "vacancy_fraction": list of fractions (e.g., [0.0, 0.1])
        work_dir : str or Path
            Base working directory
        calculator : str
            Calculator type
        n_configs_per_combo : int
            Number of configurations per parameter combination
        **kwargs : dict
            Additional parameters passed to ExsolutionWorkflow
        """
        self.base_structure = base_structure
        self.parameter_space = parameter_space
        self.work_dir = Path(work_dir)
        self.calculator = calculator
        self.n_configs = n_configs_per_combo
        self.kwargs = kwargs

        self.work_dir.mkdir(parents=True, exist_ok=True)

        self.workflows = []
        self.configurations = []

    def generate_all(self) -> List[Dict[str, Any]]:
        """
        Generate all parameter combinations.

        Returns
        -------
        list
            List of configuration dictionaries
        """
        metals = self.parameter_space.get("metal", ["Ni"])
        sizes = self.parameter_space.get("particle_size", [13])
        vac_fracs = self.parameter_space.get("vacancy_fraction", [0.1])

        self.configurations = []
        config_id = 0

        for metal in metals:
            for size in sizes:
                for vac_frac in vac_fracs:
                    for i in range(self.n_configs):
                        config = {
                            "config_id": config_id,
                            "metal": metal,
                            "particle_size": size,
                            "vacancy_fraction": vac_frac,
                            "config_index": i,
                            "name": f"{metal}_n{size}_vac{vac_frac:.2f}_{i:02d}",
                        }
                        self.configurations.append(config)
                        config_id += 1

        print(f"Generated {len(self.configurations)} configurations to screen")
        return self.configurations

    def setup_all(self) -> Dict[str, Dict[str, Path]]:
        """
        Set up calculations for all configurations.

        Returns
        -------
        dict
            Dictionary of configuration names to file paths
        """
        if not self.configurations:
            self.generate_all()

        all_files = {}

        for config in self.configurations:
            name = config["name"]
            config_dir = self.work_dir / name

            wf = ExsolutionWorkflow(
                atoms=self.base_structure,
                work_dir=config_dir,
                metal=config["metal"],
                particle_size=config["particle_size"],
                vacancy_fraction=config["vacancy_fraction"],
                calculator=self.calculator,
                **self.kwargs,
            )

            wf.generate_pathway_structures()
            files = wf.setup()
            all_files[name] = files
            self.workflows.append(wf)

        print(f"Set up {len(self.workflows)} exsolution workflows in {self.work_dir}")

        # Create master submission script
        self._create_master_script()

        return all_files

    def _create_master_script(self) -> Path:
        """Create master script to submit all configurations."""
        script_path = self.work_dir / "submit_all_screening.sh"

        lines = [
            "#!/bin/bash",
            "",
            "# Submit all exsolution screening calculations",
            "",
        ]

        for config in self.configurations:
            lines.append(f"cd {config['name']}")
            lines.append("bash submit_all.sh")
            lines.append("cd ..")
            lines.append("")

        with open(script_path, "w") as f:
            f.write("\n".join(lines))

        script_path.chmod(0o755)
        return script_path

    def parse_all(self) -> List[Dict[str, Any]]:
        """
        Parse results from all completed calculations.

        Returns
        -------
        list
            List of results for each configuration
        """
        results = []

        for wf, config in zip(self.workflows, self.configurations):
            try:
                result = wf.parse_results()
                result["config"] = config
                results.append(result)
            except Exception as e:
                results.append({
                    "config": config,
                    "error": str(e),
                })

        # Sort by exsolution energy
        valid_results = [
            r for r in results
            if "exsolution_energy" in r and r["exsolution_energy"] is not None
        ]
        if valid_results:
            valid_results.sort(key=lambda x: x["exsolution_energy"])

        return results

    def get_best_result(
        self,
        metric: str = "exsolution_energy",
        minimize: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """
        Get the best result based on specified metric.

        Parameters
        ----------
        metric : str
            Metric to optimize ("exsolution_energy", etc.)
        minimize : bool
            If True, find minimum; if False, find maximum

        Returns
        -------
        dict or None
            Best result dictionary
        """
        results = self.parse_all()
        valid = [r for r in results if metric in r and r[metric] is not None]

        if not valid:
            return None

        if minimize:
            return min(valid, key=lambda x: x[metric])
        else:
            return max(valid, key=lambda x: x[metric])


def run_exsolution_study(
    atoms: Atoms,
    work_dir: Union[str, Path],
    metal: str = "Ni",
    particle_size: int = 13,
    calculator: str = "vasp",
    **kwargs,
) -> Dict[str, Any]:
    """
    Convenience function to set up a single exsolution study.

    Parameters
    ----------
    atoms : Atoms
        Perovskite surface structure
    work_dir : str or Path
        Working directory
    metal : str
        Exsolution metal
    particle_size : int
        Nanoparticle size
    calculator : str
        Calculator type
    **kwargs : dict
        Additional workflow parameters

    Returns
    -------
    dict
        Setup results
    """
    wf = ExsolutionWorkflow(
        atoms=atoms,
        work_dir=work_dir,
        metal=metal,
        particle_size=particle_size,
        calculator=calculator,
        **kwargs,
    )
    wf.generate_pathway_structures()
    wf.setup()
    return wf.run()
