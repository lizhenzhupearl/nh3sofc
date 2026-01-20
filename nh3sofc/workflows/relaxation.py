"""Geometry optimization workflow.

Provides automated relaxation using VASP or MACE calculators.
"""

from pathlib import Path
from typing import Optional, List, Union, Dict, Any
import numpy as np
from ase import Atoms
from ase.io import write as ase_write, read as ase_read
from ase.optimize import BFGS, FIRE, LBFGS

from ..core.base import BaseWorkflow
from ..calculators.vasp.inputs import VASPInputGenerator
from ..calculators.vasp.outputs import VASPOutputParser
from ..jobs.pbs import VASPJobScript, create_vasp_job


class RelaxationWorkflow(BaseWorkflow):
    """
    Workflow for geometry optimization.

    Supports VASP (via job submission) or MACE (direct optimization).

    Examples
    --------
    >>> # VASP relaxation
    >>> wf = RelaxationWorkflow(
    ...     atoms,
    ...     work_dir="./relax",
    ...     calculator="vasp",
    ...     encut=520,
    ...     hubbard_u={"V": 3.25},
    ... )
    >>> wf.setup()
    >>> # Submit job manually, then:
    >>> results = wf.parse_results()
    >>>
    >>> # MACE relaxation (direct)
    >>> wf = RelaxationWorkflow(atoms, calculator="mace")
    >>> results = wf.run()
    """

    def __init__(
        self,
        atoms: Atoms,
        work_dir: Union[str, Path] = "./",
        calculator: str = "vasp",
        fmax: float = 0.02,
        steps: int = 300,
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
        # MACE parameters
        mace_model: Optional[str] = None,
        **kwargs,
    ):
        """
        Initialize RelaxationWorkflow.

        Parameters
        ----------
        atoms : Atoms
            Input structure
        work_dir : str or Path
            Working directory
        calculator : str
            Calculator type: "vasp" or "mace"
        fmax : float
            Force convergence criterion (eV/A)
        steps : int
            Maximum optimization steps
        encut : float
            VASP cutoff energy
        kspacing : float
            K-point spacing
        hubbard_u : dict, optional
            Hubbard U values
        vdw : str, optional
            VdW correction method
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
        mace_model : str, optional
            Path to MACE model
        **kwargs : dict
            Additional parameters
        """
        super().__init__(atoms, work_dir, calculator, **kwargs)

        self.fmax = fmax
        self.steps = steps

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

        # MACE settings
        self.mace_model = mace_model

    def setup(self) -> Dict[str, Path]:
        """
        Set up the relaxation calculation.

        For VASP: generates input files and PBS script.
        For MACE: no setup needed.

        Returns
        -------
        dict
            Dictionary of generated file paths
        """
        files = {}

        if self.calculator == "vasp":
            # Generate VASP inputs
            gen = VASPInputGenerator(
                self.atoms,
                calc_type="relax",
                work_dir=self.work_dir,
                NSW=self.steps,
                EDIFFG=-self.fmax,
            )

            vasp_files = gen.generate_all(
                encut=self.encut,
                kspacing=self.kspacing,
                hubbard_u=self.hubbard_u,
                vdw=self.vdw,
                is_surface=self.is_surface,
            )
            files.update(vasp_files)

            # Generate PBS script
            pbs_path = create_vasp_job(
                work_dir=self.work_dir,
                nodes=self.nodes,
                ppn=self.ppn,
                walltime=self.walltime,
                queue=self.queue,
            )
            files["PBS"] = pbs_path

            # Save initial structure
            init_path = self.work_dir / "initial.traj"
            ase_write(str(init_path), self.atoms)
            files["initial"] = init_path

        return files

    def run(self) -> Dict[str, Any]:
        """
        Run the relaxation.

        For VASP: returns instructions (manual submission required).
        For MACE: runs optimization directly.

        Returns
        -------
        dict
            Results dictionary
        """
        if self.calculator == "vasp":
            # VASP requires manual job submission
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
        """Run optimization with MACE calculator."""
        try:
            from mace.calculators import MACECalculator
        except ImportError:
            raise ImportError("MACE not installed. Install with: pip install mace-torch")

        # Load MACE model
        if self.mace_model is None:
            # Try to use foundation model
            try:
                from mace.calculators import mace_mp
                calc = mace_mp(model="medium", device="cpu")
            except Exception:
                raise ValueError("No MACE model specified and foundation model not available")
        else:
            calc = MACECalculator(model_path=self.mace_model, device="cpu")

        # Run optimization
        atoms = self.atoms.copy()
        atoms.calc = calc

        optimizer = BFGS(atoms, trajectory=str(self.work_dir / "opt.traj"))
        optimizer.run(fmax=self.fmax, steps=self.steps)

        # Get results
        self.results = {
            "converged": optimizer.converged(),
            "energy": atoms.get_potential_energy(),
            "forces": atoms.get_forces(),
            "max_force": np.max(np.linalg.norm(atoms.get_forces(), axis=1)),
            "n_steps": optimizer.nsteps,
            "final_structure": atoms.copy(),
        }

        # Save final structure
        ase_write(str(self.work_dir / "final.traj"), atoms)
        ase_write(str(self.work_dir / "CONTCAR"), atoms, format="vasp")

        return self.results

    def parse_results(self) -> Dict[str, Any]:
        """
        Parse results from completed VASP calculation.

        Returns
        -------
        dict
            Results dictionary
        """
        parser = VASPOutputParser(self.work_dir)

        self.results = {
            "converged": parser.is_converged(),
            "energy": parser.get_energy(),
            "energy_per_atom": parser.get_energy(per_atom=True),
            "max_force": parser.get_max_force(),
            "n_ionic_steps": parser.get_n_ionic_steps(),
            "final_structure": parser.get_final_structure(),
            "magnetic_moment": parser.get_magnetic_moment(),
            "elapsed_time": parser.get_elapsed_time(),
        }

        return self.results

    def is_complete(self) -> bool:
        """Check if calculation is complete."""
        if self.calculator == "vasp":
            outcar = self.work_dir / "OUTCAR"
            if outcar.exists():
                parser = VASPOutputParser(self.work_dir)
                return parser.is_converged()
        return False


class BatchRelaxation:
    """
    Batch relaxation of multiple structures.

    Examples
    --------
    >>> batch = BatchRelaxation(
    ...     structures=configs,
    ...     base_dir="./relaxations",
    ...     calculator="vasp",
    ... )
    >>> batch.setup_all()
    >>> # After jobs complete:
    >>> results = batch.parse_all()
    """

    def __init__(
        self,
        structures: List[Atoms],
        base_dir: Union[str, Path],
        calculator: str = "vasp",
        naming_func: Optional[callable] = None,
        **kwargs,
    ):
        """
        Initialize BatchRelaxation.

        Parameters
        ----------
        structures : list
            List of Atoms objects
        base_dir : str or Path
            Base directory for calculations
        calculator : str
            Calculator type
        naming_func : callable, optional
            Function to generate directory names from index
        **kwargs : dict
            Parameters passed to RelaxationWorkflow
        """
        self.structures = structures
        self.base_dir = Path(base_dir)
        self.calculator = calculator
        self.kwargs = kwargs

        if naming_func is None:
            self.naming_func = lambda i: f"config_{i:03d}"
        else:
            self.naming_func = naming_func

        self.workflows = []

    def setup_all(self) -> List[Dict[str, Path]]:
        """
        Set up all relaxation calculations.

        Returns
        -------
        list
            List of file dictionaries for each calculation
        """
        all_files = []

        for i, atoms in enumerate(self.structures):
            name = self.naming_func(i)
            work_dir = self.base_dir / name

            wf = RelaxationWorkflow(
                atoms,
                work_dir=work_dir,
                calculator=self.calculator,
                **self.kwargs,
            )
            files = wf.setup()
            all_files.append(files)
            self.workflows.append(wf)

        print(f"Set up {len(self.structures)} relaxation calculations in {self.base_dir}")
        return all_files

    def parse_all(self) -> List[Dict[str, Any]]:
        """
        Parse results from all completed calculations.

        Returns
        -------
        list
            List of results dictionaries
        """
        results = []

        for i, wf in enumerate(self.workflows):
            try:
                result = wf.parse_results()
                result["index"] = i
                result["name"] = self.naming_func(i)
                results.append(result)
            except Exception as e:
                results.append({
                    "index": i,
                    "name": self.naming_func(i),
                    "error": str(e),
                })

        # Sort by energy
        completed = [r for r in results if "energy" in r and r["energy"] is not None]
        if completed:
            completed.sort(key=lambda x: x["energy"])
            print(f"Lowest energy: {completed[0]['name']} ({completed[0]['energy']:.4f} eV)")

        return results

    def get_status(self) -> Dict[str, int]:
        """
        Get status of all calculations.

        Returns
        -------
        dict
            Count of completed, running, pending calculations
        """
        status = {"completed": 0, "incomplete": 0, "not_started": 0}

        for wf in self.workflows:
            outcar = wf.work_dir / "OUTCAR"
            if outcar.exists():
                if wf.is_complete():
                    status["completed"] += 1
                else:
                    status["incomplete"] += 1
            else:
                status["not_started"] += 1

        return status


def relax_structure(
    atoms: Atoms,
    work_dir: Union[str, Path],
    calculator: str = "vasp",
    **kwargs,
) -> Dict[str, Any]:
    """
    Convenience function to set up a relaxation.

    Parameters
    ----------
    atoms : Atoms
        Structure to relax
    work_dir : str or Path
        Working directory
    calculator : str
        Calculator type
    **kwargs : dict
        Additional parameters

    Returns
    -------
    dict
        Setup results or optimization results (for MACE)
    """
    wf = RelaxationWorkflow(atoms, work_dir=work_dir, calculator=calculator, **kwargs)
    wf.setup()
    return wf.run()
