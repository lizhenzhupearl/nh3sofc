"""NEB transition state search workflow.

Provides automated NEB calculations for finding transition states
and reaction barriers.
"""

from pathlib import Path
from typing import Optional, List, Union, Dict, Any, Tuple
import numpy as np
from ase import Atoms
from ase.io import write as ase_write, read as ase_read
try:
    from ase.mep import NEB
except ImportError:
    from ase.neb import NEB

from ..core.base import BaseWorkflow
from ..calculators.vasp.inputs import VASPInputGenerator
from ..calculators.vasp.outputs import VASPOutputParser, NEBOutputParser
from ..jobs.pbs import VASPJobScript, create_vasp_job


class NEBWorkflow(BaseWorkflow):
    """
    Workflow for NEB transition state calculations.

    Supports VASP NEB (via job submission) or MACE NEB (direct).

    Examples
    --------
    >>> # Set up VASP NEB
    >>> wf = NEBWorkflow(
    ...     initial=initial_atoms,
    ...     final=final_atoms,
    ...     work_dir="./neb",
    ...     n_images=5,
    ...     climbing=True,
    ... )
    >>> wf.setup()
    >>> # Submit job manually, then:
    >>> results = wf.parse_results()
    >>> barrier = results["forward_barrier"]
    """

    def __init__(
        self,
        initial: Atoms,
        final: Atoms,
        work_dir: Union[str, Path] = "./",
        calculator: str = "vasp",
        n_images: int = 5,
        climbing: bool = True,
        spring_constant: float = 5.0,
        fmax: float = 0.05,
        # VASP parameters
        encut: float = 520,
        kspacing: float = 0.03,
        hubbard_u: Optional[Dict[str, float]] = None,
        vdw: Optional[str] = None,
        is_surface: bool = True,
        # PBS parameters
        nodes: int = 1,
        ppn: int = 24,
        walltime: str = "48:00:00",
        queue: Optional[str] = None,
        # MACE parameters
        mace_model: Optional[str] = None,
        **kwargs,
    ):
        """
        Initialize NEBWorkflow.

        Parameters
        ----------
        initial : Atoms
            Initial (reactant) structure
        final : Atoms
            Final (product) structure
        work_dir : str or Path
            Working directory
        calculator : str
            Calculator type: "vasp" or "mace"
        n_images : int
            Number of intermediate images
        climbing : bool
            Use climbing image NEB
        spring_constant : float
            Spring constant for NEB (eV/A^2)
        fmax : float
            Force convergence criterion
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
        mace_model : str, optional
            Path to MACE model
        **kwargs : dict
            Additional parameters
        """
        # Use initial structure as reference
        super().__init__(initial, work_dir, calculator, **kwargs)

        self.initial = initial.copy()
        self.final = final.copy()
        self.n_images = n_images
        self.climbing = climbing
        self.spring_constant = spring_constant
        self.fmax = fmax

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

        # Storage
        self.images = []
        self.interpolated = False

    def interpolate_images(
        self,
        method: str = "idpp",
    ) -> List[Atoms]:
        """
        Generate interpolated images between initial and final.

        Parameters
        ----------
        method : str
            Interpolation method: "linear" or "idpp"

        Returns
        -------
        list
            List of images including initial and final
        """
        # Create image list
        self.images = [self.initial.copy()]
        for _ in range(self.n_images):
            self.images.append(self.initial.copy())
        self.images.append(self.final.copy())

        # Interpolate
        if method == "idpp":
            # IDPP interpolation for better initial path
            try:
                from ase.neb import idpp_interpolate
                idpp_interpolate(self.images)
            except ImportError:
                # Fallback to ASE's interpolate
                neb = NEB(self.images)
                neb.interpolate(method='idpp')
        else:
            # Linear interpolation
            neb = NEB(self.images)
            neb.interpolate()

        self.interpolated = True
        return self.images

    def setup(self) -> Dict[str, Any]:
        """
        Set up the NEB calculation.

        For VASP: generates input files for each image and PBS script.
        For MACE: generates initial images only.

        Returns
        -------
        dict
            Dictionary of generated file paths and info
        """
        if not self.interpolated:
            self.interpolate_images()

        files = {"images": []}

        if self.calculator == "vasp":
            return self._setup_vasp()
        else:
            # For MACE, just save the interpolated images
            for i, image in enumerate(self.images):
                img_path = self.work_dir / f"image_{i:02d}.traj"
                ase_write(str(img_path), image)
                files["images"].append(img_path)

            return files

    def _setup_vasp(self) -> Dict[str, Any]:
        """Set up VASP NEB calculation."""
        files = {"images": [], "dirs": []}

        # Save initial and final
        ase_write(str(self.work_dir / "initial.traj"), self.initial)
        ase_write(str(self.work_dir / "final.traj"), self.final)
        files["initial"] = self.work_dir / "initial.traj"
        files["final"] = self.work_dir / "final.traj"

        # Generate VASP NEB inputs
        gen = VASPInputGenerator(
            self.initial,
            calc_type="neb",
            work_dir=self.work_dir,
            IMAGES=self.n_images,
            LCLIMB=self.climbing,
            SPRING=-self.spring_constant,
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

        # Create image directories (00, 01, ..., N+1)
        for i, image in enumerate(self.images):
            img_dir = self.work_dir / f"{i:02d}"
            img_dir.mkdir(parents=True, exist_ok=True)

            # Write POSCAR
            poscar_path = img_dir / "POSCAR"
            ase_write(str(poscar_path), image, format="vasp")

            files["dirs"].append(img_dir)
            files["images"].append(poscar_path)

        # Generate PBS script for NEB
        pbs_script = VASPJobScript(
            work_dir=self.work_dir,
            nodes=self.nodes,
            ppn=self.ppn,
            walltime=self.walltime,
            queue=self.queue,
        )

        # Modify for NEB (need MPI across images)
        total_procs = self.nodes * self.ppn
        procs_per_image = total_procs // (self.n_images + 2)

        pbs_path = pbs_script.write()
        files["PBS"] = pbs_path

        print(f"Set up NEB with {self.n_images} images")
        print(f"Submit with: qsub {pbs_path}")

        return files

    def run(self) -> Dict[str, Any]:
        """
        Run the NEB calculation.

        For VASP: returns instructions (manual submission required).
        For MACE: runs NEB directly.

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
        """Run NEB with MACE calculator."""
        try:
            from mace.calculators import MACECalculator
        except ImportError:
            raise ImportError("MACE not installed. Install with: pip install mace-torch")

        from ase.optimize import BFGS

        # Load MACE model
        if self.mace_model is None:
            try:
                from mace.calculators import mace_mp
                calc = mace_mp(model="medium", device="cpu")
            except Exception:
                raise ValueError("No MACE model specified")
        else:
            calc = MACECalculator(model_path=self.mace_model, device="cpu")

        # Ensure images are interpolated
        if not self.interpolated:
            self.interpolate_images()

        # Attach calculators
        for image in self.images[1:-1]:  # Skip endpoints
            image.calc = calc

        # Also calculate endpoints
        self.initial.calc = calc
        self.final.calc = calc

        # Create NEB
        neb = NEB(
            self.images,
            climb=self.climbing,
            k=self.spring_constant,
        )

        # Optimize
        optimizer = BFGS(neb, trajectory=str(self.work_dir / "neb.traj"))
        optimizer.run(fmax=self.fmax, steps=500)

        # Extract results
        energies = [img.get_potential_energy() for img in self.images]
        e_initial = energies[0]
        e_final = energies[-1]
        e_max = max(energies)
        i_max = energies.index(e_max)

        self.results = {
            "converged": optimizer.converged(),
            "energies": energies,
            "forward_barrier": e_max - e_initial,
            "reverse_barrier": e_max - e_final,
            "reaction_energy": e_final - e_initial,
            "ts_index": i_max,
            "ts_structure": self.images[i_max].copy(),
            "n_steps": optimizer.nsteps,
            "images": [img.copy() for img in self.images],
        }

        # Save results
        for i, image in enumerate(self.images):
            ase_write(str(self.work_dir / f"final_{i:02d}.traj"), image)

        return self.results

    def parse_results(self) -> Dict[str, Any]:
        """
        Parse results from completed VASP NEB calculation.

        Returns
        -------
        dict
            Results dictionary including barriers and TS structure
        """
        parser = NEBOutputParser(self.work_dir)

        # Get energies and barrier
        energies = parser.get_energies()
        barrier_info = parser.get_barrier()

        self.results = {
            "converged": parser.is_converged(),
            "energies": energies,
            "forward_barrier": barrier_info["forward"],
            "reverse_barrier": barrier_info["reverse"],
            "reaction_energy": barrier_info["reaction_energy"],
            "ts_index": barrier_info["ts_index"],
            "ts_structure": parser.get_ts_structure(),
            "images": parser.get_all_images(),
        }

        return self.results

    def get_energy_path(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get energy along reaction coordinate.

        Returns
        -------
        tuple
            (reaction_coordinate, energies) arrays
        """
        if not self.results:
            self.parse_results()

        energies = np.array(self.results["energies"])
        n_points = len(energies)

        # Normalized reaction coordinate
        coord = np.linspace(0, 1, n_points)

        return coord, energies

    def plot_energy_profile(self, filename: Optional[str] = None):
        """
        Plot NEB energy profile.

        Parameters
        ----------
        filename : str, optional
            Output file path. If None, displays plot.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not installed. Cannot plot.")
            return

        coord, energies = self.get_energy_path()

        # Reference to initial state
        energies = energies - energies[0]

        fig, ax = plt.subplots(figsize=(8, 6))

        ax.plot(coord, energies, 'o-', markersize=10, linewidth=2)
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

        ax.set_xlabel("Reaction Coordinate", fontsize=12)
        ax.set_ylabel("Energy (eV)", fontsize=12)
        ax.set_title("NEB Energy Profile", fontsize=14)

        # Annotate barrier
        e_max = max(energies)
        i_max = list(energies).index(e_max)
        ax.annotate(
            f"Barrier: {e_max:.3f} eV",
            xy=(coord[i_max], e_max),
            xytext=(coord[i_max] + 0.1, e_max + 0.1),
            fontsize=10,
            arrowprops=dict(arrowstyle="->", color="red"),
        )

        plt.tight_layout()

        if filename:
            plt.savefig(filename, dpi=150)
            print(f"Saved plot to {filename}")
        else:
            plt.show()


class AutoNEB(NEBWorkflow):
    """
    Automated NEB with adaptive image insertion.

    Starts with few images and adds more near the transition state.
    """

    def __init__(
        self,
        initial: Atoms,
        final: Atoms,
        work_dir: Union[str, Path] = "./",
        n_images_start: int = 3,
        n_images_max: int = 9,
        **kwargs,
    ):
        """
        Initialize AutoNEB.

        Parameters
        ----------
        initial : Atoms
            Initial structure
        final : Atoms
            Final structure
        work_dir : str or Path
            Working directory
        n_images_start : int
            Starting number of images
        n_images_max : int
            Maximum number of images
        **kwargs : dict
            Additional parameters
        """
        super().__init__(
            initial, final, work_dir,
            n_images=n_images_start,
            **kwargs
        )

        self.n_images_start = n_images_start
        self.n_images_max = n_images_max
        self.iteration = 0

    def add_images_near_ts(self) -> int:
        """
        Add images near the transition state.

        Returns
        -------
        int
            Number of images added
        """
        if not self.results:
            raise ValueError("Run initial NEB first")

        if len(self.images) >= self.n_images_max + 2:
            print(f"Already at max images ({self.n_images_max})")
            return 0

        # Find TS location
        energies = self.results["energies"]
        ts_idx = self.results["ts_index"]

        # Add image before and after TS
        new_images = []
        for i, image in enumerate(self.images):
            new_images.append(image)

            # Add interpolated image near TS
            if i == ts_idx - 1 or i == ts_idx:
                if len(new_images) < self.n_images_max + 2:
                    # Interpolate between current and next
                    next_img = self.images[i + 1]
                    interp = self.images[i].copy()
                    interp.positions = 0.5 * (image.positions + next_img.positions)
                    new_images.append(interp)

        added = len(new_images) - len(self.images)
        self.images = new_images
        self.n_images = len(self.images) - 2

        print(f"Added {added} images near TS. Total: {self.n_images}")
        return added


def setup_neb_calculation(
    initial: Atoms,
    final: Atoms,
    work_dir: Union[str, Path],
    n_images: int = 5,
    climbing: bool = True,
    calculator: str = "vasp",
    **kwargs,
) -> NEBWorkflow:
    """
    Convenience function to set up an NEB calculation.

    Parameters
    ----------
    initial : Atoms
        Initial structure
    final : Atoms
        Final structure
    work_dir : str or Path
        Working directory
    n_images : int
        Number of intermediate images
    climbing : bool
        Use climbing image NEB
    calculator : str
        Calculator type
    **kwargs : dict
        Additional parameters

    Returns
    -------
    NEBWorkflow
        Configured workflow
    """
    wf = NEBWorkflow(
        initial=initial,
        final=final,
        work_dir=work_dir,
        n_images=n_images,
        climbing=climbing,
        calculator=calculator,
        **kwargs,
    )
    wf.setup()
    return wf


def create_neb_from_pathway(
    pathway: Dict[str, List[Atoms]],
    step_pairs: List[Tuple[str, str]],
    base_dir: Union[str, Path],
    **kwargs,
) -> Dict[str, NEBWorkflow]:
    """
    Create NEB calculations from decomposition pathway.

    Parameters
    ----------
    pathway : dict
        Dictionary of step configurations
    step_pairs : list
        List of (initial_step, final_step) tuples
    base_dir : str or Path
        Base directory for NEB calculations
    **kwargs : dict
        Parameters for NEBWorkflow

    Returns
    -------
    dict
        Dictionary of NEBWorkflow objects
    """
    base_dir = Path(base_dir)
    workflows = {}

    for initial_step, final_step in step_pairs:
        if initial_step not in pathway or final_step not in pathway:
            print(f"Skipping {initial_step} -> {final_step}: missing structures")
            continue

        # Use lowest energy configs
        initial_configs = pathway[initial_step]
        final_configs = pathway[final_step]

        # For now, use first config of each
        # In practice, should use optimized structures
        initial = initial_configs[0]
        final = final_configs[0]

        neb_name = f"{initial_step}_to_{final_step}"
        neb_dir = base_dir / neb_name

        wf = NEBWorkflow(
            initial=initial,
            final=final,
            work_dir=neb_dir,
            **kwargs,
        )

        workflows[neb_name] = wf
        print(f"Created NEB: {neb_name}")

    return workflows
