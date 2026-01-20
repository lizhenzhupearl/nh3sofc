"""High-throughput screening workflow.

Automates screening of multiple materials, compositions,
and configurations for catalyst discovery.
"""

from pathlib import Path
from typing import Optional, List, Dict, Any, Union, Callable
import itertools
import json
import numpy as np
from ase import Atoms

from ..structure.surface import SurfaceBuilder
from ..structure.defects import DefectBuilder
from ..structure.adsorbates import AdsorbatePlacer
from .relaxation import RelaxationWorkflow, BatchRelaxation
from ..database.naming import NamingConvention, DirectoryStructure


class ScreeningWorkflow:
    """
    High-throughput screening workflow.

    Generates and sets up calculations for systematic parameter scans.

    Examples
    --------
    >>> wf = ScreeningWorkflow(
    ...     base_structure=bulk,
    ...     parameter_space={
    ...         "vacancy_concentration": [0.0, 0.05, 0.10],
    ...         "miller": [(0,0,1), (1,1,0)],
    ...     }
    ... )
    >>> wf.generate_structures()
    >>> wf.setup_calculations()
    """

    def __init__(
        self,
        base_structure: Atoms,
        parameter_space: Dict[str, List[Any]],
        work_dir: Union[str, Path] = "./screening",
        calculator: str = "vasp",
        material_name: str = "material",
        **calc_params,
    ):
        """
        Initialize ScreeningWorkflow.

        Parameters
        ----------
        base_structure : Atoms
            Base structure (bulk or slab)
        parameter_space : dict
            Parameters to scan: {param_name: [values]}
        work_dir : str or Path
            Base working directory
        calculator : str
            Calculator type ("vasp" or "mace")
        material_name : str
            Material identifier
        **calc_params : dict
            Calculation parameters
        """
        self.base_structure = base_structure.copy()
        self.parameter_space = parameter_space
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        self.calculator = calculator
        self.material_name = material_name
        self.calc_params = calc_params

        # Storage
        self.configurations = []
        self.config_params = []
        self.workflows = []
        self.results = []

    def _get_parameter_combinations(self) -> List[Dict[str, Any]]:
        """Generate all parameter combinations."""
        keys = list(self.parameter_space.keys())
        values = list(self.parameter_space.values())

        combinations = []
        for combo in itertools.product(*values):
            combinations.append(dict(zip(keys, combo)))

        return combinations

    def _generate_config_name(self, params: Dict[str, Any]) -> str:
        """Generate configuration name from parameters."""
        parts = []
        for key, value in sorted(params.items()):
            if isinstance(value, tuple):
                val_str = "".join(str(v) for v in value)
            elif isinstance(value, float):
                val_str = f"{value:.2f}".replace(".", "p")
            else:
                val_str = str(value)
            parts.append(f"{key[:3]}_{val_str}")
        return "_".join(parts)

    def generate_structures(
        self,
        structure_generator: Optional[Callable] = None,
    ) -> List[Atoms]:
        """
        Generate structures for all parameter combinations.

        Parameters
        ----------
        structure_generator : callable, optional
            Custom function: (base_structure, params) -> Atoms

        Returns
        -------
        list
            Generated structures
        """
        combinations = self._get_parameter_combinations()

        self.configurations = []
        self.config_params = []

        for params in combinations:
            if structure_generator:
                atoms = structure_generator(self.base_structure, params)
            else:
                atoms = self._default_structure_generator(params)

            if atoms is not None:
                self.configurations.append(atoms)
                self.config_params.append(params)

        print(f"Generated {len(self.configurations)} configurations")
        return self.configurations

    def _default_structure_generator(
        self,
        params: Dict[str, Any],
    ) -> Optional[Atoms]:
        """Default structure generator for common parameter types."""
        atoms = self.base_structure.copy()

        # Handle miller indices (create surface)
        if "miller" in params:
            miller = params["miller"]
            builder = SurfaceBuilder(atoms)
            layers = params.get("layers", 4)
            vacuum = params.get("vacuum", 15.0)
            atoms = builder.create_surface(miller, layers=layers, vacuum=vacuum)

        # Handle defects (oxynitride creation)
        if "vacancy_concentration" in params or "nitrogen_fraction" in params:
            vac = params.get("vacancy_concentration", 0.0)
            n_frac = params.get("nitrogen_fraction", 0.0)

            if vac > 0 or n_frac > 0:
                defect_builder = DefectBuilder(atoms)
                atoms = defect_builder.create_oxynitride(
                    nitrogen_fraction=n_frac,
                    vacancy_concentration=vac,
                )

        # Handle adsorbate
        if "adsorbate" in params:
            adsorbate = params["adsorbate"]
            placer = AdsorbatePlacer(atoms)
            height = params.get("ads_height", 2.0)
            atoms = placer.add_random(adsorbate, height=height, n_configs=1)[0]

        return atoms

    def setup_calculations(self) -> List[Dict[str, Path]]:
        """
        Set up calculations for all configurations.

        Returns
        -------
        list
            List of file dictionaries
        """
        if not self.configurations:
            self.generate_structures()

        all_files = []
        self.workflows = []

        for i, (atoms, params) in enumerate(zip(self.configurations, self.config_params)):
            config_name = self._generate_config_name(params)
            config_dir = self.work_dir / config_name

            wf = RelaxationWorkflow(
                atoms,
                work_dir=config_dir,
                calculator=self.calculator,
                **self.calc_params,
            )

            files = wf.setup()
            files["config_name"] = config_name
            files["params"] = params

            all_files.append(files)
            self.workflows.append(wf)

        print(f"Set up {len(self.workflows)} calculations in {self.work_dir}")
        return all_files

    def parse_results(self) -> List[Dict[str, Any]]:
        """
        Parse results from completed calculations.

        Returns
        -------
        list
            Results for each configuration
        """
        self.results = []

        for wf, params in zip(self.workflows, self.config_params):
            try:
                result = wf.parse_results()
                result["params"] = params
                result["config_name"] = self._generate_config_name(params)
                self.results.append(result)
            except Exception as e:
                self.results.append({
                    "params": params,
                    "config_name": self._generate_config_name(params),
                    "error": str(e),
                })

        return self.results

    def get_best_result(
        self,
        metric: str = "energy",
        minimize: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """
        Get best result by metric.

        Parameters
        ----------
        metric : str
            Metric to optimize
        minimize : bool
            Whether to minimize (True) or maximize (False)

        Returns
        -------
        dict or None
            Best result
        """
        if not self.results:
            self.parse_results()

        valid_results = [
            r for r in self.results
            if metric in r and r[metric] is not None
        ]

        if not valid_results:
            return None

        if minimize:
            return min(valid_results, key=lambda x: x[metric])
        else:
            return max(valid_results, key=lambda x: x[metric])

    def get_parameter_trend(
        self,
        param_name: str,
        metric: str = "energy",
    ) -> Dict[Any, List[float]]:
        """
        Get metric values grouped by parameter.

        Parameters
        ----------
        param_name : str
            Parameter name
        metric : str
            Metric to analyze

        Returns
        -------
        dict
            {param_value: [metric_values]}
        """
        if not self.results:
            self.parse_results()

        trends = {}
        for result in self.results:
            if metric not in result or result[metric] is None:
                continue

            params = result.get("params", {})
            if param_name not in params:
                continue

            param_val = params[param_name]
            if param_val not in trends:
                trends[param_val] = []
            trends[param_val].append(result[metric])

        return trends

    def save_results(self, filename: Optional[str] = None) -> None:
        """
        Save results to JSON file.

        Parameters
        ----------
        filename : str, optional
            Output file path
        """
        if not self.results:
            self.parse_results()

        if filename is None:
            filename = self.work_dir / "screening_results.json"

        # Convert to serializable format
        data = []
        for result in self.results:
            entry = {}
            for key, value in result.items():
                if isinstance(value, np.ndarray):
                    entry[key] = value.tolist()
                elif isinstance(value, Atoms):
                    entry[key] = value.get_chemical_formula()
                elif isinstance(value, Path):
                    entry[key] = str(value)
                else:
                    try:
                        json.dumps(value)
                        entry[key] = value
                    except (TypeError, ValueError):
                        entry[key] = str(value)
            data.append(entry)

        with open(filename, "w") as f:
            json.dump(data, f, indent=2)

        print(f"Saved results to {filename}")

    def get_status(self) -> Dict[str, int]:
        """
        Get status of all calculations.

        Returns
        -------
        dict
            Status counts
        """
        status = {"completed": 0, "incomplete": 0, "not_started": 0}

        for wf in self.workflows:
            if wf.is_complete():
                status["completed"] += 1
            elif (wf.work_dir / "OUTCAR").exists():
                status["incomplete"] += 1
            else:
                status["not_started"] += 1

        return status

    def print_summary(self) -> None:
        """Print screening summary."""
        print("\n" + "=" * 60)
        print("Screening Workflow Summary")
        print("=" * 60)

        print(f"\nMaterial: {self.material_name}")
        print(f"Calculator: {self.calculator}")
        print(f"Work directory: {self.work_dir}")

        print("\nParameter space:")
        for param, values in self.parameter_space.items():
            print(f"  {param}: {values}")

        n_combos = len(self._get_parameter_combinations())
        print(f"\nTotal combinations: {n_combos}")
        print(f"Configurations generated: {len(self.configurations)}")

        if self.workflows:
            status = self.get_status()
            print(f"\nCalculation status:")
            for s, count in status.items():
                print(f"  {s}: {count}")

        if self.results:
            valid = [r for r in self.results if "energy" in r]
            if valid:
                energies = [r["energy"] for r in valid if r["energy"] is not None]
                if energies:
                    print(f"\nEnergy statistics:")
                    print(f"  Min: {min(energies):.4f} eV")
                    print(f"  Max: {max(energies):.4f} eV")
                    print(f"  Mean: {np.mean(energies):.4f} eV")

        print("=" * 60)


class CompositionScreening(ScreeningWorkflow):
    """
    Specialized screening for composition variations.

    Screens oxynitride compositions with varying N/O ratios
    and vacancy concentrations.
    """

    def __init__(
        self,
        base_structure: Atoms,
        nitrogen_fractions: List[float] = [0.0, 0.33, 0.5, 0.67, 1.0],
        vacancy_concentrations: List[float] = [0.0, 0.05, 0.1],
        **kwargs,
    ):
        """
        Initialize CompositionScreening.

        Parameters
        ----------
        base_structure : Atoms
            Base oxide structure
        nitrogen_fractions : list
            N/(N+O) fractions to screen
        vacancy_concentrations : list
            Vacancy concentrations to screen
        **kwargs : dict
            Additional parameters
        """
        parameter_space = {
            "nitrogen_fraction": nitrogen_fractions,
            "vacancy_concentration": vacancy_concentrations,
        }

        super().__init__(
            base_structure,
            parameter_space=parameter_space,
            **kwargs,
        )


class AdsorbateScreening(ScreeningWorkflow):
    """
    Specialized screening for adsorbate configurations.

    Screens different adsorption sites and orientations.
    """

    def __init__(
        self,
        slab: Atoms,
        adsorbate: str = "NH3",
        n_configs: int = 10,
        site_types: List[str] = ["ontop", "bridge", "hollow"],
        **kwargs,
    ):
        """
        Initialize AdsorbateScreening.

        Parameters
        ----------
        slab : Atoms
            Surface slab
        adsorbate : str
            Adsorbate molecule
        n_configs : int
            Configurations per site type
        site_types : list
            Site types to screen
        **kwargs : dict
            Additional parameters
        """
        parameter_space = {
            "site_type": site_types,
            "config_id": list(range(n_configs)),
        }

        super().__init__(
            slab,
            parameter_space=parameter_space,
            **kwargs,
        )

        self.adsorbate = adsorbate

    def _default_structure_generator(
        self,
        params: Dict[str, Any],
    ) -> Optional[Atoms]:
        """Generate adsorbate configuration."""
        site_type = params.get("site_type", "ontop")
        config_id = params.get("config_id", 0)

        placer = AdsorbatePlacer(self.base_structure)

        # Use different seed for each config
        np.random.seed(config_id)

        try:
            if site_type == "ontop":
                configs = placer.add_on_site(self.adsorbate, site_type="ontop")
            elif site_type == "bridge":
                configs = placer.add_on_site(self.adsorbate, site_type="bridge")
            elif site_type == "hollow":
                configs = placer.add_on_site(self.adsorbate, site_type="hollow")
            else:
                configs = placer.add_random(self.adsorbate, n_configs=1)

            if configs and len(configs) > config_id:
                return configs[config_id]
            elif configs:
                return configs[0]
        except Exception:
            # Fall back to random placement
            configs = placer.add_random(self.adsorbate, n_configs=1)
            return configs[0] if configs else None

        return None


def run_screening(
    base_structure: Atoms,
    parameter_space: Dict[str, List[Any]],
    work_dir: str = "./screening",
    calculator: str = "vasp",
    **kwargs,
) -> ScreeningWorkflow:
    """
    Convenience function to run a screening workflow.

    Parameters
    ----------
    base_structure : Atoms
        Base structure
    parameter_space : dict
        Parameters to scan
    work_dir : str
        Working directory
    calculator : str
        Calculator type
    **kwargs : dict
        Additional parameters

    Returns
    -------
    ScreeningWorkflow
        Configured workflow
    """
    wf = ScreeningWorkflow(
        base_structure,
        parameter_space=parameter_space,
        work_dir=work_dir,
        calculator=calculator,
        **kwargs,
    )

    wf.generate_structures()
    wf.setup_calculations()

    return wf
