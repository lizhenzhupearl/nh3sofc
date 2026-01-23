"""Defect creation for oxynitride materials.

Creates oxynitrides by:
1. O → N substitution (nitrogen_fraction parameter)
2. Vacancy creation (vacancy_concentration parameter x)

For materials like LaVO3 → LaVON₂₋ₓ where x controls vacancy concentration.
"""

from pathlib import Path
from typing import Optional, List, Union, Tuple, Dict, Any
import numpy as np
from ase import Atoms
from ase.io import write

from ..core.base import BaseStructure
from .surface import SlabStructure


class DefectBuilder:
    """
    Builder class for creating defects in structures.

    Supports O→N substitution and vacancy creation for oxynitride materials.

    Examples
    --------
    >>> slab = SlabStructure.from_file("LaVO3_001.cif")
    >>> builder = DefectBuilder(slab)
    >>> oxynitride = builder.create_oxynitride(
    ...     nitrogen_fraction=0.67,      # 2/3 of O → N
    ...     vacancy_concentration=0.1,   # 10% vacancies
    ... )
    """

    def __init__(self, structure: Union[SlabStructure, Atoms]):
        """
        Initialize DefectBuilder.

        Parameters
        ----------
        structure : SlabStructure or Atoms
            Structure to create defects in
        """
        if isinstance(structure, SlabStructure):
            self.atoms = structure.atoms.copy()
            self._structure_obj = structure
        elif isinstance(structure, BaseStructure):
            self.atoms = structure.atoms.copy()
            self._structure_obj = structure
        else:
            self.atoms = structure.copy()
            self._structure_obj = None

    def get_element_indices(self, element: str) -> List[int]:
        """Get indices of atoms of a specific element."""
        symbols = self.atoms.get_chemical_symbols()
        return [i for i, s in enumerate(symbols) if s == element]

    def substitute(
        self,
        from_element: str,
        to_element: str,
        fraction: float = 1.0,
        indices: Optional[List[int]] = None,
        random_seed: Optional[int] = None,
    ) -> Atoms:
        """
        Substitute one element with another.

        Parameters
        ----------
        from_element : str
            Element to replace (e.g., "O")
        to_element : str
            Element to substitute (e.g., "N")
        fraction : float
            Fraction of atoms to substitute (0.0 to 1.0)
        indices : list, optional
            Specific indices to substitute. If None, random selection.
        random_seed : int, optional
            Random seed for reproducibility

        Returns
        -------
        Atoms
            Structure with substitutions
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        new_atoms = self.atoms.copy()
        symbols = list(new_atoms.get_chemical_symbols())

        # Get indices of target element
        if indices is None:
            target_indices = self.get_element_indices(from_element)
            n_to_substitute = int(len(target_indices) * fraction)
            indices = np.random.choice(
                target_indices, size=n_to_substitute, replace=False
            ).tolist()

        # Perform substitution
        for idx in indices:
            symbols[idx] = to_element

        new_atoms.set_chemical_symbols(symbols)
        return new_atoms

    def create_vacancy(
        self,
        element: str,
        concentration: float,
        indices: Optional[List[int]] = None,
        random_seed: Optional[int] = None,
    ) -> Atoms:
        """
        Create vacancies by removing atoms.

        Parameters
        ----------
        element : str
            Element to remove (e.g., "O" for oxygen vacancies)
        concentration : float
            Fraction of atoms to remove (0.0 to 1.0)
        indices : list, optional
            Specific indices to remove. If None, random selection.
        random_seed : int, optional
            Random seed

        Returns
        -------
        Atoms
            Structure with vacancies
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        # Get indices of target element
        if indices is None:
            target_indices = self.get_element_indices(element)
            n_to_remove = int(len(target_indices) * concentration)
            if n_to_remove == 0 and concentration > 0:
                n_to_remove = 1  # Remove at least one if concentration > 0
            indices = np.random.choice(
                target_indices, size=n_to_remove, replace=False
            ).tolist()

        # Create new atoms without the vacancies
        keep_indices = [i for i in range(len(self.atoms)) if i not in indices]
        new_atoms = self.atoms[keep_indices]

        return new_atoms

    def _get_surface_indices(
        self,
        indices: List[int],
        z_threshold: float = 0.3,
    ) -> Tuple[List[int], List[int]]:
        """
        Split indices into surface and bulk based on z-position.

        Parameters
        ----------
        indices : list
            Atom indices to split
        z_threshold : float
            Fraction of slab height considered as surface region (0.3 = top 30%)

        Returns
        -------
        tuple
            (surface_indices, bulk_indices)
        """
        z_positions = self.atoms.positions[:, 2]
        z_min, z_max = z_positions.min(), z_positions.max()
        z_cutoff = z_max - z_threshold * (z_max - z_min)

        surface_indices = [i for i in indices if z_positions[i] >= z_cutoff]
        bulk_indices = [i for i in indices if z_positions[i] < z_cutoff]

        return surface_indices, bulk_indices

    def _calculate_selection_weights(
        self,
        candidates: List[int],
        placement: str = "random",
        preference: float = 0.7,
        n_positions: Optional[np.ndarray] = None,
        z_threshold: float = 0.3,
    ) -> np.ndarray:
        """
        Calculate selection weights for each candidate based on placement strategy.

        Parameters
        ----------
        candidates : list
            Candidate atom indices
        placement : str
            Placement strategy: "random", "surface", or "near_N"
        preference : float
            Preference strength (0.5 = random, 1.0 = strongly prefer target region)
            For "surface": fraction of surface atoms to be selected type (e.g., 0.8 = 80% N at surface)
            For "near_N": how strongly vacancies prefer to be near N
        n_positions : ndarray, optional
            Positions of N atoms (required for "near_N" placement)
        z_threshold : float
            Fraction of slab height considered as surface region (0.3 = top 30%)

        Returns
        -------
        ndarray
            Normalized weights for each candidate
        """
        n_candidates = len(candidates)
        if n_candidates == 0:
            return np.array([])

        if placement == "random" or preference <= 0.5:
            # Uniform weights
            return np.ones(n_candidates) / n_candidates

        elif placement == "surface":
            # Weight by z-position: higher z = higher weight
            z_positions = self.atoms.positions[:, 2]
            z_min = z_positions.min()
            z_max = z_positions.max()
            z_range = z_max - z_min

            if z_range < 0.1:
                return np.ones(n_candidates) / n_candidates

            # Normalize z to [0, 1] for candidates
            z_values = np.array([z_positions[i] for i in candidates])
            z_normalized = (z_values - z_min) / z_range

            # Apply preference: weights increase exponentially with z
            # preference controls steepness: 0.5 = flat, 1.0 = strongly surface-biased
            steepness = (preference - 0.5) * 10  # 0 to 5
            weights = np.exp(steepness * z_normalized)

            return weights / weights.sum()

        elif placement == "near_N" and n_positions is not None and len(n_positions) > 0:
            # Weight by inverse distance to nearest N
            positions = self.atoms.positions
            distances = []
            for idx in candidates:
                pos = positions[idx]
                min_dist = np.min(np.linalg.norm(n_positions - pos, axis=1))
                distances.append(max(min_dist, 0.1))  # Avoid division by zero

            distances = np.array(distances)

            # Apply preference: closer to N = higher weight
            steepness = (preference - 0.5) * 10
            weights = np.exp(-steepness * distances / distances.max())

            return weights / weights.sum()

        else:
            return np.ones(n_candidates) / n_candidates

    def _select_with_placement(
        self,
        candidates: List[int],
        n_select: int,
        placement: str = "random",
        preference: float = 0.7,
        n_positions: Optional[np.ndarray] = None,
        z_threshold: float = 0.3,
    ) -> List[int]:
        """
        Select atoms based on placement strategy with probability weighting.

        Uses weighted random sampling where atoms in preferred regions
        have higher probability of selection, but selection is still stochastic.

        Parameters
        ----------
        candidates : list
            Candidate atom indices
        n_select : int
            Number of atoms to select
        placement : str
            Placement strategy: "random", "surface", or "near_N"
        preference : float
            Preference strength (0.5 = random, 1.0 = strongly prefer target region)
        n_positions : ndarray, optional
            Positions of N atoms (required for "near_N" placement)
        z_threshold : float
            Fraction of slab height considered as surface region (0.3 = top 30%)

        Returns
        -------
        list
            Selected atom indices
        """
        if n_select <= 0 or not candidates:
            return []

        n_select = min(n_select, len(candidates))

        # Calculate selection weights
        weights = self._calculate_selection_weights(
            candidates, placement, preference, n_positions, z_threshold
        )

        # Weighted random selection without replacement
        selected = np.random.choice(
            candidates,
            size=n_select,
            replace=False,
            p=weights
        ).tolist()

        return selected

    def create_oxynitride(
        self,
        nitrogen_fraction: float = 0.67,
        vacancy_concentration: float = 0.0,
        vacancy_element: str = "N",
        placement: str = "random",
        surface_n_preference: float = 0.7,
        vacancy_preference: float = 0.7,
        z_threshold: float = 0.3,
        random_seed: Optional[int] = None,
    ) -> Atoms:
        """
        Create oxynitride by O→N substitution and optional vacancies.

        For LaVO3 → LaVON₂₋ₓ:
        - nitrogen_fraction=0.67 replaces 2/3 of O with N
        - vacancy_concentration=0.1 creates 10% vacancies in N sites

        Parameters
        ----------
        nitrogen_fraction : float
            Fraction of O atoms to replace with N (0.0 to 1.0)
        vacancy_concentration : float
            Fraction of anion sites to leave vacant (parameter x)
        vacancy_element : str
            Element to create vacancies in ("N" or "O")
        placement : str
            Defect placement strategy:
            - "random": Random placement (default)
            - "surface": Preferentially place N and vacancies near surface
            - "near_N": Place vacancies near existing N atoms
        surface_n_preference : float
            For "surface" placement: preference for N at surface (0.5-1.0)
            0.5 = random distribution, 0.7 = ~70% of surface anions are N,
            1.0 = strongly prefer N at surface. Default: 0.7
        vacancy_preference : float
            Preference strength for vacancy placement (0.5-1.0)
            For "surface": preference for vacancies at surface
            For "near_N": preference for vacancies near N atoms
            0.5 = random, 1.0 = strongly prefer target region. Default: 0.7
        z_threshold : float
            Fraction of slab height considered as surface region (0.3 = top 30%).
            Must match the value used in analyze_defect_distribution() for
            consistent analysis. Default: 0.3
        random_seed : int, optional
            Random seed for reproducibility

        Returns
        -------
        Atoms
            Oxynitride structure

        Examples
        --------
        >>> # Create LaVON₂ (no vacancies, random N distribution)
        >>> oxynitride = builder.create_oxynitride(nitrogen_fraction=0.67)
        >>>
        >>> # Create LaVON₁.₉ with N preferring surface (70% surface anions are N)
        >>> oxynitride = builder.create_oxynitride(
        ...     nitrogen_fraction=0.67,
        ...     vacancy_concentration=0.1,
        ...     placement="surface",
        ...     surface_n_preference=0.7,  # 70% of surface anions will be N
        ...     vacancy_preference=0.6,    # Vacancies slightly prefer surface
        ... )
        >>>
        >>> # Near-N vacancies with moderate preference
        >>> oxynitride = builder.create_oxynitride(
        ...     nitrogen_fraction=0.67,
        ...     vacancy_concentration=0.1,
        ...     placement="near_N",
        ...     vacancy_preference=0.7,  # Vacancies moderately prefer being near N
        ... )
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        # Step 1: O → N substitution
        o_indices = self.get_element_indices("O")
        n_oxygen = len(o_indices)
        n_to_substitute = int(n_oxygen * nitrogen_fraction)

        # Select O atoms to convert to N based on placement strategy
        # For N placement, use surface_n_preference
        substitute_indices = self._select_with_placement(
            o_indices, n_to_substitute, placement, preference=surface_n_preference,
            z_threshold=z_threshold
        )

        new_atoms = self.atoms.copy()
        symbols = list(new_atoms.get_chemical_symbols())

        for idx in substitute_indices:
            symbols[idx] = "N"

        new_atoms.set_chemical_symbols(symbols)

        # Step 2: Create vacancies
        if vacancy_concentration > 0:
            # Get indices of the element to create vacancies in
            if vacancy_element == "N":
                vacancy_candidates = substitute_indices  # The newly created N atoms
            else:
                vacancy_candidates = [
                    i for i in o_indices if i not in substitute_indices
                ]

            n_vacancies = int(len(vacancy_candidates) * vacancy_concentration)
            if n_vacancies == 0 and vacancy_concentration > 0:
                n_vacancies = 1

            # For near_N placement, get N positions
            n_positions = None
            if placement == "near_N":
                n_positions = new_atoms.positions[substitute_indices]

            # For vacancy placement, use vacancy_preference
            vacancy_indices = self._select_with_placement(
                vacancy_candidates, n_vacancies, placement,
                preference=vacancy_preference, n_positions=n_positions,
                z_threshold=z_threshold
            )

            # Remove vacancy atoms
            keep_indices = [i for i in range(len(new_atoms)) if i not in vacancy_indices]
            new_atoms = new_atoms[keep_indices]

        return new_atoms

    def create_oxynitride_pool(
        self,
        nitrogen_fraction: float = 0.67,
        vacancy_concentration: float = 0.0,
        vacancy_element: str = "N",
        n_configs_per_strategy: int = 3,
        strategies: Optional[List[str]] = None,
        surface_n_preference: float = 0.7,
        vacancy_preference: float = 0.7,
        z_threshold: float = 0.3,
        random_seed: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        """
        Generate a pool of oxynitride configurations with different placement strategies.

        Creates multiple configurations for each placement strategy to explore
        different defect distributions. Uses probability-weighted selection so
        each configuration is unique even with the same strategy.

        Parameters
        ----------
        nitrogen_fraction : float
            Fraction of O atoms to replace with N (0.0 to 1.0)
        vacancy_concentration : float
            Fraction of anion sites to leave vacant
        vacancy_element : str
            Element to create vacancies in ("N" or "O")
        n_configs_per_strategy : int
            Number of random configurations per strategy
        strategies : list, optional
            List of strategies to use. Default: ["random", "surface", "near_N"]
        surface_n_preference : float
            For "surface" strategy: preference for N at surface (0.5-1.0)
            0.5 = random, 0.7 = ~70% surface anions are N, 1.0 = strongly prefer surface
        vacancy_preference : float
            Preference strength for vacancy placement (0.5-1.0)
        z_threshold : float
            Fraction of slab height considered as surface region (0.3 = top 30%).
            Must match the value used in analyze_oxynitride_pool() for consistent
            analysis. Default: 0.3
        random_seed : int, optional
            Base random seed for reproducibility

        Returns
        -------
        list of dict
            List of configurations:
            [{
                "atoms": Atoms,
                "placement": str,
                "nitrogen_fraction": float,
                "vacancy_concentration": float,
                "surface_n_preference": float,
                "vacancy_preference": float,
                "z_threshold": float,
                "config_id": int,
            }, ...]

        Examples
        --------
        >>> builder = DefectBuilder(surface)
        >>> pool = builder.create_oxynitride_pool(
        ...     nitrogen_fraction=0.67,
        ...     vacancy_concentration=0.1,
        ...     n_configs_per_strategy=5,
        ...     surface_n_preference=0.8,  # 80% of surface anions will be N
        ...     vacancy_preference=0.6,    # Vacancies slightly prefer surface
        ... )
        >>> print(f"Generated {len(pool)} configurations")
        """
        if strategies is None:
            strategies = ["random", "surface", "near_N"]

        results = []
        config_id = 0

        for strategy in strategies:
            for i in range(n_configs_per_strategy):
                seed = None
                if random_seed is not None:
                    seed = random_seed + config_id

                atoms = self.create_oxynitride(
                    nitrogen_fraction=nitrogen_fraction,
                    vacancy_concentration=vacancy_concentration,
                    vacancy_element=vacancy_element,
                    placement=strategy,
                    surface_n_preference=surface_n_preference,
                    vacancy_preference=vacancy_preference,
                    z_threshold=z_threshold,
                    random_seed=seed,
                )

                results.append({
                    "atoms": atoms,
                    "placement": strategy,
                    "nitrogen_fraction": nitrogen_fraction,
                    "vacancy_concentration": vacancy_concentration,
                    "surface_n_preference": surface_n_preference,
                    "vacancy_preference": vacancy_preference,
                    "z_threshold": z_threshold,
                    "config_id": config_id,
                })
                config_id += 1

        return results

    def create_multiple_oxynitrides(
        self,
        nitrogen_fractions: List[float],
        vacancy_concentrations: List[float],
        n_configs_per_combo: int = 3,
        random_seed: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        """
        Generate multiple oxynitride configurations for screening.

        Parameters
        ----------
        nitrogen_fractions : list
            List of N fractions to try
        vacancy_concentrations : list
            List of vacancy concentrations (x values)
        n_configs_per_combo : int
            Number of random configurations per combination
        random_seed : int, optional
            Base random seed

        Returns
        -------
        list
            List of dicts with 'atoms', 'n_fraction', 'vac_conc', 'config_id'
        """
        results = []
        config_id = 0

        for n_frac in nitrogen_fractions:
            for vac_conc in vacancy_concentrations:
                for i in range(n_configs_per_combo):
                    seed = None
                    if random_seed is not None:
                        seed = random_seed + config_id

                    atoms = self.create_oxynitride(
                        nitrogen_fraction=n_frac,
                        vacancy_concentration=vac_conc,
                        random_seed=seed,
                    )

                    results.append({
                        "atoms": atoms,
                        "nitrogen_fraction": n_frac,
                        "vacancy_concentration": vac_conc,
                        "config_id": config_id,
                    })
                    config_id += 1

        return results

    def create_surface_vacancy(
        self,
        element: str,
        z_threshold: float = 0.2,
        random_seed: Optional[int] = None,
    ) -> Atoms:
        """
        Create a vacancy only in the surface layer.

        Parameters
        ----------
        element : str
            Element to remove
        z_threshold : float
            Fraction of slab considered as surface
        random_seed : int, optional
            Random seed

        Returns
        -------
        Atoms
            Structure with surface vacancy
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        # Find surface atoms
        z_positions = self.atoms.positions[:, 2]
        z_max = z_positions.max()
        z_min = z_positions.min()
        z_cutoff = z_max - z_threshold * (z_max - z_min)

        # Get element indices in surface region
        element_indices = self.get_element_indices(element)
        surface_element_indices = [
            i for i in element_indices
            if self.atoms.positions[i, 2] > z_cutoff
        ]

        if not surface_element_indices:
            raise ValueError(f"No {element} atoms found in surface region")

        # Remove one random surface atom
        vacancy_idx = np.random.choice(surface_element_indices)
        keep_indices = [i for i in range(len(self.atoms)) if i != vacancy_idx]

        return self.atoms[keep_indices]

    def get_defect_info(self) -> Dict[str, Any]:
        """
        Get information about the current structure composition.

        Returns
        -------
        dict
            Element counts and ratios
        """
        symbols = self.atoms.get_chemical_symbols()
        counts = {}
        for s in symbols:
            counts[s] = counts.get(s, 0) + 1

        total = len(symbols)
        ratios = {k: v / total for k, v in counts.items()}

        # Calculate anion ratio (O:N) if both present
        anion_ratio = None
        if "O" in counts and "N" in counts:
            total_anions = counts["O"] + counts["N"]
            anion_ratio = {
                "O": counts["O"] / total_anions,
                "N": counts["N"] / total_anions,
            }

        return {
            "n_atoms": total,
            "element_counts": counts,
            "element_ratios": ratios,
            "anion_ratio": anion_ratio,
            "formula": self.atoms.get_chemical_formula(),
        }


class OxynitrideStructure(BaseStructure):
    """
    Structure class specifically for oxynitride materials.

    Tracks nitrogen fraction and vacancy concentration.
    """

    def __init__(
        self,
        atoms: Atoms,
        nitrogen_fraction: Optional[float] = None,
        vacancy_concentration: Optional[float] = None,
        parent_formula: Optional[str] = None,
    ):
        """
        Initialize OxynitrideStructure.

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object
        nitrogen_fraction : float, optional
            Fraction of O replaced by N
        vacancy_concentration : float, optional
            Vacancy concentration (x parameter)
        parent_formula : str, optional
            Original oxide formula (e.g., "LaVO3")
        """
        super().__init__(atoms)
        self.nitrogen_fraction = nitrogen_fraction
        self.vacancy_concentration = vacancy_concentration
        self.parent_formula = parent_formula

    @classmethod
    def from_oxide(
        cls,
        oxide: Union[Atoms, BaseStructure],
        nitrogen_fraction: float = 0.67,
        vacancy_concentration: float = 0.0,
        random_seed: Optional[int] = None,
    ) -> "OxynitrideStructure":
        """
        Create oxynitride from an oxide structure.

        Parameters
        ----------
        oxide : Atoms or BaseStructure
            Parent oxide structure
        nitrogen_fraction : float
            Fraction of O to replace with N
        vacancy_concentration : float
            Vacancy concentration
        random_seed : int, optional
            Random seed

        Returns
        -------
        OxynitrideStructure
            New oxynitride structure
        """
        if isinstance(oxide, BaseStructure):
            atoms = oxide.atoms.copy()
            parent_formula = oxide.formula
        else:
            atoms = oxide.copy()
            parent_formula = atoms.get_chemical_formula()

        builder = DefectBuilder(atoms)
        oxynitride_atoms = builder.create_oxynitride(
            nitrogen_fraction=nitrogen_fraction,
            vacancy_concentration=vacancy_concentration,
            random_seed=random_seed,
        )

        return cls(
            oxynitride_atoms,
            nitrogen_fraction=nitrogen_fraction,
            vacancy_concentration=vacancy_concentration,
            parent_formula=parent_formula,
        )

    def get_stoichiometry(self) -> Dict[str, float]:
        """
        Calculate the stoichiometry relative to the cation.

        For LaVON₂₋ₓ, returns ratios relative to V.

        Returns
        -------
        dict
            Stoichiometric ratios
        """
        counts = {}
        for symbol in self.atoms.get_chemical_symbols():
            counts[symbol] = counts.get(symbol, 0) + 1

        # Find the transition metal (likely the catalyst site)
        transition_metals = ["V", "Ti", "Mn", "Fe", "Co", "Ni", "Cu", "Zr"]
        reference = None
        for tm in transition_metals:
            if tm in counts:
                reference = tm
                break

        if reference is None:
            # Fall back to first non-O, non-N element
            for symbol in counts:
                if symbol not in ["O", "N"]:
                    reference = symbol
                    break

        if reference is None:
            reference = list(counts.keys())[0]

        ref_count = counts[reference]
        stoich = {k: v / ref_count for k, v in counts.items()}

        return stoich

    def __repr__(self) -> str:
        return (
            f"OxynitrideStructure({self.formula}, "
            f"N_frac={self.nitrogen_fraction}, "
            f"vac_x={self.vacancy_concentration})"
        )


def generate_vacancy_series(
    structure: Union[Atoms, BaseStructure],
    vacancy_concentrations: List[float],
    element: str = "O",
    n_configs: int = 3,
    random_seed: Optional[int] = None,
) -> List[Dict[str, Any]]:
    """
    Generate a series of structures with different vacancy concentrations.

    Parameters
    ----------
    structure : Atoms or BaseStructure
        Base structure
    vacancy_concentrations : list
        List of vacancy concentrations to generate
    element : str
        Element to create vacancies in
    n_configs : int
        Number of random configurations per concentration
    random_seed : int, optional
        Base random seed

    Returns
    -------
    list
        List of dicts with 'atoms', 'concentration', 'config_id'
    """
    builder = DefectBuilder(structure)
    results = []
    config_id = 0

    for conc in vacancy_concentrations:
        for i in range(n_configs):
            seed = None
            if random_seed is not None:
                seed = random_seed + config_id

            atoms = builder.create_vacancy(
                element=element,
                concentration=conc,
                random_seed=seed,
            )

            results.append({
                "atoms": atoms,
                "concentration": conc,
                "config_id": config_id,
                "n_vacancies": builder.atoms.get_chemical_symbols().count(element)
                             - atoms.get_chemical_symbols().count(element),
            })
            config_id += 1

    return results


def analyze_defect_distribution(
    atoms: Atoms,
    reference_atoms: Optional[Atoms] = None,
    z_threshold: float = 0.3,
    near_n_cutoff: float = 3.0,
) -> Dict[str, Any]:
    """
    Analyze the distribution of N atoms and vacancies in a structure.

    Parameters
    ----------
    atoms : Atoms
        Structure to analyze (oxynitride with possible vacancies)
    reference_atoms : Atoms, optional
        Original structure before vacancy creation (to count vacancies)
    z_threshold : float
        Fraction of slab height considered as surface region (default: 0.3 = top 30%)
    near_n_cutoff : float
        Distance cutoff for "near N" analysis in Angstroms (default: 3.0)

    Returns
    -------
    dict
        Analysis results with keys:
        - n_total: total number of N atoms
        - n_surface: N atoms in surface region
        - n_surface_fraction: fraction of N in surface region
        - o_total: total number of O atoms
        - o_surface: O atoms in surface region
        - surface_n_ratio: N/(N+O) ratio in surface region
        - bulk_n_ratio: N/(N+O) ratio in bulk region
        - vacancy_total: total vacancies (if reference provided)
        - vacancy_surface: vacancies in surface region
        - vacancy_near_n: vacancies within cutoff of N atoms
        - vacancy_near_n_fraction: fraction of vacancies near N

    Examples
    --------
    >>> stats = analyze_defect_distribution(oxynitride, reference_atoms=original)
    >>> print(f"Surface N ratio: {stats['surface_n_ratio']:.1%}")
    >>> print(f"Vacancies near N: {stats['vacancy_near_n_fraction']:.1%}")
    """
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    z_positions = positions[:, 2]

    # Define surface region (top z_threshold fraction of slab)
    z_min, z_max = z_positions.min(), z_positions.max()
    z_cutoff = z_max - z_threshold * (z_max - z_min)

    # Get indices by element
    n_indices = [i for i, s in enumerate(symbols) if s == "N"]
    o_indices = [i for i, s in enumerate(symbols) if s == "O"]

    # Count N and O in surface vs bulk
    n_surface = [i for i in n_indices if z_positions[i] >= z_cutoff]
    n_bulk = [i for i in n_indices if z_positions[i] < z_cutoff]
    o_surface = [i for i in o_indices if z_positions[i] >= z_cutoff]
    o_bulk = [i for i in o_indices if z_positions[i] < z_cutoff]

    # Calculate ratios
    total_surface_anions = len(n_surface) + len(o_surface)
    total_bulk_anions = len(n_bulk) + len(o_bulk)

    surface_n_ratio = len(n_surface) / total_surface_anions if total_surface_anions > 0 else 0
    bulk_n_ratio = len(n_bulk) / total_bulk_anions if total_bulk_anions > 0 else 0

    results = {
        "n_total": len(n_indices),
        "n_surface": len(n_surface),
        "n_bulk": len(n_bulk),
        "n_surface_fraction": len(n_surface) / len(n_indices) if n_indices else 0,
        "o_total": len(o_indices),
        "o_surface": len(o_surface),
        "o_bulk": len(o_bulk),
        "surface_n_ratio": surface_n_ratio,
        "bulk_n_ratio": bulk_n_ratio,
        "z_threshold": z_threshold,
        "z_cutoff": z_cutoff,
    }

    # Vacancy analysis (requires reference structure)
    if reference_atoms is not None:
        ref_symbols = reference_atoms.get_chemical_symbols()
        ref_positions = reference_atoms.get_positions()

        # Find vacancy positions by comparing with reference
        # Vacancies are anion sites in reference that are missing in current
        ref_anion_positions = []
        for i, s in enumerate(ref_symbols):
            if s in ["O", "N"]:
                ref_anion_positions.append(ref_positions[i])

        current_anion_positions = []
        for i, s in enumerate(symbols):
            if s in ["O", "N"]:
                current_anion_positions.append(positions[i])

        # Find missing positions (vacancies)
        vacancy_positions = []
        for ref_pos in ref_anion_positions:
            found = False
            for cur_pos in current_anion_positions:
                if np.linalg.norm(ref_pos - cur_pos) < 0.5:  # Within 0.5 Å
                    found = True
                    break
            if not found:
                vacancy_positions.append(ref_pos)

        vacancy_positions = np.array(vacancy_positions) if vacancy_positions else np.array([]).reshape(0, 3)

        # Count vacancies in surface region
        if len(vacancy_positions) > 0:
            vacancy_z = vacancy_positions[:, 2]
            vacancy_surface = np.sum(vacancy_z >= z_cutoff)
            vacancy_bulk = np.sum(vacancy_z < z_cutoff)

            # Count vacancies near N atoms
            n_positions = positions[n_indices]
            vacancy_near_n = 0
            for vac_pos in vacancy_positions:
                if len(n_positions) > 0:
                    min_dist = np.min(np.linalg.norm(n_positions - vac_pos, axis=1))
                    if min_dist <= near_n_cutoff:
                        vacancy_near_n += 1

            results.update({
                "vacancy_total": len(vacancy_positions),
                "vacancy_surface": int(vacancy_surface),
                "vacancy_bulk": int(vacancy_bulk),
                "vacancy_surface_fraction": vacancy_surface / len(vacancy_positions) if len(vacancy_positions) > 0 else 0,
                "vacancy_near_n": vacancy_near_n,
                "vacancy_near_n_fraction": vacancy_near_n / len(vacancy_positions) if len(vacancy_positions) > 0 else 0,
                "near_n_cutoff": near_n_cutoff,
            })
        else:
            results.update({
                "vacancy_total": 0,
                "vacancy_surface": 0,
                "vacancy_bulk": 0,
                "vacancy_surface_fraction": 0,
                "vacancy_near_n": 0,
                "vacancy_near_n_fraction": 0,
                "near_n_cutoff": near_n_cutoff,
            })

    return results


def analyze_oxynitride_pool(
    pool: List[Dict[str, Any]],
    reference_atoms: Optional[Atoms] = None,
    z_threshold: float = 0.3,
    near_n_cutoff: float = 3.0,
) -> Dict[str, Any]:
    """
    Analyze a pool of oxynitride configurations and compute statistics.

    Parameters
    ----------
    pool : list of dict
        Pool from create_oxynitride_pool(), each with "atoms" and "placement" keys
    reference_atoms : Atoms, optional
        Original structure before defects (for vacancy counting)
    z_threshold : float
        Fraction of slab considered as surface (default: 0.3 = top 30%)
    near_n_cutoff : float
        Distance cutoff for "near N" analysis (default: 3.0 Å)

    Returns
    -------
    dict
        Summary statistics:
        - by_strategy: dict with stats for each placement strategy
        - overall: aggregated stats across all configurations

    Examples
    --------
    >>> pool = defect.create_oxynitride_pool(...)
    >>> stats = analyze_oxynitride_pool(pool, reference_atoms=slab.atoms)
    >>>
    >>> # Print summary by strategy
    >>> for strategy, data in stats["by_strategy"].items():
    ...     print(f"{strategy}:")
    ...     print(f"  Surface N ratio: {data['surface_n_ratio_mean']:.1%} ± {data['surface_n_ratio_std']:.1%}")
    """
    # Group by placement strategy
    by_strategy = {}

    for config in pool:
        strategy = config.get("placement", "unknown")
        atoms = config["atoms"]

        analysis = analyze_defect_distribution(
            atoms,
            reference_atoms=reference_atoms,
            z_threshold=z_threshold,
            near_n_cutoff=near_n_cutoff,
        )

        if strategy not in by_strategy:
            by_strategy[strategy] = []
        by_strategy[strategy].append(analysis)

    # Compute statistics for each strategy
    results = {"by_strategy": {}, "overall": {}}

    all_analyses = []

    for strategy, analyses in by_strategy.items():
        all_analyses.extend(analyses)

        strategy_stats = {
            "n_configs": len(analyses),
            "surface_n_ratio_mean": np.mean([a["surface_n_ratio"] for a in analyses]),
            "surface_n_ratio_std": np.std([a["surface_n_ratio"] for a in analyses]),
            "bulk_n_ratio_mean": np.mean([a["bulk_n_ratio"] for a in analyses]),
            "bulk_n_ratio_std": np.std([a["bulk_n_ratio"] for a in analyses]),
            "n_surface_fraction_mean": np.mean([a["n_surface_fraction"] for a in analyses]),
            "n_surface_fraction_std": np.std([a["n_surface_fraction"] for a in analyses]),
        }

        # Add vacancy stats if available
        if "vacancy_total" in analyses[0]:
            strategy_stats.update({
                "vacancy_surface_fraction_mean": np.mean([a["vacancy_surface_fraction"] for a in analyses]),
                "vacancy_surface_fraction_std": np.std([a["vacancy_surface_fraction"] for a in analyses]),
                "vacancy_near_n_fraction_mean": np.mean([a["vacancy_near_n_fraction"] for a in analyses]),
                "vacancy_near_n_fraction_std": np.std([a["vacancy_near_n_fraction"] for a in analyses]),
            })

        results["by_strategy"][strategy] = strategy_stats

    # Overall statistics
    if all_analyses:
        results["overall"] = {
            "n_configs": len(all_analyses),
            "surface_n_ratio_mean": np.mean([a["surface_n_ratio"] for a in all_analyses]),
            "surface_n_ratio_std": np.std([a["surface_n_ratio"] for a in all_analyses]),
            "n_surface_fraction_mean": np.mean([a["n_surface_fraction"] for a in all_analyses]),
        }

        if "vacancy_total" in all_analyses[0]:
            results["overall"].update({
                "vacancy_surface_fraction_mean": np.mean([a["vacancy_surface_fraction"] for a in all_analyses]),
                "vacancy_near_n_fraction_mean": np.mean([a["vacancy_near_n_fraction"] for a in all_analyses]),
            })

    return results


def print_defect_analysis(
    stats: Dict[str, Any],
    title: str = "Defect Distribution Analysis",
) -> None:
    """
    Print a formatted summary of defect distribution analysis.

    Parameters
    ----------
    stats : dict
        Results from analyze_defect_distribution() or analyze_oxynitride_pool()
    title : str
        Title for the report

    Examples
    --------
    >>> stats = analyze_defect_distribution(oxynitride)
    >>> print_defect_analysis(stats)

    >>> pool_stats = analyze_oxynitride_pool(pool)
    >>> print_defect_analysis(pool_stats, title="Pool Analysis")
    """
    print("=" * 60)
    print(title)
    print("=" * 60)

    # Single structure analysis
    if "n_total" in stats:
        print(f"\nNitrogen Distribution:")
        print(f"  Total N atoms: {stats['n_total']}")
        print(f"  N in surface (top {stats['z_threshold']*100:.0f}%): {stats['n_surface']} ({stats['n_surface_fraction']*100:.1f}%)")
        print(f"  N in bulk: {stats['n_bulk']}")
        print(f"\nAnion Composition by Region:")
        print(f"  Surface N/(N+O) ratio: {stats['surface_n_ratio']*100:.1f}%")
        print(f"  Bulk N/(N+O) ratio: {stats['bulk_n_ratio']*100:.1f}%")

        if "vacancy_total" in stats:
            print(f"\nVacancy Distribution:")
            print(f"  Total vacancies: {stats['vacancy_total']}")
            print(f"  Vacancies in surface: {stats['vacancy_surface']} ({stats['vacancy_surface_fraction']*100:.1f}%)")
            print(f"  Vacancies in bulk: {stats['vacancy_bulk']}")
            print(f"  Vacancies near N (< {stats['near_n_cutoff']} Å): {stats['vacancy_near_n']} ({stats['vacancy_near_n_fraction']*100:.1f}%)")

    # Pool analysis
    elif "by_strategy" in stats:
        print(f"\nTotal configurations: {stats['overall']['n_configs']}")

        for strategy, data in stats["by_strategy"].items():
            print(f"\n{strategy.upper()} Strategy ({data['n_configs']} configs):")
            print(f"  Surface N ratio: {data['surface_n_ratio_mean']*100:.1f}% ± {data['surface_n_ratio_std']*100:.1f}%")
            print(f"  Bulk N ratio: {data['bulk_n_ratio_mean']*100:.1f}% ± {data['bulk_n_ratio_std']*100:.1f}%")
            print(f"  N in surface: {data['n_surface_fraction_mean']*100:.1f}% ± {data['n_surface_fraction_std']*100:.1f}%")

            if "vacancy_surface_fraction_mean" in data:
                print(f"  Vacancies in surface: {data['vacancy_surface_fraction_mean']*100:.1f}% ± {data['vacancy_surface_fraction_std']*100:.1f}%")
                print(f"  Vacancies near N: {data['vacancy_near_n_fraction_mean']*100:.1f}% ± {data['vacancy_near_n_fraction_std']*100:.1f}%")

    print("=" * 60)
