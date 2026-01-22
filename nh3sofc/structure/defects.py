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
        surface_fraction: float = 0.3,
    ) -> Tuple[List[int], List[int]]:
        """
        Split indices into surface and bulk based on z-position.

        Parameters
        ----------
        indices : list
            Atom indices to split
        surface_fraction : float
            Fraction of slab height considered as surface region

        Returns
        -------
        tuple
            (surface_indices, bulk_indices)
        """
        z_positions = self.atoms.positions[:, 2]
        z_min, z_max = z_positions.min(), z_positions.max()
        z_cutoff = z_max - surface_fraction * (z_max - z_min)

        surface_indices = [i for i in indices if z_positions[i] >= z_cutoff]
        bulk_indices = [i for i in indices if z_positions[i] < z_cutoff]

        return surface_indices, bulk_indices

    def _calculate_selection_weights(
        self,
        candidates: List[int],
        placement: str = "random",
        preference: float = 0.7,
        n_positions: Optional[np.ndarray] = None,
        surface_fraction: float = 0.3,
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
        surface_fraction : float
            Fraction of slab height considered as surface region

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
            candidates, placement, preference, n_positions
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
            o_indices, n_to_substitute, placement, preference=surface_n_preference
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
                preference=vacancy_preference, n_positions=n_positions
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
                    random_seed=seed,
                )

                results.append({
                    "atoms": atoms,
                    "placement": strategy,
                    "nitrogen_fraction": nitrogen_fraction,
                    "vacancy_concentration": vacancy_concentration,
                    "surface_n_preference": surface_n_preference,
                    "vacancy_preference": vacancy_preference,
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
