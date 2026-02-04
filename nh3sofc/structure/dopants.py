"""Acceptor dopant module for ceria-based materials.

Creates doped ceria structures (GDC, SDC, PDC) with charge-compensating
oxygen vacancies for solid oxide fuel cell applications.

Chemistry:
- Trivalent dopants (Sm³⁺, Gd³⁺, Pr³⁺) substitute Ce⁴⁺ sites as acceptor dopants
- Charge compensation requires oxygen vacancy formation
- 2 Ce⁴⁺ → 2 M³⁺ + V_O^{2+} (2 dopants : 1 vacancy ratio)
- Formula: Ce₁₋₂ₓM₂ₓO₂₋ₓ (where x = dopant mole fraction from M₂O₃)

Examples
--------
>>> from nh3sofc.structure import SlabStructure, DopantBuilder
>>> slab = SlabStructure.from_file("CeO2_111.cif")
>>> builder = DopantBuilder(slab)
>>> gdc = builder.create_doped_structure(
...     dopant="Gd",
...     dopant_fraction=0.10,  # 10% Gd on Ce sites
... )
"""

from pathlib import Path
from typing import Optional, List, Union, Tuple, Dict, Any, TYPE_CHECKING

import numpy as np
from ase import Atoms

from ..core.base import BaseStructure
from ..core.constants import ACCEPTOR_DOPANTS, HOST_CATIONS, CHARGE_COMPENSATION
from .surface import SlabStructure


class DopantBuilder:
    """
    Builder class for creating acceptor-doped oxide structures with charge compensation.

    Supports doping of CeO2 with trivalent dopants (Gd, Sm, Pr, Y, La, Nd) and
    automatic creation of charge-compensating oxygen vacancies.

    Examples
    --------
    >>> slab = SlabStructure.from_file("CeO2_111.cif")
    >>> builder = DopantBuilder(slab)

    >>> # Create 10% Gd-doped ceria with random placement
    >>> gdc = builder.create_doped_structure(
    ...     dopant="Gd",
    ...     dopant_fraction=0.10,
    ... )

    >>> # Create with surface-preferring vacancies
    >>> gdc_surface = builder.create_doped_structure(
    ...     dopant="Gd",
    ...     dopant_fraction=0.10,
    ...     vacancy_placement="surface",
    ...     vacancy_preference=0.8,
    ... )

    >>> # Generate pool of configurations for screening
    >>> pool = builder.create_doped_pool(
    ...     dopant="Sm",
    ...     dopant_fraction=0.15,
    ...     n_configs_per_strategy=5,
    ... )
    """

    def __init__(self, structure: Union[SlabStructure, Atoms]):
        """
        Initialize DopantBuilder.

        Parameters
        ----------
        structure : SlabStructure or Atoms
            CeO2 structure (bulk or slab) to dope
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
        dopant_positions: Optional[np.ndarray] = None,
        z_threshold: float = 0.3,
    ) -> np.ndarray:
        """
        Calculate selection weights for each candidate based on placement strategy.

        Parameters
        ----------
        candidates : list
            Candidate atom indices
        placement : str
            Placement strategy: "random", "surface", "bulk", or "near_dopant"
        preference : float
            Preference strength (0.5 = random, 1.0 = strongly prefer target region)
        dopant_positions : ndarray, optional
            Positions of dopant atoms (required for "near_dopant" placement)
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
            steepness = (preference - 0.5) * 10  # 0 to 5
            weights = np.exp(steepness * z_normalized)

            return weights / weights.sum()

        elif placement == "bulk":
            # Weight by z-position: lower z = higher weight (inverse of surface)
            z_positions = self.atoms.positions[:, 2]
            z_min = z_positions.min()
            z_max = z_positions.max()
            z_range = z_max - z_min

            if z_range < 0.1:
                return np.ones(n_candidates) / n_candidates

            # Normalize z to [0, 1] for candidates, then invert
            z_values = np.array([z_positions[i] for i in candidates])
            z_normalized = 1.0 - (z_values - z_min) / z_range

            # Apply preference
            steepness = (preference - 0.5) * 10
            weights = np.exp(steepness * z_normalized)

            return weights / weights.sum()

        elif placement == "near_dopant" and dopant_positions is not None and len(dopant_positions) > 0:
            # Weight by inverse distance to nearest dopant
            positions = self.atoms.positions
            distances = []
            for idx in candidates:
                pos = positions[idx]
                min_dist = np.min(np.linalg.norm(dopant_positions - pos, axis=1))
                distances.append(max(min_dist, 0.1))  # Avoid division by zero

            distances = np.array(distances)

            # Apply preference: closer to dopant = higher weight
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
        dopant_positions: Optional[np.ndarray] = None,
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
            Placement strategy: "random", "surface", "bulk", or "near_dopant"
        preference : float
            Preference strength (0.5 = random, 1.0 = strongly prefer target region)
        dopant_positions : ndarray, optional
            Positions of dopant atoms (required for "near_dopant" placement)
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
            candidates, placement, preference, dopant_positions, z_threshold
        )

        # Weighted random selection without replacement
        selected = np.random.choice(
            candidates,
            size=n_select,
            replace=False,
            p=weights
        ).tolist()

        return selected

    def create_doped_structure(
        self,
        dopant: str,
        dopant_fraction: float,
        host_cation: str = "Ce",
        auto_compensate: bool = True,
        vacancy_placement: str = "random",
        dopant_placement: str = "random",
        dopant_preference: float = 0.7,
        vacancy_preference: float = 0.7,
        pr_trivalent_fraction: float = 1.0,
        z_threshold: float = 0.3,
        random_seed: Optional[int] = None,
    ) -> Atoms:
        """
        Create doped structure with charge-compensating vacancies.

        Parameters
        ----------
        dopant : str
            Dopant element ("Sm", "Gd", "Pr", "Y", "La", "Nd")
        dopant_fraction : float
            Fraction of host cations to replace (0.0 to 1.0)
        host_cation : str
            Cation to replace (default: "Ce")
        auto_compensate : bool
            Create compensating vacancies automatically (default: True)
        vacancy_placement : str
            "random" ignores preference parameter (uniform distribution)
            "surface" uses preference to weight toward high-z atoms
            "near_dopant" uses preference to weight toward dopant sites
            "bulk" uses preference to weight toward low-z atoms
        dopant_placement : str
            "random" ignores preference parameter (uniform distribution)
            "surface" uses preference to weight toward high-z atoms
            "bulk" uses preference to weight toward low-z atoms
        dopant_preference : float
            Bias strength for dopant placement: 0.5 = uniform random,
            0.7 = moderate bias, 1.0 = strong bias
            Only applies when dopant_placement != "random"
        vacancy_preference : float
            Bias strength for vacancy placement: 0.5 = uniform random,
            0.7 = moderate bias, 1.0 = strong bias
            Only applies when vacancy_placement != "random"
        pr_trivalent_fraction : float
            For Pr only: fraction that is Pr³⁺ (rest is Pr⁴⁺, no vacancy contribution)
            Vacancy count = floor(n_Pr * pr_trivalent_fraction / 2)
        z_threshold : float
            Fraction of slab considered "surface" (default: 0.3 = top 30%)
        random_seed : int, optional
            Random seed for reproducibility

        Returns
        -------
        Atoms
            Doped structure with charge-compensating vacancies

        Examples
        --------
        >>> # 10% Gd-doped ceria with random placement
        >>> gdc = builder.create_doped_structure(
        ...     dopant="Gd",
        ...     dopant_fraction=0.10,
        ... )

        >>> # GDC with vacancies preferring dopant sites
        >>> gdc_associated = builder.create_doped_structure(
        ...     dopant="Gd",
        ...     dopant_fraction=0.10,
        ...     vacancy_placement="near_dopant",
        ...     vacancy_preference=0.8,
        ... )

        >>> # Pr-doped ceria with mixed valence
        >>> pdc = builder.create_doped_structure(
        ...     dopant="Pr",
        ...     dopant_fraction=0.20,
        ...     pr_trivalent_fraction=0.5,  # 50% Pr³⁺, 50% Pr⁴⁺
        ... )
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        # Validate dopant
        if dopant not in ACCEPTOR_DOPANTS:
            available = ", ".join(ACCEPTOR_DOPANTS.keys())
            raise ValueError(f"Unknown dopant '{dopant}'. Available: {available}")

        # Get host cation indices
        host_indices = self.get_element_indices(host_cation)
        if not host_indices:
            raise ValueError(f"No {host_cation} atoms found in structure")

        n_host = len(host_indices)
        n_dopants = int(n_host * dopant_fraction)

        if n_dopants == 0 and dopant_fraction > 0:
            n_dopants = 1  # At least one dopant if fraction > 0

        # Step 1: Select and substitute host cations with dopant
        dopant_indices = self._select_with_placement(
            host_indices,
            n_dopants,
            placement=dopant_placement,
            preference=dopant_preference,
            z_threshold=z_threshold,
        )

        new_atoms = self.atoms.copy()
        symbols = list(new_atoms.get_chemical_symbols())

        for idx in dopant_indices:
            symbols[idx] = dopant

        new_atoms.set_chemical_symbols(symbols)

        # Step 2: Create charge-compensating vacancies
        if auto_compensate and n_dopants > 0:
            dopant_charge = ACCEPTOR_DOPANTS[dopant]["charge"]
            n_dopants_per_vacancy, n_vacancies_per_unit = CHARGE_COMPENSATION.get(
                dopant_charge, (2, 1)
            )

            # Calculate effective dopants for vacancy calculation
            effective_dopants = n_dopants
            if dopant == "Pr":
                # Pr can be Pr³⁺ or Pr⁴⁺; only Pr³⁺ contributes to vacancies
                effective_dopants = int(n_dopants * pr_trivalent_fraction)

            # Calculate number of vacancies needed
            n_vacancies = (effective_dopants * n_vacancies_per_unit) // n_dopants_per_vacancy

            if n_vacancies > 0:
                # Get oxygen indices
                o_indices = [
                    i for i, s in enumerate(new_atoms.get_chemical_symbols())
                    if s == "O"
                ]

                if not o_indices:
                    raise ValueError("No oxygen atoms found for vacancy creation")

                # Get dopant positions for near_dopant placement
                dopant_positions = None
                if vacancy_placement == "near_dopant":
                    dopant_positions = new_atoms.positions[dopant_indices]

                # Select vacancy positions
                vacancy_indices = self._select_with_placement(
                    o_indices,
                    n_vacancies,
                    placement=vacancy_placement,
                    preference=vacancy_preference,
                    dopant_positions=dopant_positions,
                    z_threshold=z_threshold,
                )

                # Remove vacancy atoms
                keep_indices = [
                    i for i in range(len(new_atoms))
                    if i not in vacancy_indices
                ]
                new_atoms = new_atoms[keep_indices]

        return new_atoms

    def create_doped_pool(
        self,
        dopant: str,
        dopant_fraction: Union[float, List[float]],
        n_configs: int = 3,
        strategies: Optional[List[str]] = None,
        host_cation: str = "Ce",
        dopant_placement: str = "random",
        dopant_preference: float = 0.7,
        vacancy_preference: float = 0.7,
        pr_trivalent_fraction: float = 1.0,
        z_threshold: float = 0.3,
        random_seed: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        """
        Generate multiple doped configurations for screening.

        Creates a pool of doped structures by varying:
        - Dopant fraction (if a list is provided)
        - Vacancy placement strategy
        - Random configurations per combination

        Parameters
        ----------
        dopant : str
            Dopant element ("Sm", "Gd", "Pr", "Y", "La", "Nd")
        dopant_fraction : float or list of float
            Fraction(s) of host cations to replace. Can be a single value
            or a list for concentration series (e.g., [0.05, 0.10, 0.15, 0.20])
        n_configs : int
            Number of random configurations per (fraction, strategy) combination.
            Default: 3
        strategies : list, optional
            List of vacancy placement strategies.
            Default: ["random"] if dopant_fraction is a list, else
            ["random", "surface", "near_dopant"]
        host_cation : str
            Cation to replace (default: "Ce")
        dopant_placement : str
            Dopant placement strategy (default: "random")
        dopant_preference : float
            Dopant placement preference strength (default: 0.7)
        vacancy_preference : float
            Vacancy placement preference strength (default: 0.7)
        pr_trivalent_fraction : float
            For Pr: fraction that is Pr³⁺ (default: 1.0)
        z_threshold : float
            Fraction of slab considered surface (default: 0.3)
        random_seed : int, optional
            Base random seed for reproducibility

        Returns
        -------
        list of dict
            List of configurations:
            [{
                "atoms": Atoms,
                "dopant": str,
                "dopant_fraction": float,
                "vacancy_placement": str,
                "dopant_placement": str,
                "config_id": int,
            }, ...]

        Examples
        --------
        >>> # Single fraction, multiple strategies
        >>> pool = builder.create_doped_pool(
        ...     dopant="Sm",
        ...     dopant_fraction=0.15,
        ...     n_configs=5,
        ...     strategies=["random", "surface", "near_dopant"],
        ... )
        >>> print(f"Generated {len(pool)} configurations")  # 15

        >>> # Multiple fractions (concentration series)
        >>> pool = builder.create_doped_pool(
        ...     dopant="Gd",
        ...     dopant_fraction=[0.05, 0.10, 0.15, 0.20],
        ...     n_configs=3,
        ... )
        >>> print(f"Generated {len(pool)} configurations")  # 12
        """
        # Normalize dopant_fraction to list
        if isinstance(dopant_fraction, (int, float)):
            fractions = [float(dopant_fraction)]
        else:
            fractions = list(dopant_fraction)

        # Default strategies: if multiple fractions, just use random
        # (concentration series typically don't need strategy variation)
        if strategies is None:
            if len(fractions) > 1:
                strategies = ["random"]
            else:
                strategies = ["random", "surface", "near_dopant"]

        results = []
        config_id = 0

        for frac in fractions:
            for strategy in strategies:
                for i in range(n_configs):
                    seed = None
                    if random_seed is not None:
                        seed = random_seed + config_id

                    atoms = self.create_doped_structure(
                        dopant=dopant,
                        dopant_fraction=frac,
                        host_cation=host_cation,
                        vacancy_placement=strategy,
                        dopant_placement=dopant_placement,
                        dopant_preference=dopant_preference,
                        vacancy_preference=vacancy_preference,
                        pr_trivalent_fraction=pr_trivalent_fraction,
                        z_threshold=z_threshold,
                        random_seed=seed,
                    )

                    results.append({
                        "atoms": atoms,
                        "dopant": dopant,
                        "dopant_fraction": frac,
                        "vacancy_placement": strategy,
                        "dopant_placement": dopant_placement,
                        "dopant_preference": dopant_preference,
                        "vacancy_preference": vacancy_preference,
                        "z_threshold": z_threshold,
                        "config_id": config_id,
                    })
                    config_id += 1

        return results

    def get_dopant_info(self) -> Dict[str, Any]:
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

        # Identify dopants present
        dopants_present = [s for s in counts.keys() if s in ACCEPTOR_DOPANTS]

        return {
            "n_atoms": total,
            "element_counts": counts,
            "element_ratios": ratios,
            "dopants_present": dopants_present,
            "formula": self.atoms.get_chemical_formula(),
        }


class DopedCeriaStructure(BaseStructure):
    """
    Structure class for doped ceria materials.

    Tracks dopant and vacancy information for GDC, SDC, PDC, etc.

    Examples
    --------
    >>> doped = DopedCeriaStructure.from_ceria(
    ...     ceria_slab,
    ...     dopant="Gd",
    ...     dopant_fraction=0.10,
    ... )
    >>> print(doped.get_stoichiometry())
    """

    def __init__(
        self,
        atoms: Atoms,
        dopant: str,
        dopant_fraction: float,
        n_vacancies: int,
        parent_formula: Optional[str] = None,
    ):
        """
        Initialize DopedCeriaStructure.

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object
        dopant : str
            Dopant element symbol
        dopant_fraction : float
            Fraction of host cations replaced
        n_vacancies : int
            Number of oxygen vacancies created
        parent_formula : str, optional
            Original ceria formula (e.g., "Ce32O64")
        """
        super().__init__(atoms)
        self.dopant = dopant
        self.dopant_fraction = dopant_fraction
        self.n_vacancies = n_vacancies
        self.parent_formula = parent_formula

    @classmethod
    def from_ceria(
        cls,
        ceria: Union[Atoms, BaseStructure],
        dopant: str,
        dopant_fraction: float,
        host_cation: str = "Ce",
        vacancy_placement: str = "random",
        dopant_placement: str = "random",
        dopant_preference: float = 0.7,
        vacancy_preference: float = 0.7,
        pr_trivalent_fraction: float = 1.0,
        z_threshold: float = 0.3,
        random_seed: Optional[int] = None,
    ) -> "DopedCeriaStructure":
        """
        Create doped ceria from pristine CeO2 structure.

        Parameters
        ----------
        ceria : Atoms or BaseStructure
            Parent CeO2 structure
        dopant : str
            Dopant element
        dopant_fraction : float
            Fraction of Ce to replace
        host_cation : str
            Cation to replace (default: "Ce")
        vacancy_placement : str
            Vacancy placement strategy
        dopant_placement : str
            Dopant placement strategy
        dopant_preference : float
            Dopant placement preference strength
        vacancy_preference : float
            Vacancy placement preference strength
        pr_trivalent_fraction : float
            For Pr: fraction that is Pr³⁺
        z_threshold : float
            Fraction of slab considered surface
        random_seed : int, optional
            Random seed

        Returns
        -------
        DopedCeriaStructure
            New doped ceria structure
        """
        if isinstance(ceria, BaseStructure):
            atoms = ceria.atoms.copy()
            parent_formula = ceria.formula
        else:
            atoms = ceria.copy()
            parent_formula = atoms.get_chemical_formula()

        # Count original O atoms for vacancy tracking
        original_o_count = sum(1 for s in atoms.get_chemical_symbols() if s == "O")

        builder = DopantBuilder(atoms)
        doped_atoms = builder.create_doped_structure(
            dopant=dopant,
            dopant_fraction=dopant_fraction,
            host_cation=host_cation,
            vacancy_placement=vacancy_placement,
            dopant_placement=dopant_placement,
            dopant_preference=dopant_preference,
            vacancy_preference=vacancy_preference,
            pr_trivalent_fraction=pr_trivalent_fraction,
            z_threshold=z_threshold,
            random_seed=random_seed,
        )

        # Count vacancies created
        new_o_count = sum(1 for s in doped_atoms.get_chemical_symbols() if s == "O")
        n_vacancies = original_o_count - new_o_count

        return cls(
            doped_atoms,
            dopant=dopant,
            dopant_fraction=dopant_fraction,
            n_vacancies=n_vacancies,
            parent_formula=parent_formula,
        )

    def get_stoichiometry(self) -> Dict[str, float]:
        """
        Calculate the stoichiometry relative to total cation sites.

        For Ce₁₋ₓGdₓO₂₋δ, returns ratios relative to cation sites.

        Returns
        -------
        dict
            Stoichiometric ratios
        """
        counts = {}
        for symbol in self.atoms.get_chemical_symbols():
            counts[symbol] = counts.get(symbol, 0) + 1

        # Sum cation counts (Ce + dopants)
        cation_elements = ["Ce", "Zr"] + list(ACCEPTOR_DOPANTS.keys())
        total_cations = sum(counts.get(c, 0) for c in cation_elements)

        if total_cations == 0:
            total_cations = 1  # Avoid division by zero

        stoich = {k: v / total_cations for k, v in counts.items()}

        return stoich

    def get_dopant_name(self) -> str:
        """
        Get the full name of the doped material.

        Returns
        -------
        str
            Material name (e.g., "GDC", "SDC", "PDC")
        """
        dopant_names = {
            "Gd": "GDC",  # Gadolinium-doped Ceria
            "Sm": "SDC",  # Samarium-doped Ceria
            "Pr": "PDC",  # Praseodymium-doped Ceria
            "Y": "YDC",   # Yttrium-doped Ceria
            "La": "LDC",  # Lanthanum-doped Ceria
            "Nd": "NDC",  # Neodymium-doped Ceria
        }
        return dopant_names.get(self.dopant, f"{self.dopant}-doped Ceria")

    def __repr__(self) -> str:
        return (
            f"DopedCeriaStructure({self.get_dopant_name()}, "
            f"{self.dopant}_frac={self.dopant_fraction:.1%}, "
            f"n_vac={self.n_vacancies})"
        )


def analyze_dopant_distribution(
    atoms: Atoms,
    dopant: str,
    reference_atoms: Optional[Atoms] = None,
    z_threshold: float = 0.3,
    near_dopant_cutoff: float = 3.5,
) -> Dict[str, Any]:
    """
    Analyze the distribution of dopants and vacancies in a structure.

    Parameters
    ----------
    atoms : Atoms
        Doped structure to analyze
    dopant : str
        Dopant element symbol (e.g., "Gd", "Sm")
    reference_atoms : Atoms, optional
        Original structure before doping (to count vacancies)
    z_threshold : float
        Fraction of slab height considered as surface region (default: 0.3 = top 30%)
    near_dopant_cutoff : float
        Distance cutoff for "near dopant" analysis in Angstroms (default: 3.5)

    Returns
    -------
    dict
        Analysis results with keys:
        - dopant_total: total number of dopant atoms
        - dopant_surface: dopants in surface region
        - dopant_surface_fraction: fraction of dopants in surface region
        - dopant_bulk: dopants in bulk region
        - ce_total: total Ce atoms remaining
        - o_total: total O atoms
        - vacancy_total: total vacancies (if reference provided)
        - vacancy_surface: vacancies in surface region
        - vacancy_near_dopant: vacancies within cutoff of dopant atoms
        - vacancy_near_dopant_fraction: fraction of vacancies near dopants

    Examples
    --------
    >>> stats = analyze_dopant_distribution(gdc, dopant="Gd", reference_atoms=ceo2)
    >>> print(f"Gd in surface: {stats['dopant_surface_fraction']:.1%}")
    >>> print(f"Vacancies near Gd: {stats['vacancy_near_dopant_fraction']:.1%}")
    """
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    z_positions = positions[:, 2]

    # Define surface region (top z_threshold fraction of slab)
    z_min, z_max = z_positions.min(), z_positions.max()
    z_cutoff = z_max - z_threshold * (z_max - z_min)

    # Get indices by element
    dopant_indices = [i for i, s in enumerate(symbols) if s == dopant]
    ce_indices = [i for i, s in enumerate(symbols) if s == "Ce"]
    o_indices = [i for i, s in enumerate(symbols) if s == "O"]

    # Count dopants in surface vs bulk
    dopant_surface = [i for i in dopant_indices if z_positions[i] >= z_cutoff]
    dopant_bulk = [i for i in dopant_indices if z_positions[i] < z_cutoff]

    results = {
        "dopant": dopant,
        "dopant_total": len(dopant_indices),
        "dopant_surface": len(dopant_surface),
        "dopant_bulk": len(dopant_bulk),
        "dopant_surface_fraction": len(dopant_surface) / len(dopant_indices) if dopant_indices else 0,
        "ce_total": len(ce_indices),
        "o_total": len(o_indices),
        "z_threshold": z_threshold,
        "z_cutoff": z_cutoff,
    }

    # Calculate dopant fraction
    total_cations = len(dopant_indices) + len(ce_indices)
    if total_cations > 0:
        results["dopant_fraction"] = len(dopant_indices) / total_cations
    else:
        results["dopant_fraction"] = 0

    # Vacancy analysis (requires reference structure)
    if reference_atoms is not None:
        ref_symbols = reference_atoms.get_chemical_symbols()
        ref_positions = reference_atoms.get_positions()

        # Find vacancy positions by comparing O sites
        ref_o_positions = []
        for i, s in enumerate(ref_symbols):
            if s == "O":
                ref_o_positions.append(ref_positions[i])

        current_o_positions = []
        for i, s in enumerate(symbols):
            if s == "O":
                current_o_positions.append(positions[i])

        # Find missing O positions (vacancies)
        vacancy_positions = []
        for ref_pos in ref_o_positions:
            found = False
            for cur_pos in current_o_positions:
                if np.linalg.norm(ref_pos - cur_pos) < 0.5:  # Within 0.5 A
                    found = True
                    break
            if not found:
                vacancy_positions.append(ref_pos)

        vacancy_positions = np.array(vacancy_positions) if vacancy_positions else np.array([]).reshape(0, 3)

        # Count vacancies by region
        if len(vacancy_positions) > 0:
            vacancy_z = vacancy_positions[:, 2]
            vacancy_surface = np.sum(vacancy_z >= z_cutoff)
            vacancy_bulk = np.sum(vacancy_z < z_cutoff)

            # Count vacancies near dopant atoms
            dopant_positions = positions[dopant_indices] if dopant_indices else np.array([]).reshape(0, 3)
            vacancy_near_dopant = 0
            for vac_pos in vacancy_positions:
                if len(dopant_positions) > 0:
                    min_dist = np.min(np.linalg.norm(dopant_positions - vac_pos, axis=1))
                    if min_dist <= near_dopant_cutoff:
                        vacancy_near_dopant += 1

            results.update({
                "vacancy_total": len(vacancy_positions),
                "vacancy_surface": int(vacancy_surface),
                "vacancy_bulk": int(vacancy_bulk),
                "vacancy_surface_fraction": vacancy_surface / len(vacancy_positions) if len(vacancy_positions) > 0 else 0,
                "vacancy_near_dopant": vacancy_near_dopant,
                "vacancy_near_dopant_fraction": vacancy_near_dopant / len(vacancy_positions) if len(vacancy_positions) > 0 else 0,
                "near_dopant_cutoff": near_dopant_cutoff,
            })
        else:
            results.update({
                "vacancy_total": 0,
                "vacancy_surface": 0,
                "vacancy_bulk": 0,
                "vacancy_surface_fraction": 0,
                "vacancy_near_dopant": 0,
                "vacancy_near_dopant_fraction": 0,
                "near_dopant_cutoff": near_dopant_cutoff,
            })

    return results


def generate_dopant_series(
    structure: Union[Atoms, BaseStructure],
    dopant: str,
    dopant_fractions: List[float],
    n_configs: int = 3,
    host_cation: str = "Ce",
    random_seed: Optional[int] = None,
) -> List[Dict[str, Any]]:
    """
    Generate structures with varying dopant concentrations.

    .. deprecated::
        Use ``DopantBuilder.create_doped_pool(dopant_fraction=[...])`` instead.
        This function is kept for backwards compatibility.

    Parameters
    ----------
    structure : Atoms or BaseStructure
        Base CeO2 structure
    dopant : str
        Dopant element
    dopant_fractions : list
        List of dopant fractions to generate (e.g., [0.05, 0.10, 0.15, 0.20])
    n_configs : int
        Number of random configurations per concentration (default: 3)
    host_cation : str
        Cation to replace (default: "Ce")
    random_seed : int, optional
        Base random seed

    Returns
    -------
    list of dict
        List of configurations with 'atoms', 'dopant', 'dopant_fraction', 'config_id'

    Examples
    --------
    >>> # Preferred: use create_doped_pool with list of fractions
    >>> builder = DopantBuilder(ceo2_slab)
    >>> pool = builder.create_doped_pool(
    ...     dopant="Gd",
    ...     dopant_fraction=[0.05, 0.10, 0.15, 0.20],
    ...     n_configs=5,
    ... )
    """
    import warnings
    warnings.warn(
        "generate_dopant_series() is deprecated. Use "
        "DopantBuilder.create_doped_pool(dopant_fraction=[...]) instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    builder = DopantBuilder(structure)
    return builder.create_doped_pool(
        dopant=dopant,
        dopant_fraction=dopant_fractions,
        n_configs=n_configs,
        strategies=["random"],
        host_cation=host_cation,
        random_seed=random_seed,
    )


def print_dopant_analysis(
    stats: Dict[str, Any],
    title: str = "Dopant Distribution Analysis",
) -> None:
    """
    Print a formatted summary of dopant distribution analysis.

    Parameters
    ----------
    stats : dict
        Results from analyze_dopant_distribution()
    title : str
        Title for the report

    Examples
    --------
    >>> stats = analyze_dopant_distribution(gdc, "Gd", reference_atoms=ceo2)
    >>> print_dopant_analysis(stats)
    """
    print("=" * 60)
    print(title)
    print("=" * 60)

    dopant = stats.get("dopant", "?")
    print(f"\nDopant: {dopant}")
    print(f"  Total {dopant} atoms: {stats['dopant_total']}")
    print(f"  {dopant} in surface (top {stats['z_threshold']*100:.0f}%): {stats['dopant_surface']} ({stats['dopant_surface_fraction']*100:.1f}%)")
    print(f"  {dopant} in bulk: {stats['dopant_bulk']}")

    if "dopant_fraction" in stats:
        print(f"  Dopant fraction: {stats['dopant_fraction']*100:.1f}%")

    print(f"\nHost Cations:")
    print(f"  Ce atoms remaining: {stats['ce_total']}")
    print(f"  O atoms: {stats['o_total']}")

    if "vacancy_total" in stats:
        print(f"\nVacancy Distribution:")
        print(f"  Total vacancies: {stats['vacancy_total']}")
        print(f"  Vacancies in surface: {stats['vacancy_surface']} ({stats['vacancy_surface_fraction']*100:.1f}%)")
        print(f"  Vacancies in bulk: {stats['vacancy_bulk']}")
        print(f"  Vacancies near {dopant} (< {stats['near_dopant_cutoff']} A): {stats['vacancy_near_dopant']} ({stats['vacancy_near_dopant_fraction']*100:.1f}%)")

    # Print charge balance check
    if stats['dopant_total'] > 0 and "vacancy_total" in stats:
        expected_vacancies = stats['dopant_total'] // 2
        actual_vacancies = stats['vacancy_total']
        print(f"\nCharge Balance Check:")
        print(f"  Expected vacancies (n_dopant/2): {expected_vacancies}")
        print(f"  Actual vacancies: {actual_vacancies}")
        if actual_vacancies == expected_vacancies:
            print("  Status: Charge balanced")
        else:
            print(f"  Status: Imbalanced by {abs(actual_vacancies - expected_vacancies)} vacancies")

    print("=" * 60)
