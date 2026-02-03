"""Energetics analysis for adsorption and surface calculations.

Provides tools for calculating adsorption energies, surface energies,
and related thermodynamic quantities.
"""

from pathlib import Path
from typing import Optional, List, Union, Dict, Any, Tuple
import numpy as np
from ase import Atoms
from ase.io import read as ase_read

from ..core.constants import GAS_PHASE_ENERGIES


class AdsorptionEnergyCalculator:
    """
    Calculate adsorption energies for adsorbate-surface systems.

    E_ads = E(adsorbate/surface) - E(surface) - E(adsorbate_gas)

    Examples
    --------
    >>> calc = AdsorptionEnergyCalculator()
    >>> calc.set_surface_energy(-150.0)
    >>> calc.set_gas_reference("NH3", -19.54)
    >>> E_ads = calc.calculate(-170.0)  # E_ads = -170 - (-150) - (-19.54) = -0.46 eV
    """

    def __init__(
        self,
        e_surface: Optional[float] = None,
        gas_references: Optional[Dict[str, float]] = None,
    ):
        """
        Initialize AdsorptionEnergyCalculator.

        Parameters
        ----------
        e_surface : float, optional
            Clean surface energy in eV
        gas_references : dict, optional
            Gas phase reference energies
        """
        self.e_surface = e_surface
        self.gas_references = gas_references or GAS_PHASE_ENERGIES.copy()

    def set_surface_energy(self, energy: float) -> None:
        """Set clean surface energy."""
        self.e_surface = energy

    def set_gas_reference(self, molecule: str, energy: float) -> None:
        """Set gas phase reference energy for a molecule."""
        self.gas_references[molecule] = energy

    def calculate(
        self,
        e_total: float,
        adsorbate: str = "NH3",
        n_adsorbates: int = 1,
    ) -> float:
        """
        Calculate adsorption energy.

        E_ads = E(ads/surf) - E(surf) - n * E(gas)

        Parameters
        ----------
        e_total : float
            Total energy of adsorbate/surface system
        adsorbate : str
            Adsorbate molecule name
        n_adsorbates : int
            Number of adsorbate molecules

        Returns
        -------
        float
            Adsorption energy in eV
        """
        if self.e_surface is None:
            raise ValueError("Surface energy not set")

        if adsorbate not in self.gas_references:
            raise ValueError(f"No gas reference for {adsorbate}")

        e_gas = self.gas_references[adsorbate]
        e_ads = e_total - self.e_surface - n_adsorbates * e_gas

        return e_ads

    def calculate_per_adsorbate(
        self,
        e_total: float,
        adsorbate: str = "NH3",
        n_adsorbates: int = 1,
    ) -> float:
        """Calculate adsorption energy per adsorbate molecule."""
        e_ads_total = self.calculate(e_total, adsorbate, n_adsorbates)
        return e_ads_total / n_adsorbates

    def calculate_differential(
        self,
        e_n: float,
        e_n_minus_1: float,
        adsorbate: str = "NH3",
    ) -> float:
        """
        Calculate differential adsorption energy.

        dE_ads = E(n ads) - E((n-1) ads) - E(gas)

        Parameters
        ----------
        e_n : float
            Energy with n adsorbates
        e_n_minus_1 : float
            Energy with n-1 adsorbates
        adsorbate : str
            Adsorbate molecule

        Returns
        -------
        float
            Differential adsorption energy
        """
        e_gas = self.gas_references[adsorbate]
        return e_n - e_n_minus_1 - e_gas


class SurfaceEnergyCalculator:
    """
    Calculate surface formation energies.

    γ = (E_slab - n * E_bulk) / (2 * A)

    Examples
    --------
    >>> calc = SurfaceEnergyCalculator()
    >>> gamma = calc.calculate(
    ...     e_slab=-200.0,
    ...     e_bulk=-5.0,
    ...     n_atoms=40,
    ...     area=50.0,
    ... )
    """

    def __init__(self, e_bulk_per_atom: Optional[float] = None):
        """
        Initialize SurfaceEnergyCalculator.

        Parameters
        ----------
        e_bulk_per_atom : float, optional
            Bulk energy per atom
        """
        self.e_bulk_per_atom = e_bulk_per_atom

    def set_bulk_energy(self, energy: float, per_atom: bool = True) -> None:
        """Set bulk reference energy."""
        self.e_bulk_per_atom = energy if per_atom else None

    def calculate(
        self,
        e_slab: float,
        n_atoms: int,
        area: float,
        e_bulk_per_atom: Optional[float] = None,
        symmetric: bool = True,
    ) -> float:
        """
        Calculate surface energy.

        Parameters
        ----------
        e_slab : float
            Slab total energy in eV
        n_atoms : int
            Number of atoms in slab
        area : float
            Surface area in A^2
        e_bulk_per_atom : float, optional
            Bulk energy per atom
        symmetric : bool
            Whether slab has two equivalent surfaces

        Returns
        -------
        float
            Surface energy in eV/A^2
        """
        e_bulk = e_bulk_per_atom or self.e_bulk_per_atom
        if e_bulk is None:
            raise ValueError("Bulk energy not set")

        n_surfaces = 2 if symmetric else 1
        gamma = (e_slab - n_atoms * e_bulk) / (n_surfaces * area)

        return gamma

    def calculate_j_m2(self, *args, **kwargs) -> float:
        """Calculate surface energy in J/m^2."""
        gamma_ev_a2 = self.calculate(*args, **kwargs)
        # 1 eV/A^2 = 16.0218 J/m^2
        return gamma_ev_a2 * 16.0218


def plot_surface_stability(
    surface_energies: Dict[str, float],
    title: str = "Surface Stability",
    filename: Optional[str] = None,
    unit: str = "J/m²",
    figsize: Tuple[float, float] = (8, 5),
) -> Tuple:
    """
    Plot surface energies as a bar chart to compare stability.

    Lower surface energy = more thermodynamically stable.

    Parameters
    ----------
    surface_energies : dict
        Surface energies by Miller index, e.g., {"111": 0.8, "110": 1.2}
    title : str
        Plot title
    filename : str, optional
        Save figure to file if provided
    unit : str
        Energy unit for y-axis label
    figsize : tuple
        Figure size (width, height)

    Returns
    -------
    tuple
        (fig, ax) matplotlib objects

    Examples
    --------
    >>> energies = {"111": 0.85, "110": 1.25, "100": 2.05}
    >>> fig, ax = plot_surface_stability(energies, title="CeO2 Surface Stability")
    """
    import matplotlib.pyplot as plt

    # Sort by surface energy (most stable first)
    sorted_data = sorted(surface_energies.items(), key=lambda x: x[1])
    millers = [f"({m})" for m, _ in sorted_data]
    gammas = [g for _, g in sorted_data]

    # Color gradient: green (stable) to red (unstable)
    colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(millers)))

    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar(millers, gammas, color=colors, edgecolor='black', linewidth=1.2)

    # Add value labels on bars
    for bar, gamma in zip(bars, gammas):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height + 0.05 * max(gammas),
                f'{gamma:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax.set_xlabel('Surface', fontsize=12)
    ax.set_ylabel(f'Surface Energy ({unit})', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')

    # Add stability annotation
    ax.annotate('← More stable', xy=(0, 0), xytext=(0.02, 0.98),
                xycoords='axes fraction', fontsize=10, color='green',
                ha='left', va='top')

    ax.set_ylim(0, max(gammas) * 1.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"Saved: {filename}")

    return fig, ax


class DecompositionEnergetics:
    """
    Analyze energetics of NH3 decomposition pathway.

    Tracks energies for: NH3* → NH2*+H* → NH*+2H* → N*+3H*
    """

    STEPS = ["NH3*", "NH2*+H*", "NH*+2H*", "N*+3H*"]

    def __init__(
        self,
        e_surface: Optional[float] = None,
    ):
        """
        Initialize DecompositionEnergetics.

        Parameters
        ----------
        e_surface : float, optional
            Clean surface energy
        """
        self.e_surface = e_surface
        self.energies = {}
        self.structures = {}

    def set_surface_energy(self, energy: float) -> None:
        """Set clean surface energy."""
        self.e_surface = energy

    def add_step_energy(
        self,
        step: str,
        energy: float,
        structure: Optional[Atoms] = None,
    ) -> None:
        """
        Add energy for a decomposition step.

        Parameters
        ----------
        step : str
            Step name (e.g., "NH3*", "NH2*+H*")
        energy : float
            Total energy in eV
        structure : Atoms, optional
            Optimized structure
        """
        self.energies[step] = energy
        if structure is not None:
            self.structures[step] = structure

    def get_relative_energies(
        self,
        reference: str = "NH3*",
    ) -> Dict[str, float]:
        """
        Get energies relative to reference state.

        Parameters
        ----------
        reference : str
            Reference state for energy zero

        Returns
        -------
        dict
            Relative energies for each step
        """
        if reference not in self.energies:
            raise ValueError(f"Reference {reference} not found")

        e_ref = self.energies[reference]
        return {step: e - e_ref for step, e in self.energies.items()}

    def get_reaction_energies(self) -> Dict[str, float]:
        """
        Get step-by-step reaction energies.

        Returns
        -------
        dict
            Reaction energy for each elementary step
        """
        relative = self.get_relative_energies()

        reactions = {}
        step_pairs = [
            ("NH3*→NH2*+H*", "NH3*", "NH2*+H*"),
            ("NH2*+H*→NH*+2H*", "NH2*+H*", "NH*+2H*"),
            ("NH*+2H*→N*+3H*", "NH*+2H*", "N*+3H*"),
        ]

        for name, initial, final in step_pairs:
            if initial in relative and final in relative:
                reactions[name] = relative[final] - relative[initial]

        return reactions

    def get_energy_span(self) -> Tuple[float, str, str]:
        """
        Calculate energy span (approximate activation energy).

        Returns highest energy difference in the pathway.

        Returns
        -------
        tuple
            (energy_span, summit_state, trough_state)
        """
        relative = self.get_relative_energies()

        # Find global maximum and minimum
        max_state = max(relative, key=relative.get)
        min_state = min(relative, key=relative.get)

        span = relative[max_state] - relative[min_state]

        return span, max_state, min_state

    def is_thermodynamically_favorable(self) -> bool:
        """Check if overall decomposition is exothermic."""
        relative = self.get_relative_energies()
        final_states = ["N*+3H*"]

        for state in final_states:
            if state in relative and relative[state] < 0:
                return True
        return False

    def print_summary(self) -> None:
        """Print summary of decomposition energetics."""
        print("\n" + "=" * 60)
        print("NH3 Decomposition Energetics")
        print("=" * 60)

        # Relative energies
        relative = self.get_relative_energies()
        print("\nRelative Energies (vs NH3*):")
        for step in self.STEPS:
            if step in relative:
                print(f"  {step:15s}: {relative[step]:+.3f} eV")

        # Reaction energies
        reactions = self.get_reaction_energies()
        print("\nReaction Energies:")
        for name, dE in reactions.items():
            sign = "exo" if dE < 0 else "endo"
            print(f"  {name:25s}: {dE:+.3f} eV ({sign}thermic)")

        # Energy span
        span, summit, trough = self.get_energy_span()
        print(f"\nEnergy Span: {span:.3f} eV")
        print(f"  Summit: {summit}")
        print(f"  Trough: {trough}")

        # Overall favorability
        favorable = self.is_thermodynamically_favorable()
        print(f"\nThermodynamically favorable: {favorable}")

        print("=" * 60)


class BindingEnergyAnalyzer:
    """
    Analyze binding energies and trends.

    Useful for comparing adsorption across different surfaces
    or different adsorbates.
    """

    def __init__(self):
        """Initialize BindingEnergyAnalyzer."""
        self.data = {}

    def add_system(
        self,
        name: str,
        e_binding: float,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Add a system's binding energy.

        Parameters
        ----------
        name : str
            System identifier
        e_binding : float
            Binding energy in eV
        metadata : dict, optional
            Additional information
        """
        self.data[name] = {
            "e_binding": e_binding,
            "metadata": metadata or {},
        }

    def get_ranking(self, ascending: bool = True) -> List[Tuple[str, float]]:
        """
        Get systems ranked by binding energy.

        Parameters
        ----------
        ascending : bool
            If True, strongest binding first (most negative)

        Returns
        -------
        list
            List of (name, binding_energy) tuples
        """
        items = [(k, v["e_binding"]) for k, v in self.data.items()]
        return sorted(items, key=lambda x: x[1], reverse=not ascending)

    def get_statistics(self) -> Dict[str, float]:
        """
        Get statistical summary of binding energies.

        Returns
        -------
        dict
            Statistics (mean, std, min, max)
        """
        energies = [v["e_binding"] for v in self.data.values()]

        if not energies:
            return {}

        return {
            "mean": np.mean(energies),
            "std": np.std(energies),
            "min": np.min(energies),
            "max": np.max(energies),
            "range": np.max(energies) - np.min(energies),
        }

    def find_optimal(
        self,
        target: float,
        tolerance: float = 0.1,
    ) -> List[str]:
        """
        Find systems with binding energy near target.

        Parameters
        ----------
        target : float
            Target binding energy
        tolerance : float
            Acceptable deviation from target

        Returns
        -------
        list
            Names of systems within tolerance
        """
        optimal = []
        for name, info in self.data.items():
            if abs(info["e_binding"] - target) <= tolerance:
                optimal.append(name)
        return optimal


def calculate_adsorption_energy(
    e_adsorbate_surface: float,
    e_surface: float,
    e_adsorbate_gas: float,
) -> float:
    """
    Calculate adsorption energy.

    E_ads = E(ads/surf) - E(surf) - E(gas)

    Parameters
    ----------
    e_adsorbate_surface : float
        Energy of adsorbate on surface
    e_surface : float
        Clean surface energy
    e_adsorbate_gas : float
        Gas phase adsorbate energy

    Returns
    -------
    float
        Adsorption energy in eV
    """
    return e_adsorbate_surface - e_surface - e_adsorbate_gas


def calculate_coverage_dependent_energy(
    energies: List[float],
    coverages: List[float],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate coverage-dependent adsorption energies.

    Parameters
    ----------
    energies : list
        Adsorption energies at different coverages
    coverages : list
        Coverage values (ML)

    Returns
    -------
    tuple
        (coverages_array, differential_energies)
    """
    cov = np.array(coverages)
    e = np.array(energies)

    # Sort by coverage
    idx = np.argsort(cov)
    cov = cov[idx]
    e = e[idx]

    # Differential energies
    de = np.gradient(e, cov)

    return cov, de
