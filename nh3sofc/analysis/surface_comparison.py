"""Surface comparison and ranking for catalyst screening.

Provides tools for comparing catalytic surfaces based on
energetics, kinetics, and activity descriptors.
"""

from typing import Optional, List, Dict, Any, Tuple
import numpy as np

from .rds import EnergySpanModel, BEPRelation


class SurfaceComparator:
    """
    Compare and rank catalyst surfaces.

    Supports comparison based on:
    - Energy profiles
    - Activation barriers
    - Activity descriptors
    - Volcano relationships

    Examples
    --------
    >>> surfaces = {
    ...     'LaO-term': {'NH3*': 0.0, 'NH2*+H*': 0.8, 'NH*+2H*': 1.5, 'N*+3H*': 1.2},
    ...     'VO2-term': {'NH3*': 0.0, 'NH2*+H*': 0.6, 'NH*+2H*': 1.2, 'N*+3H*': 0.9},
    ... }
    >>> comp = SurfaceComparator(surfaces)
    >>> ranking = comp.rank_surfaces()
    """

    def __init__(
        self,
        surfaces: Dict[str, Dict[str, float]],
        barriers: Optional[Dict[str, Dict[str, float]]] = None,
    ):
        """
        Initialize SurfaceComparator.

        Parameters
        ----------
        surfaces : dict
            Surface name -> {intermediate: energy}
        barriers : dict, optional
            Surface name -> {step: barrier}
        """
        self.surfaces = surfaces
        self.barriers = barriers or {}
        self.bep = BEPRelation()

    def get_energy_span(self, surface: str) -> float:
        """
        Get energy span for a surface.

        Parameters
        ----------
        surface : str
            Surface name

        Returns
        -------
        float
            Energy span in eV
        """
        if surface not in self.surfaces:
            raise ValueError(f"Unknown surface: {surface}")

        pathway = self.surfaces[surface]

        # Get transition state energies (use barriers if available, else BEP)
        ts_energies = {}
        steps = list(pathway.keys())

        for i in range(len(steps) - 1):
            initial, final = steps[i], steps[i+1]
            ts_name = f"TS_{initial}_{final}"
            dE = pathway[final] - pathway[initial]

            if surface in self.barriers and f"{initial}→{final}" in self.barriers[surface]:
                barrier = self.barriers[surface][f"{initial}→{final}"]
            else:
                barrier = self.bep.estimate_barrier(dE)

            ts_energies[ts_name] = pathway[initial] + barrier

        # Calculate energy span
        model = EnergySpanModel(pathway, ts_energies)
        span, _, _ = model.get_energy_span()

        return span

    def get_max_barrier(self, surface: str) -> Tuple[str, float]:
        """
        Get maximum barrier for a surface.

        Parameters
        ----------
        surface : str
            Surface name

        Returns
        -------
        tuple
            (step_name, barrier)
        """
        pathway = self.surfaces[surface]
        steps = list(pathway.keys())

        max_barrier = 0.0
        max_step = ""

        for i in range(len(steps) - 1):
            initial, final = steps[i], steps[i+1]
            step_name = f"{initial}→{final}"
            dE = pathway[final] - pathway[initial]

            if surface in self.barriers and step_name in self.barriers[surface]:
                barrier = self.barriers[surface][step_name]
            else:
                barrier = self.bep.estimate_barrier(dE)

            if barrier > max_barrier:
                max_barrier = barrier
                max_step = step_name

        return max_step, max_barrier

    def get_overall_reaction_energy(self, surface: str) -> float:
        """
        Get overall reaction energy.

        Parameters
        ----------
        surface : str
            Surface name

        Returns
        -------
        float
            Overall reaction energy in eV
        """
        pathway = self.surfaces[surface]
        steps = list(pathway.keys())

        if len(steps) < 2:
            return 0.0

        return pathway[steps[-1]] - pathway[steps[0]]

    def rank_by_energy_span(self) -> List[Tuple[str, float]]:
        """
        Rank surfaces by energy span (lower is better).

        Returns
        -------
        list
            List of (surface_name, energy_span) sorted by span
        """
        spans = [(s, self.get_energy_span(s)) for s in self.surfaces]
        return sorted(spans, key=lambda x: x[1])

    def rank_by_max_barrier(self) -> List[Tuple[str, float]]:
        """
        Rank surfaces by maximum barrier (lower is better).

        Returns
        -------
        list
            List of (surface_name, max_barrier) sorted by barrier
        """
        barriers = [(s, self.get_max_barrier(s)[1]) for s in self.surfaces]
        return sorted(barriers, key=lambda x: x[1])

    def rank_surfaces(
        self,
        method: str = "energy_span",
    ) -> List[Tuple[str, float]]:
        """
        Rank surfaces by specified method.

        Parameters
        ----------
        method : str
            Ranking method: "energy_span", "max_barrier", "reaction_energy"

        Returns
        -------
        list
            Ranked surfaces
        """
        if method == "energy_span":
            return self.rank_by_energy_span()
        elif method == "max_barrier":
            return self.rank_by_max_barrier()
        elif method == "reaction_energy":
            energies = [(s, self.get_overall_reaction_energy(s)) for s in self.surfaces]
            return sorted(energies, key=lambda x: x[1])
        else:
            raise ValueError(f"Unknown method: {method}")

    def get_best_surface(self, method: str = "energy_span") -> str:
        """
        Get the best performing surface.

        Parameters
        ----------
        method : str
            Ranking method

        Returns
        -------
        str
            Name of best surface
        """
        ranking = self.rank_surfaces(method)
        return ranking[0][0] if ranking else ""

    def compare_surfaces(
        self,
        surface1: str,
        surface2: str,
    ) -> Dict[str, Any]:
        """
        Compare two surfaces in detail.

        Parameters
        ----------
        surface1 : str
            First surface name
        surface2 : str
            Second surface name

        Returns
        -------
        dict
            Comparison results
        """
        span1 = self.get_energy_span(surface1)
        span2 = self.get_energy_span(surface2)

        step1, barrier1 = self.get_max_barrier(surface1)
        step2, barrier2 = self.get_max_barrier(surface2)

        dG1 = self.get_overall_reaction_energy(surface1)
        dG2 = self.get_overall_reaction_energy(surface2)

        return {
            surface1: {
                "energy_span": span1,
                "max_barrier": barrier1,
                "rds_step": step1,
                "reaction_energy": dG1,
            },
            surface2: {
                "energy_span": span2,
                "max_barrier": barrier2,
                "rds_step": step2,
                "reaction_energy": dG2,
            },
            "better_surface": surface1 if span1 < span2 else surface2,
            "span_difference": abs(span1 - span2),
        }

    def get_descriptor_values(
        self,
        descriptor_step: str = "NH3*",
    ) -> Dict[str, float]:
        """
        Get descriptor values (e.g., NH3 binding energy).

        Parameters
        ----------
        descriptor_step : str
            Step to use as descriptor

        Returns
        -------
        dict
            Surface -> descriptor value
        """
        descriptors = {}
        for surface, pathway in self.surfaces.items():
            if descriptor_step in pathway:
                descriptors[surface] = pathway[descriptor_step]
        return descriptors

    def plot_energy_profiles(
        self,
        filename: Optional[str] = None,
        surfaces: Optional[List[str]] = None,
    ) -> None:
        """
        Plot energy profiles for surfaces.

        Parameters
        ----------
        filename : str, optional
            Output file path
        surfaces : list, optional
            Surfaces to include (default: all)
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not installed. Cannot plot.")
            return

        surfaces_to_plot = surfaces or list(self.surfaces.keys())

        fig, ax = plt.subplots(figsize=(10, 6))

        for surface in surfaces_to_plot:
            pathway = self.surfaces[surface]
            steps = list(pathway.keys())
            energies = [pathway[s] for s in steps]

            x = np.arange(len(steps))
            ax.plot(x, energies, 'o-', label=surface, linewidth=2, markersize=8)

        ax.set_xticks(np.arange(len(steps)))
        ax.set_xticklabels(steps, rotation=45, ha='right')
        ax.set_ylabel("Energy (eV)", fontsize=12)
        ax.set_title("Energy Profiles Comparison", fontsize=14)
        ax.legend()
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

        plt.tight_layout()

        if filename:
            plt.savefig(filename, dpi=150)
            print(f"Saved to {filename}")
        else:
            plt.show()

    def plot_volcano(
        self,
        descriptor_step: str = "NH3*",
        activity_metric: str = "energy_span",
        filename: Optional[str] = None,
    ) -> None:
        """
        Plot volcano curve.

        Parameters
        ----------
        descriptor_step : str
            Descriptor (x-axis)
        activity_metric : str
            Activity measure (y-axis)
        filename : str, optional
            Output file path
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not installed. Cannot plot.")
            return

        descriptors = self.get_descriptor_values(descriptor_step)

        activities = {}
        for surface in descriptors:
            if activity_metric == "energy_span":
                activities[surface] = -self.get_energy_span(surface)  # Negative: lower span = higher activity
            elif activity_metric == "max_barrier":
                _, barrier = self.get_max_barrier(surface)
                activities[surface] = -barrier

        # Plot
        fig, ax = plt.subplots(figsize=(8, 6))

        x = [descriptors[s] for s in descriptors]
        y = [activities[s] for s in descriptors]
        labels = list(descriptors.keys())

        ax.scatter(x, y, s=100)

        for i, label in enumerate(labels):
            ax.annotate(label, (x[i], y[i]), textcoords="offset points",
                       xytext=(5, 5), fontsize=9)

        ax.set_xlabel(f"{descriptor_step} Energy (eV)", fontsize=12)
        ax.set_ylabel(f"Activity (-{activity_metric})", fontsize=12)
        ax.set_title("Volcano Plot", fontsize=14)

        plt.tight_layout()

        if filename:
            plt.savefig(filename, dpi=150)
            print(f"Saved to {filename}")
        else:
            plt.show()

    def print_summary(self) -> None:
        """Print comparison summary."""
        print("\n" + "=" * 70)
        print("Surface Comparison Summary")
        print("=" * 70)

        # Energy spans
        print("\nRanking by Energy Span:")
        for i, (surface, span) in enumerate(self.rank_by_energy_span(), 1):
            print(f"  {i}. {surface:20s}: δE = {span:.3f} eV")

        # Max barriers
        print("\nRate-Determining Steps:")
        for surface in self.surfaces:
            step, barrier = self.get_max_barrier(surface)
            print(f"  {surface:20s}: {step} (E_a = {barrier:.3f} eV)")

        # Overall thermodynamics
        print("\nOverall Reaction Energies:")
        for surface in self.surfaces:
            dG = self.get_overall_reaction_energy(surface)
            sign = "exo" if dG < 0 else "endo"
            print(f"  {surface:20s}: ΔG = {dG:+.3f} eV ({sign}thermic)")

        print("=" * 70)


class ActivityDescriptor:
    """
    Activity descriptor analysis for catalyst screening.

    Uses scaling relations to predict activity from simple descriptors.
    """

    def __init__(self):
        """Initialize ActivityDescriptor."""
        self.data = {}  # descriptor_value -> activity

    def add_point(
        self,
        surface: str,
        descriptor: float,
        activity: float,
    ) -> None:
        """
        Add a data point.

        Parameters
        ----------
        surface : str
            Surface identifier
        descriptor : float
            Descriptor value
        activity : float
            Activity measure
        """
        self.data[surface] = {
            "descriptor": descriptor,
            "activity": activity,
        }

    def fit_scaling_relation(self) -> Tuple[float, float, float]:
        """
        Fit linear scaling relation.

        activity = a * descriptor + b

        Returns
        -------
        tuple
            (slope, intercept, r_squared)
        """
        if len(self.data) < 2:
            return (0.0, 0.0, 0.0)

        x = np.array([v["descriptor"] for v in self.data.values()])
        y = np.array([v["activity"] for v in self.data.values()])

        coeffs = np.polyfit(x, y, 1)
        slope, intercept = coeffs

        # R-squared
        y_pred = slope * x + intercept
        ss_res = np.sum((y - y_pred) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        return slope, intercept, r_sq

    def predict_activity(self, descriptor: float) -> float:
        """
        Predict activity from descriptor.

        Parameters
        ----------
        descriptor : float
            Descriptor value

        Returns
        -------
        float
            Predicted activity
        """
        slope, intercept, _ = self.fit_scaling_relation()
        return slope * descriptor + intercept

    def find_optimal_descriptor(
        self,
        descriptor_range: Tuple[float, float] = (-2.0, 0.0),
    ) -> float:
        """
        Find optimal descriptor value.

        For volcano-type relations, finds the peak.

        Parameters
        ----------
        descriptor_range : tuple
            Range to search

        Returns
        -------
        float
            Optimal descriptor value
        """
        # Simple: return mean of best performers
        if len(self.data) < 2:
            return 0.0

        sorted_data = sorted(self.data.items(),
                            key=lambda x: x[1]["activity"],
                            reverse=True)

        # Take top 20% or at least 1
        n_top = max(1, len(sorted_data) // 5)
        top_descriptors = [v[1]["descriptor"] for v in sorted_data[:n_top]]

        return np.mean(top_descriptors)


def compare_surfaces(
    surfaces: Dict[str, Dict[str, float]],
    method: str = "energy_span",
) -> List[Tuple[str, float]]:
    """
    Compare and rank surfaces.

    Parameters
    ----------
    surfaces : dict
        Surface energetics data
    method : str
        Comparison method

    Returns
    -------
    list
        Ranked surfaces
    """
    comparator = SurfaceComparator(surfaces)
    return comparator.rank_surfaces(method)


def get_best_catalyst(
    surfaces: Dict[str, Dict[str, float]],
    method: str = "energy_span",
) -> str:
    """
    Get the best catalyst surface.

    Parameters
    ----------
    surfaces : dict
        Surface energetics data
    method : str
        Ranking method

    Returns
    -------
    str
        Name of best surface
    """
    comparator = SurfaceComparator(surfaces)
    return comparator.get_best_surface(method)
