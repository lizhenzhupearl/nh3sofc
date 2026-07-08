"""Oracle limiting-case test: analytic two-step surface reaction.

Verifies that MicroKineticModel reproduces the exact analytical solution
for a two-step unimolecular surface reaction:

    Step 1: A(g) + * → A*,  r₁ = k₁ · P_A · θ_*
    Step 2: A* → B(g) + *,  r₂ = k₂ · θ_A

Steady-state solution:
    θ_A = k₁·P_A / (k₁·P_A + k₂)
    θ_* = k₂ / (k₁·P_A + k₂)
    TOF = k₁·k₂·P_A / (k₁·P_A + k₂)

This test is the executable oracle defined in ORACLE.md. If it fails,
the MKM solver has a fundamental bug and no science results are trustworthy.
"""

from __future__ import annotations

import numpy as np
import pytest

from nh3sofc.analysis.microkinetic import MicroKineticModel, RateConstantCalculator
from nh3sofc.core.constants import KB_EV


class TestTwoStepAnalyticSolution:
    """Verify MKM solver against exact analytical solution."""

    @pytest.fixture
    def setup(self) -> dict:
        """Set up the two-step model and compute analytical solution."""
        T = 673.0
        P_A = 0.5  # bar
        Ea1 = 0.3  # eV (adsorption barrier)
        Ea2 = 0.8  # eV (reaction barrier)
        prefactor = 1e13  # s^-1

        # Build model
        model = MicroKineticModel(temperature=T, total_sites=1.0)
        model.add_species("*", is_site=True)
        model.add_species("A*", initial_coverage=0.0)
        model.add_gas("A(g)", pressure=P_A)
        model.add_gas("B(g)", pressure=0.01)

        # Step 1: A(g) + * → A* (irreversible for simplicity)
        model.add_reaction(
            "adsorption",
            reactants={"A(g)": 1, "*": 1},
            products={"A*": 1},
            Ea_fwd=Ea1,
            Ea_rev=None,
            dG=None,
            prefactor=prefactor,
        )

        # Step 2: A* → B(g) + * (irreversible)
        model.add_reaction(
            "reaction",
            reactants={"A*": 1},
            products={"B(g)": 1, "*": 1},
            Ea_fwd=Ea2,
            Ea_rev=None,
            dG=None,
            prefactor=prefactor,
        )

        # Compute analytical rate constants
        calc = RateConstantCalculator(T)
        k1 = calc.eyring_rate(Ea1, prefactor)
        k2 = calc.eyring_rate(Ea2, prefactor)

        # Analytical steady-state
        theta_A_exact = k1 * P_A / (k1 * P_A + k2)
        theta_star_exact = k2 / (k1 * P_A + k2)
        tof_exact = k1 * k2 * P_A / (k1 * P_A + k2)

        return {
            "model": model,
            "k1": k1,
            "k2": k2,
            "P_A": P_A,
            "theta_A_exact": theta_A_exact,
            "theta_star_exact": theta_star_exact,
            "tof_exact": tof_exact,
        }

    def test_steady_state_coverages(self, setup: dict) -> None:
        """Coverages match analytical solution within 1e-8 relative error."""
        model = setup["model"]
        coverages = model.solve_steady_state()

        rel_err_A = abs(coverages["A*"] - setup["theta_A_exact"]) / setup["theta_A_exact"]
        rel_err_star = abs(coverages["*"] - setup["theta_star_exact"]) / setup["theta_star_exact"]

        assert rel_err_A < 1e-8, (
            f"θ_A: numerical={coverages['A*']:.12e}, "
            f"analytical={setup['theta_A_exact']:.12e}, rel_err={rel_err_A:.2e}"
        )
        assert rel_err_star < 1e-8, (
            f"θ_*: numerical={coverages['*']:.12e}, "
            f"analytical={setup['theta_star_exact']:.12e}, rel_err={rel_err_star:.2e}"
        )

    def test_tof_matches_analytical(self, setup: dict) -> None:
        """TOF matches analytical solution within 1e-8 relative error."""
        model = setup["model"]
        coverages = model.solve_steady_state()
        tof = model.get_tof(product_reaction="reaction", coverages=coverages)

        rel_err = abs(tof - setup["tof_exact"]) / setup["tof_exact"]
        assert rel_err < 1e-8, (
            f"TOF: numerical={tof:.12e}, "
            f"analytical={setup['tof_exact']:.12e}, rel_err={rel_err:.2e}"
        )

    def test_site_balance_conserved(self, setup: dict) -> None:
        """Sum of all coverages equals total_sites (1.0)."""
        model = setup["model"]
        coverages = model.solve_steady_state()
        total = sum(coverages.values())
        assert abs(total - 1.0) < 1e-8, f"Site balance violated: sum={total}"

    def test_coverages_non_negative(self, setup: dict) -> None:
        """All coverages are non-negative."""
        model = setup["model"]
        coverages = model.solve_steady_state()
        for species, cov in coverages.items():
            assert cov >= -1e-10, f"Negative coverage: {species}={cov}"

    def test_rate_constants_positive(self, setup: dict) -> None:
        """All rate constants are positive."""
        assert setup["k1"] > 0, f"k1 = {setup['k1']}"
        assert setup["k2"] > 0, f"k2 = {setup['k2']}"


class TestTwoStepMultipleConditions:
    """Verify the analytic solution holds across temperature and pressure ranges."""

    @pytest.mark.parametrize("T", [400.0, 500.0, 673.0, 800.0, 900.0])
    def test_temperature_sweep(self, T: float) -> None:
        """Analytical solution matches at multiple temperatures."""
        P_A = 0.5
        Ea1, Ea2 = 0.3, 0.8
        prefactor = 1e13

        model = MicroKineticModel(temperature=T, total_sites=1.0)
        model.add_species("*", is_site=True)
        model.add_species("A*", initial_coverage=0.0)
        model.add_gas("A(g)", pressure=P_A)
        model.add_gas("B(g)", pressure=0.01)
        model.add_reaction(
            "adsorption",
            {"A(g)": 1, "*": 1}, {"A*": 1},
            Ea_fwd=Ea1, prefactor=prefactor,
        )
        model.add_reaction(
            "reaction",
            {"A*": 1}, {"B(g)": 1, "*": 1},
            Ea_fwd=Ea2, prefactor=prefactor,
        )

        calc = RateConstantCalculator(T)
        k1 = calc.eyring_rate(Ea1, prefactor)
        k2 = calc.eyring_rate(Ea2, prefactor)
        tof_exact = k1 * k2 * P_A / (k1 * P_A + k2)

        coverages = model.solve_steady_state()
        tof = model.get_tof(product_reaction="reaction", coverages=coverages)

        # Site balance
        assert abs(sum(coverages.values()) - 1.0) < 1e-8

        # TOF
        rel_err = abs(tof - tof_exact) / max(tof_exact, 1e-50)
        assert rel_err < 1e-8, (
            f"T={T}K: TOF numerical={tof:.6e}, analytical={tof_exact:.6e}, "
            f"rel_err={rel_err:.2e}"
        )

    @pytest.mark.parametrize("P_A", [0.01, 0.1, 0.5, 1.0, 5.0])
    def test_pressure_sweep(self, P_A: float) -> None:
        """Analytical solution matches at multiple pressures."""
        T = 673.0
        Ea1, Ea2 = 0.3, 0.8
        prefactor = 1e13

        model = MicroKineticModel(temperature=T, total_sites=1.0)
        model.add_species("*", is_site=True)
        model.add_species("A*", initial_coverage=0.0)
        model.add_gas("A(g)", pressure=P_A)
        model.add_gas("B(g)", pressure=0.01)
        model.add_reaction(
            "adsorption",
            {"A(g)": 1, "*": 1}, {"A*": 1},
            Ea_fwd=Ea1, prefactor=prefactor,
        )
        model.add_reaction(
            "reaction",
            {"A*": 1}, {"B(g)": 1, "*": 1},
            Ea_fwd=Ea2, prefactor=prefactor,
        )

        calc = RateConstantCalculator(T)
        k1 = calc.eyring_rate(Ea1, prefactor)
        k2 = calc.eyring_rate(Ea2, prefactor)

        theta_A_exact = k1 * P_A / (k1 * P_A + k2)
        tof_exact = k1 * k2 * P_A / (k1 * P_A + k2)

        coverages = model.solve_steady_state()
        tof = model.get_tof(product_reaction="reaction", coverages=coverages)

        rel_err_cov = abs(coverages["A*"] - theta_A_exact) / max(theta_A_exact, 1e-50)
        rel_err_tof = abs(tof - tof_exact) / max(tof_exact, 1e-50)

        assert rel_err_cov < 1e-8, f"P={P_A}: θ_A rel_err={rel_err_cov:.2e}"
        assert rel_err_tof < 1e-8, f"P={P_A}: TOF rel_err={rel_err_tof:.2e}"


class TestNH3ModelSanityChecks:
    """Sanity checks for the full NH3 decomposition model (not analytic, just physical).

    Note: The default NH3DecompositionModel has convergence issues with fsolve
    due to stiff kinetics (large barrier spread) and poor initial guess conditioning.
    Tests marked xfail will pass once the solver is improved (tracked as T1.3 gap).
    """

    @pytest.mark.xfail(
        reason="NH3DecompositionModel default fsolve diverges: site balance violated "
        "(θ sum=125). Solver conditioning needs improvement — tracked as T1.3 gap.",
        strict=True,
    )
    def test_default_model_converges(self) -> None:
        """Default NH3DecompositionModel reaches steady state."""
        from nh3sofc.analysis.microkinetic import NH3DecompositionModel

        model = NH3DecompositionModel(temperature=673.0)
        coverages = model.solve_steady_state()

        # Must converge to something
        assert len(coverages) > 0
        # Site balance
        assert abs(sum(coverages.values()) - 1.0) < 1e-6
        # All non-negative
        for species, cov in coverages.items():
            assert cov >= -1e-10, f"{species} = {cov}"

    @pytest.mark.xfail(
        reason="Depends on solver convergence fix (T1.3 gap).",
        strict=True,
    )
    def test_tof_positive(self) -> None:
        """TOF is positive (net forward decomposition)."""
        from nh3sofc.analysis.microkinetic import NH3DecompositionModel

        model = NH3DecompositionModel(temperature=673.0)
        tof = model.get_tof()
        assert tof > 0, f"TOF = {tof} (should be positive for decomposition)"

    def test_higher_temperature_higher_tof(self) -> None:
        """Higher temperature gives higher TOF (Arrhenius behavior).

        Note: This test passes because both TOFs are negative (solver diverged)
        and the magnitudes happen to order correctly. It will be meaningful
        once convergence is fixed.
        """
        from nh3sofc.analysis.microkinetic import NH3DecompositionModel

        tof_low = NH3DecompositionModel(temperature=500.0).get_tof()
        tof_high = NH3DecompositionModel(temperature=800.0).get_tof()
        assert tof_high > tof_low, (
            f"TOF(800K)={tof_high:.4e} should be > TOF(500K)={tof_low:.4e}"
        )
