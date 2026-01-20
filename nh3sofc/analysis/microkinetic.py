"""Microkinetic modeling for catalytic reactions.

Provides tools for steady-state analysis and turnover frequency
calculations based on DFT-derived energetics.
"""

from typing import Optional, List, Dict, Any, Tuple, Callable
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve

from ..core.constants import KB_EV, H_EV_S


class RateConstantCalculator:
    """
    Calculate rate constants from activation energies.

    Uses transition state theory (Eyring equation):
    k = (kT/h) * exp(-Ea/kT)
    """

    def __init__(self, temperature: float = 673.0):
        """
        Initialize RateConstantCalculator.

        Parameters
        ----------
        temperature : float
            Temperature in K
        """
        self.T = temperature

    def eyring_rate(
        self,
        barrier: float,
        prefactor: Optional[float] = None,
    ) -> float:
        """
        Calculate rate constant using Eyring equation.

        k = A * exp(-Ea/kT)

        Parameters
        ----------
        barrier : float
            Activation energy in eV
        prefactor : float, optional
            Pre-exponential factor (s^-1). Default: kT/h

        Returns
        -------
        float
            Rate constant in s^-1
        """
        if prefactor is None:
            # kT/h prefactor
            prefactor = KB_EV * self.T / H_EV_S

        if barrier < 0:
            barrier = 0.0  # Barrierless

        return prefactor * np.exp(-barrier / (KB_EV * self.T))

    def equilibrium_constant(self, dG: float) -> float:
        """
        Calculate equilibrium constant.

        K = exp(-ΔG/kT)

        Parameters
        ----------
        dG : float
            Gibbs free energy change in eV

        Returns
        -------
        float
            Equilibrium constant (dimensionless)
        """
        return np.exp(-dG / (KB_EV * self.T))

    def reverse_rate(
        self,
        k_forward: float,
        dG: float,
    ) -> float:
        """
        Calculate reverse rate constant from forward rate and ΔG.

        k_rev = k_fwd / K_eq

        Parameters
        ----------
        k_forward : float
            Forward rate constant
        dG : float
            Reaction Gibbs energy

        Returns
        -------
        float
            Reverse rate constant
        """
        K_eq = self.equilibrium_constant(dG)
        return k_forward / K_eq if K_eq > 0 else 0.0


class MicroKineticModel:
    """
    Microkinetic model for surface reactions.

    Solves steady-state coverage equations and calculates
    turnover frequencies.

    Examples
    --------
    >>> model = MicroKineticModel(temperature=673)
    >>> model.add_species("*", is_site=True)
    >>> model.add_species("NH3*")
    >>> model.add_reaction("NH3_ads", {"NH3(g)": 1, "*": 1}, {"NH3*": 1}, Ea_fwd=0.0)
    >>> coverages = model.solve_steady_state()
    >>> tof = model.get_tof()
    """

    def __init__(
        self,
        temperature: float = 673.0,
        total_sites: float = 1.0,
    ):
        """
        Initialize MicroKineticModel.

        Parameters
        ----------
        temperature : float
            Temperature in K
        total_sites : float
            Total site concentration (normalized to 1)
        """
        self.T = temperature
        self.total_sites = total_sites
        self.rate_calc = RateConstantCalculator(temperature)

        self.species = {}  # name -> {"is_site": bool, "initial": float}
        self.reactions = {}  # name -> reaction info
        self.gas_pressures = {}  # gas species -> pressure (bar)

    def add_species(
        self,
        name: str,
        is_site: bool = False,
        initial_coverage: float = 0.0,
    ) -> None:
        """
        Add a surface species.

        Parameters
        ----------
        name : str
            Species name
        is_site : bool
            Whether this is an empty site
        initial_coverage : float
            Initial coverage
        """
        self.species[name] = {
            "is_site": is_site,
            "initial": initial_coverage,
        }

        if is_site:
            self.species[name]["initial"] = self.total_sites

    def add_gas(self, name: str, pressure: float = 1.0) -> None:
        """
        Add a gas phase species.

        Parameters
        ----------
        name : str
            Gas species name
        pressure : float
            Partial pressure in bar
        """
        self.gas_pressures[name] = pressure

    def add_reaction(
        self,
        name: str,
        reactants: Dict[str, int],
        products: Dict[str, int],
        Ea_fwd: float,
        Ea_rev: Optional[float] = None,
        dG: Optional[float] = None,
        prefactor: float = 1e13,
    ) -> None:
        """
        Add a reaction.

        Parameters
        ----------
        name : str
            Reaction name
        reactants : dict
            Reactant species and stoichiometry
        products : dict
            Product species and stoichiometry
        Ea_fwd : float
            Forward activation energy (eV)
        Ea_rev : float, optional
            Reverse activation energy (eV)
        dG : float, optional
            Reaction Gibbs energy (eV)
        prefactor : float
            Pre-exponential factor (s^-1)
        """
        # Calculate rate constants
        k_fwd = self.rate_calc.eyring_rate(Ea_fwd, prefactor)

        if Ea_rev is not None:
            k_rev = self.rate_calc.eyring_rate(Ea_rev, prefactor)
        elif dG is not None:
            k_rev = self.rate_calc.reverse_rate(k_fwd, dG)
        else:
            # Assume irreversible
            k_rev = 0.0

        self.reactions[name] = {
            "reactants": reactants,
            "products": products,
            "Ea_fwd": Ea_fwd,
            "Ea_rev": Ea_rev,
            "dG": dG,
            "k_fwd": k_fwd,
            "k_rev": k_rev,
            "prefactor": prefactor,
        }

    def get_rate(
        self,
        reaction_name: str,
        coverages: Dict[str, float],
    ) -> Tuple[float, float]:
        """
        Calculate forward and reverse rates for a reaction.

        Parameters
        ----------
        reaction_name : str
            Reaction name
        coverages : dict
            Current coverages

        Returns
        -------
        tuple
            (forward_rate, reverse_rate)
        """
        rxn = self.reactions[reaction_name]

        # Forward rate
        r_fwd = rxn["k_fwd"]
        for species, stoich in rxn["reactants"].items():
            if species in coverages:
                r_fwd *= coverages[species] ** stoich
            elif species in self.gas_pressures:
                r_fwd *= self.gas_pressures[species] ** stoich

        # Reverse rate
        r_rev = rxn["k_rev"]
        for species, stoich in rxn["products"].items():
            if species in coverages:
                r_rev *= coverages[species] ** stoich
            elif species in self.gas_pressures:
                r_rev *= self.gas_pressures[species] ** stoich

        return r_fwd, r_rev

    def get_net_rate(
        self,
        reaction_name: str,
        coverages: Dict[str, float],
    ) -> float:
        """Get net rate (forward - reverse)."""
        r_fwd, r_rev = self.get_rate(reaction_name, coverages)
        return r_fwd - r_rev

    def _coverage_odes(
        self,
        y: np.ndarray,
        t: float,
        species_list: List[str],
    ) -> np.ndarray:
        """ODE system for coverage evolution."""
        coverages = dict(zip(species_list, y))

        dydt = np.zeros(len(y))

        for rxn_name, rxn in self.reactions.items():
            r_net = self.get_net_rate(rxn_name, coverages)

            # Update species based on stoichiometry
            for i, species in enumerate(species_list):
                # Consumption in reactants
                if species in rxn["reactants"]:
                    dydt[i] -= rxn["reactants"][species] * r_net

                # Production in products
                if species in rxn["products"]:
                    dydt[i] += rxn["products"][species] * r_net

        return dydt

    def _steady_state_equations(
        self,
        y: np.ndarray,
        species_list: List[str],
        site_species: str,
    ) -> np.ndarray:
        """Steady-state equations (dydt = 0)."""
        # Get time derivatives
        dydt = self._coverage_odes(y, 0, species_list)

        # Replace site conservation equation
        site_idx = species_list.index(site_species)

        # Site conservation: sum of all coverages = total_sites
        total_coverage = sum(y)
        dydt[site_idx] = total_coverage - self.total_sites

        return dydt

    def solve_steady_state(
        self,
        method: str = "fsolve",
        **kwargs,
    ) -> Dict[str, float]:
        """
        Solve for steady-state coverages.

        Parameters
        ----------
        method : str
            Solution method ("fsolve" or "ode")
        **kwargs : dict
            Additional solver parameters

        Returns
        -------
        dict
            Steady-state coverages
        """
        # Get surface species (not gas)
        species_list = [s for s in self.species if s not in self.gas_pressures]

        # Find site species
        site_species = None
        for s, info in self.species.items():
            if info.get("is_site", False):
                site_species = s
                break

        if site_species is None:
            raise ValueError("No site species defined")

        # Initial guess
        y0 = np.array([
            self.species[s]["initial"] for s in species_list
        ])

        # Ensure initial guess is valid
        y0 = np.maximum(y0, 1e-10)
        y0 = y0 / y0.sum() * self.total_sites

        if method == "fsolve":
            solution = fsolve(
                self._steady_state_equations,
                y0,
                args=(species_list, site_species),
                full_output=True,
                **kwargs,
            )
            y_ss = solution[0]
        else:
            # Integrate to steady state
            t = np.linspace(0, 1e6, 10000)
            solution = odeint(self._coverage_odes, y0, t, args=(species_list,))
            y_ss = solution[-1]

        # Ensure non-negative
        y_ss = np.maximum(y_ss, 0.0)

        return dict(zip(species_list, y_ss))

    def get_tof(
        self,
        product_reaction: Optional[str] = None,
        coverages: Optional[Dict[str, float]] = None,
    ) -> float:
        """
        Calculate turnover frequency.

        Parameters
        ----------
        product_reaction : str, optional
            Reaction that produces the product
        coverages : dict, optional
            Coverages (calculates steady-state if not provided)

        Returns
        -------
        float
            TOF in s^-1
        """
        if coverages is None:
            coverages = self.solve_steady_state()

        if product_reaction is None:
            # Use last reaction
            product_reaction = list(self.reactions.keys())[-1]

        return self.get_net_rate(product_reaction, coverages)

    def sensitivity_analysis(
        self,
        parameter: str,
        perturbation: float = 0.01,
    ) -> Dict[str, float]:
        """
        Perform sensitivity analysis.

        Parameters
        ----------
        parameter : str
            Parameter to perturb (reaction name for barrier)
        perturbation : float
            Relative perturbation

        Returns
        -------
        dict
            Sensitivity coefficients
        """
        # Get baseline TOF
        tof_base = self.get_tof()

        sensitivities = {}

        for rxn_name in self.reactions:
            # Perturb forward barrier
            rxn = self.reactions[rxn_name]
            Ea_orig = rxn["Ea_fwd"]

            # Increase barrier
            rxn["Ea_fwd"] = Ea_orig * (1 + perturbation)
            rxn["k_fwd"] = self.rate_calc.eyring_rate(rxn["Ea_fwd"], rxn["prefactor"])

            tof_plus = self.get_tof()

            # Restore
            rxn["Ea_fwd"] = Ea_orig
            rxn["k_fwd"] = self.rate_calc.eyring_rate(Ea_orig, rxn["prefactor"])

            # Calculate sensitivity
            if tof_base > 0:
                sens = (tof_plus - tof_base) / tof_base / perturbation
            else:
                sens = 0.0

            sensitivities[rxn_name] = sens

        return sensitivities

    def print_summary(
        self,
        coverages: Optional[Dict[str, float]] = None,
    ) -> None:
        """Print model summary."""
        if coverages is None:
            coverages = self.solve_steady_state()

        print("\n" + "=" * 60)
        print("Microkinetic Model Summary")
        print("=" * 60)

        print(f"\nTemperature: {self.T} K")

        print("\nSteady-State Coverages:")
        for species, cov in sorted(coverages.items()):
            print(f"  {species:15s}: {cov:.6f}")

        print("\nReaction Rates:")
        for rxn_name in self.reactions:
            r_net = self.get_net_rate(rxn_name, coverages)
            print(f"  {rxn_name:20s}: {r_net:.4e} s^-1")

        tof = self.get_tof(coverages=coverages)
        print(f"\nTurnover Frequency: {tof:.4e} s^-1")

        print("=" * 60)


class NH3DecompositionModel(MicroKineticModel):
    """
    Microkinetic model for NH3 decomposition.

    NH3* → NH2* + H* → NH* + 2H* → N* + 3H* → N2 + 3H2

    Simplified model with key elementary steps.
    """

    def __init__(
        self,
        temperature: float = 673.0,
        barriers: Optional[Dict[str, float]] = None,
        reaction_energies: Optional[Dict[str, float]] = None,
    ):
        """
        Initialize NH3DecompositionModel.

        Parameters
        ----------
        temperature : float
            Temperature in K
        barriers : dict, optional
            Activation barriers for each step
        reaction_energies : dict, optional
            Reaction energies for each step
        """
        super().__init__(temperature)

        # Default barriers (eV) - should be from DFT
        default_barriers = {
            "NH3_ads": 0.0,
            "NH3_NH2": 1.2,
            "NH2_NH": 1.0,
            "NH_N": 0.8,
            "N_N2": 1.5,
            "H_H2": 0.6,
        }

        # Default reaction energies (eV)
        default_dG = {
            "NH3_ads": -0.5,
            "NH3_NH2": 0.3,
            "NH2_NH": 0.2,
            "NH_N": -0.1,
            "N_N2": -1.0,
            "H_H2": 0.2,
        }

        self.barriers = barriers or default_barriers
        self.reaction_energies = reaction_energies or default_dG

        self._setup_model()

    def _setup_model(self) -> None:
        """Set up species and reactions."""
        # Surface species
        self.add_species("*", is_site=True)
        self.add_species("NH3*", initial_coverage=0.0)
        self.add_species("NH2*", initial_coverage=0.0)
        self.add_species("NH*", initial_coverage=0.0)
        self.add_species("N*", initial_coverage=0.0)
        self.add_species("H*", initial_coverage=0.0)

        # Gas species
        self.add_gas("NH3(g)", pressure=0.1)
        self.add_gas("N2(g)", pressure=0.01)
        self.add_gas("H2(g)", pressure=0.01)

        # Reactions
        # NH3 adsorption
        self.add_reaction(
            "NH3_ads",
            {"NH3(g)": 1, "*": 1},
            {"NH3*": 1},
            Ea_fwd=self.barriers["NH3_ads"],
            dG=self.reaction_energies["NH3_ads"],
        )

        # NH3* → NH2* + H*
        self.add_reaction(
            "NH3_NH2",
            {"NH3*": 1, "*": 1},
            {"NH2*": 1, "H*": 1},
            Ea_fwd=self.barriers["NH3_NH2"],
            dG=self.reaction_energies["NH3_NH2"],
        )

        # NH2* → NH* + H*
        self.add_reaction(
            "NH2_NH",
            {"NH2*": 1, "*": 1},
            {"NH*": 1, "H*": 1},
            Ea_fwd=self.barriers["NH2_NH"],
            dG=self.reaction_energies["NH2_NH"],
        )

        # NH* → N* + H*
        self.add_reaction(
            "NH_N",
            {"NH*": 1, "*": 1},
            {"N*": 1, "H*": 1},
            Ea_fwd=self.barriers["NH_N"],
            dG=self.reaction_energies["NH_N"],
        )

        # 2N* → N2 + 2*
        self.add_reaction(
            "N_N2",
            {"N*": 2},
            {"N2(g)": 1, "*": 2},
            Ea_fwd=self.barriers["N_N2"],
            dG=self.reaction_energies["N_N2"],
        )

        # 2H* → H2 + 2*
        self.add_reaction(
            "H_H2",
            {"H*": 2},
            {"H2(g)": 1, "*": 2},
            Ea_fwd=self.barriers["H_H2"],
            dG=self.reaction_energies["H_H2"],
        )


def calculate_tof(
    barriers: Dict[str, float],
    reaction_energies: Dict[str, float],
    temperature: float = 673.0,
) -> float:
    """
    Calculate TOF for NH3 decomposition.

    Parameters
    ----------
    barriers : dict
        Activation barriers
    reaction_energies : dict
        Reaction energies
    temperature : float
        Temperature in K

    Returns
    -------
    float
        Turnover frequency in s^-1
    """
    model = NH3DecompositionModel(
        temperature=temperature,
        barriers=barriers,
        reaction_energies=reaction_energies,
    )

    return model.get_tof()


def arrhenius_analysis(
    barriers: Dict[str, float],
    reaction_energies: Dict[str, float],
    temperatures: List[float] = None,
) -> Dict[str, Any]:
    """
    Perform Arrhenius analysis over temperature range.

    Parameters
    ----------
    barriers : dict
        Activation barriers
    reaction_energies : dict
        Reaction energies
    temperatures : list, optional
        Temperatures to analyze

    Returns
    -------
    dict
        Arrhenius analysis results
    """
    if temperatures is None:
        temperatures = np.linspace(400, 900, 20)

    tofs = []
    for T in temperatures:
        model = NH3DecompositionModel(
            temperature=T,
            barriers=barriers,
            reaction_energies=reaction_energies,
        )
        tofs.append(model.get_tof())

    tofs = np.array(tofs)

    # Arrhenius fit: ln(TOF) = ln(A) - Ea/kT
    inv_T = 1.0 / np.array(temperatures)
    ln_tof = np.log(np.maximum(tofs, 1e-50))

    # Linear fit
    valid = np.isfinite(ln_tof)
    if np.sum(valid) > 2:
        coeffs = np.polyfit(inv_T[valid], ln_tof[valid], 1)
        Ea_apparent = -coeffs[0] * KB_EV
        prefactor = np.exp(coeffs[1])
    else:
        Ea_apparent = 0.0
        prefactor = 0.0

    return {
        "temperatures": temperatures,
        "tofs": tofs,
        "apparent_Ea": Ea_apparent,
        "prefactor": prefactor,
    }
