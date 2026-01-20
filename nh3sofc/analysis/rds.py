"""Rate-determining step (RDS) identification.

Provides tools for identifying the rate-determining step in catalytic
reactions using various approaches:
- Thermodynamic approximation (highest ΔE step)
- Energy span model (Kozuch-Shaik)
- Brønsted-Evans-Polanyi (BEP) relations
"""

from typing import Optional, List, Dict, Any, Tuple
import numpy as np

from ..core.constants import KB_EV


class EnergySpanModel:
    """
    Energy span model for RDS identification.

    Based on Kozuch and Shaik's energy span model:
    δE = T_TDTS - I_TDI  (if TDTS appears after TDI)
    δE = T_TDTS - I_TDI + ΔGr  (if TDTS appears before TDI)

    where TDTS = TOF-determining transition state
          TDI = TOF-determining intermediate
    """

    def __init__(
        self,
        intermediates: Dict[str, float],
        transition_states: Optional[Dict[str, float]] = None,
        reaction_energy: float = 0.0,
    ):
        """
        Initialize EnergySpanModel.

        Parameters
        ----------
        intermediates : dict
            Intermediate energies {name: energy}
        transition_states : dict, optional
            Transition state energies {name: energy}
        reaction_energy : float
            Overall reaction energy (for cyclic catalysis)
        """
        self.intermediates = intermediates
        self.transition_states = transition_states or {}
        self.dG_rxn = reaction_energy

        # Store ordered species
        self._order_species()

    def _order_species(self) -> None:
        """Order species by assumed reaction coordinate."""
        # Assume intermediates are named in order or use provided order
        self.int_names = list(self.intermediates.keys())
        self.ts_names = list(self.transition_states.keys())

    def get_energy_span(self) -> Tuple[float, str, str]:
        """
        Calculate energy span.

        Returns
        -------
        tuple
            (energy_span, TDTS_name, TDI_name)
        """
        # Find TDI (lowest intermediate)
        tdi_name = min(self.intermediates, key=self.intermediates.get)
        tdi_energy = self.intermediates[tdi_name]

        # Find TDTS (highest TS or intermediate if no TS provided)
        if self.transition_states:
            tdts_name = max(self.transition_states, key=self.transition_states.get)
            tdts_energy = self.transition_states[tdts_name]
        else:
            # Use highest intermediate as approximate TDTS
            tdts_name = max(self.intermediates, key=self.intermediates.get)
            tdts_energy = self.intermediates[tdts_name]

        # Calculate energy span
        # For cyclic catalysis, need to consider position in cycle
        tdi_idx = self.int_names.index(tdi_name) if tdi_name in self.int_names else 0

        if self.ts_names and tdts_name in self.ts_names:
            # Estimate TS position (between intermediates)
            ts_idx = self.ts_names.index(tdts_name)
            tdts_after_tdi = ts_idx >= tdi_idx
        else:
            tdts_idx = self.int_names.index(tdts_name) if tdts_name in self.int_names else 0
            tdts_after_tdi = tdts_idx >= tdi_idx

        if tdts_after_tdi:
            delta_e = tdts_energy - tdi_energy
        else:
            delta_e = tdts_energy - tdi_energy + self.dG_rxn

        return delta_e, tdts_name, tdi_name

    def get_apparent_activation_energy(self) -> float:
        """
        Get apparent activation energy (energy span).

        Returns
        -------
        float
            Apparent activation energy in eV
        """
        delta_e, _, _ = self.get_energy_span()
        return delta_e

    def get_degree_of_rate_control(self, step: str) -> float:
        """
        Calculate degree of rate control for a step.

        X_RC = (∂ln(r)/∂(-G/RT)) at constant G of other states

        Parameters
        ----------
        step : str
            Step name (intermediate or TS)

        Returns
        -------
        float
            Degree of rate control (0 to 1)
        """
        delta_e, tdts, tdi = self.get_energy_span()

        # Full DRC requires numerical derivatives
        # Simplified: TDTS and TDI have highest DRC
        if step == tdts:
            return 1.0
        elif step == tdi:
            return 1.0
        else:
            return 0.0

    def estimate_tof(self, temperature: float, prefactor: float = 1e13) -> float:
        """
        Estimate turnover frequency.

        TOF ≈ A * exp(-δE / kT)

        Parameters
        ----------
        temperature : float
            Temperature in K
        prefactor : float
            Pre-exponential factor (s^-1)

        Returns
        -------
        float
            Estimated TOF in s^-1
        """
        delta_e = self.get_apparent_activation_energy()
        return prefactor * np.exp(-delta_e / (KB_EV * temperature))


class BEPRelation:
    """
    Brønsted-Evans-Polanyi (BEP) relations.

    E_a = α * ΔE + β

    where E_a = activation energy
          ΔE = reaction energy
          α, β = BEP parameters
    """

    # Common BEP parameters from literature
    BEP_PARAMETERS = {
        # (alpha, beta) for different reaction types
        "N-H_dissociation": (0.87, 0.95),  # NH3 → NH2 + H
        "C-H_dissociation": (0.75, 0.90),
        "O-H_dissociation": (0.65, 0.85),
        "N-N_formation": (0.50, 1.20),
        "H-H_formation": (0.60, 0.80),
        "default": (0.80, 1.00),
    }

    def __init__(
        self,
        alpha: Optional[float] = None,
        beta: Optional[float] = None,
        reaction_type: str = "default",
    ):
        """
        Initialize BEPRelation.

        Parameters
        ----------
        alpha : float, optional
            BEP slope parameter
        beta : float, optional
            BEP intercept parameter
        reaction_type : str
            Reaction type for default parameters
        """
        if alpha is not None and beta is not None:
            self.alpha = alpha
            self.beta = beta
        elif reaction_type in self.BEP_PARAMETERS:
            self.alpha, self.beta = self.BEP_PARAMETERS[reaction_type]
        else:
            self.alpha, self.beta = self.BEP_PARAMETERS["default"]

    def estimate_barrier(self, reaction_energy: float) -> float:
        """
        Estimate activation barrier from reaction energy.

        E_a = α * ΔE + β

        Parameters
        ----------
        reaction_energy : float
            Reaction energy in eV

        Returns
        -------
        float
            Estimated activation energy in eV
        """
        e_a = self.alpha * reaction_energy + self.beta

        # Barrier cannot be negative
        return max(0.0, e_a)

    def estimate_reverse_barrier(self, reaction_energy: float) -> float:
        """
        Estimate reverse activation barrier.

        E_a_rev = E_a - ΔE

        Parameters
        ----------
        reaction_energy : float
            Reaction energy in eV

        Returns
        -------
        float
            Reverse activation energy in eV
        """
        e_a_fwd = self.estimate_barrier(reaction_energy)
        e_a_rev = e_a_fwd - reaction_energy

        return max(0.0, e_a_rev)

    @classmethod
    def fit_from_data(
        cls,
        reaction_energies: List[float],
        activation_energies: List[float],
    ) -> "BEPRelation":
        """
        Fit BEP parameters from data.

        Parameters
        ----------
        reaction_energies : list
            Reaction energies
        activation_energies : list
            Activation energies

        Returns
        -------
        BEPRelation
            Fitted BEP relation
        """
        dE = np.array(reaction_energies)
        Ea = np.array(activation_energies)

        # Linear fit
        coeffs = np.polyfit(dE, Ea, 1)
        alpha, beta = coeffs[0], coeffs[1]

        return cls(alpha=alpha, beta=beta)


class RDSAnalyzer:
    """
    Comprehensive RDS analysis for catalytic reactions.

    Combines multiple approaches to identify the rate-determining step.
    """

    def __init__(self):
        """Initialize RDSAnalyzer."""
        self.pathway = {}  # step_name -> energy
        self.barriers = {}  # step_name -> barrier
        self.bep = BEPRelation()

    def set_pathway(self, pathway: Dict[str, float]) -> None:
        """
        Set reaction pathway energies.

        Parameters
        ----------
        pathway : dict
            Step energies relative to reference
        """
        self.pathway = pathway.copy()

    def set_barriers(self, barriers: Dict[str, float]) -> None:
        """
        Set known activation barriers.

        Parameters
        ----------
        barriers : dict
            Activation barriers for each step
        """
        self.barriers = barriers.copy()

    def find_rds_thermodynamic(self) -> Tuple[str, float]:
        """
        Find RDS using thermodynamic approximation.

        Returns step with highest reaction energy.

        Returns
        -------
        tuple
            (step_name, reaction_energy)
        """
        steps = list(self.pathway.keys())
        reaction_energies = {}

        for i in range(len(steps) - 1):
            step_name = f"{steps[i]}→{steps[i+1]}"
            dE = self.pathway[steps[i+1]] - self.pathway[steps[i]]
            reaction_energies[step_name] = dE

        if not reaction_energies:
            return ("", 0.0)

        rds = max(reaction_energies.items(), key=lambda x: x[1])
        return rds

    def find_rds_bep(self) -> Tuple[str, float]:
        """
        Find RDS using BEP-estimated barriers.

        Returns step with highest estimated barrier.

        Returns
        -------
        tuple
            (step_name, estimated_barrier)
        """
        steps = list(self.pathway.keys())
        estimated_barriers = {}

        for i in range(len(steps) - 1):
            step_name = f"{steps[i]}→{steps[i+1]}"
            dE = self.pathway[steps[i+1]] - self.pathway[steps[i]]

            # Use known barrier if available, otherwise estimate
            if step_name in self.barriers:
                barrier = self.barriers[step_name]
            else:
                barrier = self.bep.estimate_barrier(dE)

            estimated_barriers[step_name] = barrier

        if not estimated_barriers:
            return ("", 0.0)

        rds = max(estimated_barriers.items(), key=lambda x: x[1])
        return rds

    def find_rds_energy_span(self) -> Tuple[str, float, str, str]:
        """
        Find RDS using energy span model.

        Returns
        -------
        tuple
            (energy_span, TDTS, TDI)
        """
        # Estimate TS energies using BEP if not provided
        ts_energies = {}
        steps = list(self.pathway.keys())

        for i in range(len(steps) - 1):
            ts_name = f"TS_{steps[i]}_{steps[i+1]}"
            initial_e = self.pathway[steps[i]]
            dE = self.pathway[steps[i+1]] - initial_e

            if f"{steps[i]}→{steps[i+1]}" in self.barriers:
                barrier = self.barriers[f"{steps[i]}→{steps[i+1]}"]
            else:
                barrier = self.bep.estimate_barrier(dE)

            ts_energies[ts_name] = initial_e + barrier

        model = EnergySpanModel(
            self.pathway,
            ts_energies,
            reaction_energy=self.pathway.get(steps[-1], 0) - self.pathway.get(steps[0], 0),
        )

        span, tdts, tdi = model.get_energy_span()
        return span, tdts, tdi

    def get_all_barriers(self) -> Dict[str, Dict[str, float]]:
        """
        Get all barriers (known and estimated).

        Returns
        -------
        dict
            Barriers for each step with source
        """
        steps = list(self.pathway.keys())
        all_barriers = {}

        for i in range(len(steps) - 1):
            step_name = f"{steps[i]}→{steps[i+1]}"
            dE = self.pathway[steps[i+1]] - self.pathway[steps[i]]

            if step_name in self.barriers:
                all_barriers[step_name] = {
                    "barrier": self.barriers[step_name],
                    "source": "calculated",
                    "reaction_energy": dE,
                }
            else:
                all_barriers[step_name] = {
                    "barrier": self.bep.estimate_barrier(dE),
                    "source": "BEP_estimated",
                    "reaction_energy": dE,
                }

        return all_barriers

    def print_analysis(self) -> None:
        """Print comprehensive RDS analysis."""
        print("\n" + "=" * 60)
        print("Rate-Determining Step Analysis")
        print("=" * 60)

        # Thermodynamic RDS
        thermo_rds, thermo_dE = self.find_rds_thermodynamic()
        print(f"\n1. Thermodynamic RDS (highest ΔE):")
        print(f"   Step: {thermo_rds}")
        print(f"   ΔE = {thermo_dE:.3f} eV")

        # BEP RDS
        bep_rds, bep_barrier = self.find_rds_bep()
        print(f"\n2. BEP-estimated RDS (highest E_a):")
        print(f"   Step: {bep_rds}")
        print(f"   E_a ≈ {bep_barrier:.3f} eV")

        # Energy span
        span, tdts, tdi = self.find_rds_energy_span()
        print(f"\n3. Energy Span Model:")
        print(f"   Energy span: {span:.3f} eV")
        print(f"   TDTS: {tdts}")
        print(f"   TDI: {tdi}")

        # All barriers
        print("\n4. Step-by-step barriers:")
        barriers = self.get_all_barriers()
        for step, info in barriers.items():
            source = "calc" if info["source"] == "calculated" else "BEP"
            print(f"   {step}: E_a = {info['barrier']:.3f} eV ({source}), "
                  f"ΔE = {info['reaction_energy']:.3f} eV")

        print("=" * 60)


def find_thermodynamic_rds(
    pathway: Dict[str, float],
) -> Tuple[str, float]:
    """
    Find rate-determining step using thermodynamic approximation.

    Parameters
    ----------
    pathway : dict
        Step energies {step_name: energy}

    Returns
    -------
    tuple
        (rds_step, reaction_energy)
    """
    analyzer = RDSAnalyzer()
    analyzer.set_pathway(pathway)
    return analyzer.find_rds_thermodynamic()


def estimate_barrier_bep(
    reaction_energy: float,
    reaction_type: str = "default",
) -> float:
    """
    Estimate activation barrier using BEP relation.

    Parameters
    ----------
    reaction_energy : float
        Reaction energy in eV
    reaction_type : str
        Type of reaction

    Returns
    -------
    float
        Estimated barrier in eV
    """
    bep = BEPRelation(reaction_type=reaction_type)
    return bep.estimate_barrier(reaction_energy)


def calculate_rate_constant(
    barrier: float,
    temperature: float,
    prefactor: float = 1e13,
) -> float:
    """
    Calculate rate constant using Eyring equation.

    k = A * exp(-E_a / kT)

    Parameters
    ----------
    barrier : float
        Activation energy in eV
    temperature : float
        Temperature in K
    prefactor : float
        Pre-exponential factor (s^-1)

    Returns
    -------
    float
        Rate constant in s^-1
    """
    return prefactor * np.exp(-barrier / (KB_EV * temperature))
