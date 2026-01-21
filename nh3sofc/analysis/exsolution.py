"""Exsolution energetics analysis.

Provides tools for analyzing exsolution processes:
- Vacancy formation energies
- Segregation energies
- Exsolution driving force
- Nanoparticle binding energies
- Comparison with clean surface catalysis
"""

from pathlib import Path
from typing import Optional, List, Union, Dict, Any, Tuple
import numpy as np

from ..core.constants import (
    KB_EV,
    GAS_PHASE_ENERGIES,
    EXSOLUTION_METALS,
    OPERATING_CONDITIONS,
)


class ExsolutionEnergetics:
    """
    Analyze energetics of exsolution processes.

    Calculates:
    - Vacancy formation energies
    - Segregation energies
    - Exsolution driving force
    - Nanoparticle formation/binding energies
    - Metal-support interaction strength

    Examples
    --------
    >>> energetics = ExsolutionEnergetics()
    >>> energetics.set_reference_energies(
    ...     e_pristine=-200.0,
    ...     e_bulk_metal_per_atom=-5.51,
    ... )
    >>>
    >>> # Add calculation results
    >>> energetics.add_stage("pristine", -200.0)
    >>> energetics.add_stage("defective", -195.0)
    >>> energetics.add_stage("exsolved", -210.0, n_particle_atoms=13)
    >>>
    >>> # Calculate energetics
    >>> e_exs = energetics.calculate_exsolution_energy()
    >>> print(f"Exsolution driving force: {e_exs:.2f} eV")
    """

    def __init__(
        self,
        metal: str = "Ni",
        temperature: float = 873.0,
        p_o2: float = 1e-20,
    ):
        """
        Initialize ExsolutionEnergetics.

        Parameters
        ----------
        metal : str
            Exsolution metal (Ni, Co, Fe)
        temperature : float
            Temperature in K (default: 873 K = 600 C)
        p_o2 : float
            Oxygen partial pressure in atm (default: 1e-20 atm, reducing)
        """
        self.metal = metal
        self.temperature = temperature
        self.p_o2 = p_o2

        # Reference energies
        self.e_pristine = None
        self.e_bulk_metal_per_atom = EXSOLUTION_METALS.get(metal, {}).get(
            "bulk_energy", -5.5
        )
        self.e_o2 = GAS_PHASE_ENERGIES.get("O2", -9.86)

        # Stage energies
        self.stages = {}

    def set_reference_energies(
        self,
        e_pristine: Optional[float] = None,
        e_bulk_metal_per_atom: Optional[float] = None,
        e_o2: Optional[float] = None,
    ) -> None:
        """
        Set reference energies for calculations.

        Parameters
        ----------
        e_pristine : float, optional
            Energy of pristine perovskite surface
        e_bulk_metal_per_atom : float, optional
            Bulk metal energy per atom
        e_o2 : float, optional
            O2 molecule energy
        """
        if e_pristine is not None:
            self.e_pristine = e_pristine
        if e_bulk_metal_per_atom is not None:
            self.e_bulk_metal_per_atom = e_bulk_metal_per_atom
        if e_o2 is not None:
            self.e_o2 = e_o2

    def add_stage(
        self,
        stage: str,
        energy: float,
        n_particle_atoms: int = 0,
        n_o_vacancies: int = 0,
        structure: Optional[Any] = None,
    ) -> None:
        """
        Add energy for a calculation stage.

        Parameters
        ----------
        stage : str
            Stage name (pristine, defective, segregated, exsolved)
        energy : float
            Total energy in eV
        n_particle_atoms : int
            Number of atoms in nanoparticle (for exsolved stage)
        n_o_vacancies : int
            Number of oxygen vacancies
        structure : Atoms, optional
            Structure for this stage
        """
        self.stages[stage] = {
            "energy": energy,
            "n_particle_atoms": n_particle_atoms,
            "n_o_vacancies": n_o_vacancies,
            "structure": structure,
        }

        # Auto-set pristine reference
        if stage == "pristine" and self.e_pristine is None:
            self.e_pristine = energy

    def get_oxygen_chemical_potential(
        self,
        temperature: Optional[float] = None,
        p_o2: Optional[float] = None,
    ) -> float:
        """
        Calculate oxygen chemical potential at given conditions.

        μ_O(T, p) = 0.5 * [E(O₂) + μ°(T) + kT*ln(p/p°)]

        Parameters
        ----------
        temperature : float, optional
            Temperature in K (default: self.temperature)
        p_o2 : float, optional
            O2 partial pressure in atm (default: self.p_o2)

        Returns
        -------
        float
            Oxygen chemical potential in eV
        """
        T = temperature if temperature is not None else self.temperature
        p = p_o2 if p_o2 is not None else self.p_o2

        # Standard chemical potential correction (approximate)
        # This should be calculated more rigorously from thermochemical tables
        # Here using simplified form: μ°(T) ≈ -T*S°/2 where S° ≈ 0.002 eV/K for O2
        mu_std_correction = -T * 0.001  # Simplified entropy term

        # Pressure correction
        p_ref = 1.0  # Reference pressure (1 atm)
        mu_pressure = KB_EV * T * np.log(max(p / p_ref, 1e-30))

        mu_O = 0.5 * (self.e_o2 + mu_std_correction + mu_pressure)

        return mu_O

    def calculate_vacancy_formation_energy(
        self,
        e_defective: Optional[float] = None,
        e_pristine: Optional[float] = None,
        n_vacancies: int = 1,
        include_chemical_potential: bool = True,
    ) -> float:
        """
        Calculate oxygen vacancy formation energy.

        E_vac = [E(defective) - E(pristine) + n * μ_O] / n

        Parameters
        ----------
        e_defective : float, optional
            Energy of defective structure
        e_pristine : float, optional
            Energy of pristine structure
        n_vacancies : int
            Number of oxygen vacancies
        include_chemical_potential : bool
            If True, include O chemical potential correction

        Returns
        -------
        float
            Vacancy formation energy per vacancy in eV
        """
        e_def = e_defective or self.stages.get("defective", {}).get("energy")
        e_pris = e_pristine or self.e_pristine

        if e_def is None or e_pris is None:
            raise ValueError("Missing energies for vacancy calculation")

        e_vac = e_def - e_pris

        if include_chemical_potential:
            mu_O = self.get_oxygen_chemical_potential()
            e_vac += n_vacancies * mu_O

        return e_vac / n_vacancies

    def calculate_segregation_energy(
        self,
        e_segregated: Optional[float] = None,
        e_bulk_distributed: Optional[float] = None,
    ) -> float:
        """
        Calculate surface segregation energy.

        E_seg = E(surface_segregated) - E(bulk_distributed)

        Parameters
        ----------
        e_segregated : float, optional
            Energy with metal segregated to surface
        e_bulk_distributed : float, optional
            Energy with metal distributed in bulk

        Returns
        -------
        float
            Segregation energy in eV (negative = favorable)
        """
        e_seg = e_segregated or self.stages.get("segregated", {}).get("energy")
        e_bulk = e_bulk_distributed or self.stages.get("defective", {}).get("energy")

        if e_seg is None or e_bulk is None:
            raise ValueError("Missing energies for segregation calculation")

        return e_seg - e_bulk

    def calculate_exsolution_energy(
        self,
        e_exsolved: Optional[float] = None,
        e_substrate: Optional[float] = None,
        n_atoms: Optional[int] = None,
        e_bulk_metal: Optional[float] = None,
        include_vacancy_correction: bool = True,
    ) -> float:
        """
        Calculate exsolution driving force.

        E_exs = E(with_particle) + E(defective_substrate) - E(pristine) - n*E_bulk_metal

        Parameters
        ----------
        e_exsolved : float, optional
            Energy of structure with exsolved particle
        e_substrate : float, optional
            Energy of defective substrate (vacancies left after exsolution)
        n_atoms : int, optional
            Number of atoms in nanoparticle
        e_bulk_metal : float, optional
            Bulk metal energy per atom
        include_vacancy_correction : bool
            Include oxygen chemical potential correction

        Returns
        -------
        float
            Exsolution energy in eV (negative = favorable)
        """
        e_exs = e_exsolved or self.stages.get("exsolved", {}).get("energy")
        e_pris = self.e_pristine

        if e_exs is None or e_pris is None:
            raise ValueError("Missing energies for exsolution calculation")

        # Get particle size
        if n_atoms is None:
            n_atoms = self.stages.get("exsolved", {}).get("n_particle_atoms", 0)

        # Get bulk metal energy
        if e_bulk_metal is None:
            e_bulk_metal = self.e_bulk_metal_per_atom

        # Simple exsolution energy (particle formation from pristine)
        e_exsolution = e_exs - e_pris - n_atoms * e_bulk_metal

        # Vacancy correction
        if include_vacancy_correction:
            n_vac = self.stages.get("exsolved", {}).get("n_o_vacancies", 0)
            if n_vac > 0:
                mu_O = self.get_oxygen_chemical_potential()
                e_exsolution += n_vac * mu_O

        return e_exsolution

    def calculate_particle_binding_energy(
        self,
        e_system: float,
        e_surface: float,
        e_isolated_particle: float,
    ) -> float:
        """
        Calculate nanoparticle-surface binding energy.

        E_bind = E(particle/surface) - E(surface) - E(isolated_particle)

        Parameters
        ----------
        e_system : float
            Energy of combined particle-surface system
        e_surface : float
            Energy of clean/defective surface
        e_isolated_particle : float
            Energy of isolated metal cluster

        Returns
        -------
        float
            Binding energy in eV (negative = stable)
        """
        return e_system - e_surface - e_isolated_particle

    def get_exsolution_driving_force(
        self,
        temperature: Optional[float] = None,
        p_o2: Optional[float] = None,
    ) -> Dict[str, float]:
        """
        Calculate temperature and pressure dependent exsolution driving force.

        Parameters
        ----------
        temperature : float, optional
            Temperature in K
        p_o2 : float, optional
            Oxygen partial pressure in atm

        Returns
        -------
        dict
            Dictionary with Gibbs free energy and components
        """
        T = temperature if temperature is not None else self.temperature
        p = p_o2 if p_o2 is not None else self.p_o2

        # Get base exsolution energy (DFT)
        try:
            e_exs = self.calculate_exsolution_energy(include_vacancy_correction=False)
        except ValueError:
            return {"error": "Missing stage energies"}

        # Oxygen chemical potential correction
        n_vac = self.stages.get("exsolved", {}).get("n_o_vacancies", 0)
        mu_O = self.get_oxygen_chemical_potential(T, p)
        vacancy_term = n_vac * mu_O

        # Configurational entropy (simplified)
        # ΔS_config ≈ k_B * ln(Ω) where Ω is number of configurations
        # This is a rough approximation
        n_atoms = self.stages.get("exsolved", {}).get("n_particle_atoms", 1)
        s_config = KB_EV * np.log(max(n_atoms, 1))  # Very simplified
        entropy_term = -T * s_config

        # Total Gibbs free energy
        g_exsolution = e_exs + vacancy_term + entropy_term

        return {
            "delta_G": g_exsolution,
            "delta_E": e_exs,
            "vacancy_term": vacancy_term,
            "entropy_term": entropy_term,
            "mu_O": mu_O,
            "temperature": T,
            "p_o2": p,
            "favorable": g_exsolution < 0,
        }

    def compare_with_clean_surface(
        self,
        exsolved_energies: Dict[str, float],
        clean_surface_energies: Dict[str, float],
    ) -> Dict[str, Any]:
        """
        Compare catalytic activity of exsolved particle vs clean surface.

        Parameters
        ----------
        exsolved_energies : dict
            Adsorption/reaction energies on exsolved particle system
            Keys: reaction intermediates (e.g., "NH3", "NH2", "NH", "N")
        clean_surface_energies : dict
            Same energies on clean perovskite surface

        Returns
        -------
        dict
            Comparison results with activity differences
        """
        comparison = {
            "exsolved": exsolved_energies,
            "clean_surface": clean_surface_energies,
            "differences": {},
            "more_favorable_on_exsolved": [],
            "more_favorable_on_clean": [],
        }

        for key in set(exsolved_energies.keys()) & set(clean_surface_energies.keys()):
            e_exs = exsolved_energies[key]
            e_clean = clean_surface_energies[key]
            diff = e_exs - e_clean

            comparison["differences"][key] = diff

            if diff < -0.1:  # Threshold for meaningful difference
                comparison["more_favorable_on_exsolved"].append(key)
            elif diff > 0.1:
                comparison["more_favorable_on_clean"].append(key)

        return comparison

    def print_summary(self) -> None:
        """Print summary of exsolution energetics."""
        print(f"\n{'='*50}")
        print(f"Exsolution Energetics Summary ({self.metal})")
        print(f"{'='*50}")
        print(f"Temperature: {self.temperature} K")
        print(f"p(O2): {self.p_o2:.2e} atm")
        print(f"μ_O: {self.get_oxygen_chemical_potential():.3f} eV")
        print()

        print("Stage Energies:")
        for stage, data in self.stages.items():
            print(f"  {stage}: {data['energy']:.4f} eV")
        print()

        try:
            e_vac = self.calculate_vacancy_formation_energy()
            print(f"Vacancy formation energy: {e_vac:.3f} eV/vacancy")
        except ValueError:
            pass

        try:
            e_seg = self.calculate_segregation_energy()
            print(f"Segregation energy: {e_seg:.3f} eV")
        except ValueError:
            pass

        try:
            result = self.get_exsolution_driving_force()
            if "error" not in result:
                print(f"Exsolution ΔG: {result['delta_G']:.3f} eV")
                print(f"  Favorable: {'Yes' if result['favorable'] else 'No'}")
        except ValueError:
            pass

        print(f"{'='*50}\n")


def calculate_exsolution_driving_force(
    e_pristine: float,
    e_exsolved: float,
    n_particle_atoms: int,
    metal: str = "Ni",
    n_o_vacancies: int = 0,
    temperature: float = 873.0,
    p_o2: float = 1e-20,
) -> Dict[str, float]:
    """
    Convenience function to calculate exsolution driving force.

    Parameters
    ----------
    e_pristine : float
        Energy of pristine perovskite surface (eV)
    e_exsolved : float
        Energy of structure with exsolved particle (eV)
    n_particle_atoms : int
        Number of atoms in nanoparticle
    metal : str
        Exsolution metal (Ni, Co, Fe)
    n_o_vacancies : int
        Number of oxygen vacancies created
    temperature : float
        Temperature in K
    p_o2 : float
        Oxygen partial pressure in atm

    Returns
    -------
    dict
        Exsolution driving force and components
    """
    energetics = ExsolutionEnergetics(
        metal=metal,
        temperature=temperature,
        p_o2=p_o2,
    )

    energetics.set_reference_energies(e_pristine=e_pristine)
    energetics.add_stage(
        "exsolved",
        e_exsolved,
        n_particle_atoms=n_particle_atoms,
        n_o_vacancies=n_o_vacancies,
    )

    return energetics.get_exsolution_driving_force()
