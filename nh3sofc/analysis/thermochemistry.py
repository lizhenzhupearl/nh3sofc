"""Thermochemistry analysis for Gibbs free energy calculations.

Provides tools for calculating thermodynamic properties including
ZPE, entropy, enthalpy, and Gibbs free energy at various temperatures.
"""

from pathlib import Path
from typing import Optional, List, Union, Dict, Any, Tuple
import numpy as np

from ..core.constants import KB_EV, H_EV_S, C_CMS, EV_TO_J, AVOGADRO


class HarmonicThermo:
    """
    Harmonic thermodynamics for adsorbates.

    Uses harmonic approximation for vibrational contributions.
    Appropriate for strongly bound adsorbates on surfaces.

    G = E_elec + ZPE + H_vib(T) - T*S_vib(T)
    """

    def __init__(
        self,
        frequencies: List[float],
        electronic_energy: float = 0.0,
    ):
        """
        Initialize HarmonicThermo.

        Parameters
        ----------
        frequencies : list
            Vibrational frequencies in cm^-1
        electronic_energy : float
            Electronic energy in eV
        """
        self.frequencies = np.array(frequencies)
        self.electronic_energy = electronic_energy

        # Filter real, positive frequencies
        self.real_frequencies = self.frequencies[
            (self.frequencies > 0) & np.isfinite(self.frequencies)
        ]

    def get_zpe(self) -> float:
        """
        Calculate zero-point energy.

        ZPE = (1/2) * Σ hν

        Returns
        -------
        float
            Zero-point energy in eV
        """
        # Convert cm^-1 to eV: E = hν = h * c * ν_cm
        # h = 4.136e-15 eV·s, c = 2.998e10 cm/s
        zpe = 0.5 * np.sum(self.real_frequencies) * H_EV_S * C_CMS
        return zpe

    def get_vibrational_energy(self, T: float) -> float:
        """
        Calculate vibrational internal energy U_vib(T).

        U_vib = Σ hν * [1/2 + 1/(exp(hν/kT) - 1)]

        Parameters
        ----------
        T : float
            Temperature in K

        Returns
        -------
        float
            Vibrational energy in eV
        """
        if T <= 0:
            return self.get_zpe()

        u_vib = 0.0
        for freq in self.real_frequencies:
            hv = freq * H_EV_S * C_CMS  # Convert to eV
            x = hv / (KB_EV * T)

            if x > 100:  # Avoid overflow
                u_vib += hv / 2  # Just ZPE contribution
            else:
                u_vib += hv * (0.5 + 1.0 / (np.exp(x) - 1))

        return u_vib

    def get_vibrational_entropy(self, T: float) -> float:
        """
        Calculate vibrational entropy S_vib(T).

        S_vib = k * Σ [x/(exp(x)-1) - ln(1-exp(-x))]
        where x = hν/kT

        Parameters
        ----------
        T : float
            Temperature in K

        Returns
        -------
        float
            Vibrational entropy in eV/K
        """
        if T <= 0:
            return 0.0

        s_vib = 0.0
        for freq in self.real_frequencies:
            hv = freq * H_EV_S * C_CMS
            x = hv / (KB_EV * T)

            if x > 100:
                continue  # Negligible contribution

            s_vib += KB_EV * (x / (np.exp(x) - 1) - np.log(1 - np.exp(-x)))

        return s_vib

    def get_vibrational_free_energy(self, T: float) -> float:
        """
        Calculate vibrational Helmholtz free energy F_vib(T).

        F_vib = U_vib - T*S_vib

        Parameters
        ----------
        T : float
            Temperature in K

        Returns
        -------
        float
            Vibrational free energy in eV
        """
        u_vib = self.get_vibrational_energy(T)
        s_vib = self.get_vibrational_entropy(T)
        return u_vib - T * s_vib

    def get_gibbs_energy(self, T: float) -> float:
        """
        Calculate Gibbs free energy for adsorbate.

        G = E_elec + F_vib(T)

        Parameters
        ----------
        T : float
            Temperature in K

        Returns
        -------
        float
            Gibbs free energy in eV
        """
        return self.electronic_energy + self.get_vibrational_free_energy(T)

    def get_heat_capacity(self, T: float) -> float:
        """
        Calculate heat capacity Cv at constant volume.

        Cv = dU/dT = k * Σ x^2 * exp(x) / (exp(x)-1)^2

        Parameters
        ----------
        T : float
            Temperature in K

        Returns
        -------
        float
            Heat capacity in eV/K
        """
        if T <= 0:
            return 0.0

        cv = 0.0
        for freq in self.real_frequencies:
            hv = freq * H_EV_S * C_CMS
            x = hv / (KB_EV * T)

            if x > 100:
                continue

            exp_x = np.exp(x)
            cv += KB_EV * x**2 * exp_x / (exp_x - 1)**2

        return cv


class IdealGasThermo:
    """
    Ideal gas thermodynamics for gas phase molecules.

    Includes translational, rotational, and vibrational contributions.

    G = E_elec + ZPE + H(T) - T*S(T,p)
    """

    def __init__(
        self,
        frequencies: List[float],
        electronic_energy: float,
        mass: float,
        geometry: str = "nonlinear",
        symmetry_number: int = 1,
        spin: int = 0,
    ):
        """
        Initialize IdealGasThermo.

        Parameters
        ----------
        frequencies : list
            Vibrational frequencies in cm^-1
        electronic_energy : float
            Electronic energy in eV
        mass : float
            Molecular mass in amu
        geometry : str
            "linear", "nonlinear", or "monatomic"
        symmetry_number : int
            Rotational symmetry number
        spin : int
            Spin multiplicity minus 1 (0 for singlet)
        """
        self.harmonic = HarmonicThermo(frequencies, electronic_energy)
        self.mass = mass  # amu
        self.geometry = geometry
        self.symmetry = symmetry_number
        self.spin = spin

    def get_translational_entropy(self, T: float, p: float = 1.0) -> float:
        """
        Calculate translational entropy (Sackur-Tetrode equation).

        Parameters
        ----------
        T : float
            Temperature in K
        p : float
            Pressure in bar

        Returns
        -------
        float
            Translational entropy in eV/K
        """
        from scipy.constants import k, h, N_A, pi

        # Mass in kg
        m = self.mass / 1000 / N_A

        # Thermal wavelength
        lambda_th = h / np.sqrt(2 * pi * m * k * T)

        # Volume per molecule at pressure p
        V = k * T / (p * 1e5)  # 1 bar = 1e5 Pa

        # Sackur-Tetrode
        q_trans = V / lambda_th**3
        S_trans = k * (np.log(q_trans) + 5/2)

        return S_trans / EV_TO_J

    def get_rotational_entropy(self, T: float) -> float:
        """
        Calculate rotational entropy.

        Parameters
        ----------
        T : float
            Temperature in K

        Returns
        -------
        float
            Rotational entropy in eV/K
        """
        from scipy.constants import k

        if self.geometry == "monatomic":
            return 0.0

        # Simplified: assumes high-T limit
        if self.geometry == "linear":
            # S_rot = k * (ln(T/σθ_rot) + 1)
            S_rot = k * (np.log(T / self.symmetry) + 1)
        else:
            # S_rot = k * (3/2 * ln(T) - ln(σ) + 3/2 + const)
            S_rot = k * (1.5 * np.log(T / self.symmetry) + 1.5)

        return S_rot / EV_TO_J

    def get_electronic_entropy(self) -> float:
        """
        Calculate electronic entropy.

        Returns
        -------
        float
            Electronic entropy in eV/K
        """
        # S_elec = k * ln(2S+1)
        multiplicity = 2 * self.spin + 1
        return KB_EV * np.log(multiplicity)

    def get_total_entropy(self, T: float, p: float = 1.0) -> float:
        """
        Calculate total entropy.

        Parameters
        ----------
        T : float
            Temperature in K
        p : float
            Pressure in bar

        Returns
        -------
        float
            Total entropy in eV/K
        """
        S_trans = self.get_translational_entropy(T, p)
        S_rot = self.get_rotational_entropy(T)
        S_vib = self.harmonic.get_vibrational_entropy(T)
        S_elec = self.get_electronic_entropy()

        return S_trans + S_rot + S_vib + S_elec

    def get_enthalpy_correction(self, T: float) -> float:
        """
        Calculate enthalpy correction H(T) - H(0).

        Parameters
        ----------
        T : float
            Temperature in K

        Returns
        -------
        float
            Enthalpy correction in eV
        """
        # Translational: 3/2 kT for 3D
        H_trans = 1.5 * KB_EV * T

        # Rotational
        if self.geometry == "monatomic":
            H_rot = 0.0
        elif self.geometry == "linear":
            H_rot = KB_EV * T
        else:
            H_rot = 1.5 * KB_EV * T

        # Vibrational: U_vib - ZPE
        H_vib = self.harmonic.get_vibrational_energy(T) - self.harmonic.get_zpe()

        # pV = kT (ideal gas)
        pV = KB_EV * T

        return H_trans + H_rot + H_vib + pV

    def get_gibbs_energy(self, T: float, p: float = 1.0) -> float:
        """
        Calculate Gibbs free energy.

        G = E_elec + ZPE + H(T) - T*S(T,p)

        Parameters
        ----------
        T : float
            Temperature in K
        p : float
            Pressure in bar

        Returns
        -------
        float
            Gibbs free energy in eV
        """
        E_elec = self.harmonic.electronic_energy
        ZPE = self.harmonic.get_zpe()
        H_corr = self.get_enthalpy_correction(T)
        S = self.get_total_entropy(T, p)

        return E_elec + ZPE + H_corr - T * S


class ReactionThermodynamics:
    """
    Calculate thermodynamics for chemical reactions.

    Handles reactions with gas phase and adsorbed species.
    """

    def __init__(self, temperature: float = 298.15, pressure: float = 1.0):
        """
        Initialize ReactionThermodynamics.

        Parameters
        ----------
        temperature : float
            Temperature in K
        pressure : float
            Pressure in bar
        """
        self.T = temperature
        self.p = pressure

        self.species = {}  # name -> thermo object or dict

    def add_adsorbate(
        self,
        name: str,
        electronic_energy: float,
        frequencies: Optional[List[float]] = None,
    ) -> None:
        """
        Add an adsorbed species.

        Parameters
        ----------
        name : str
            Species name (e.g., "NH3*")
        electronic_energy : float
            Electronic energy in eV
        frequencies : list, optional
            Vibrational frequencies in cm^-1
        """
        if frequencies:
            self.species[name] = HarmonicThermo(frequencies, electronic_energy)
        else:
            self.species[name] = {"E": electronic_energy, "type": "static"}

    def add_gas(
        self,
        name: str,
        electronic_energy: float,
        frequencies: List[float],
        mass: float,
        geometry: str = "nonlinear",
        symmetry: int = 1,
    ) -> None:
        """
        Add a gas phase species.

        Parameters
        ----------
        name : str
            Species name (e.g., "NH3(g)")
        electronic_energy : float
            Electronic energy in eV
        frequencies : list
            Vibrational frequencies in cm^-1
        mass : float
            Molecular mass in amu
        geometry : str
            Molecular geometry
        symmetry : int
            Symmetry number
        """
        self.species[name] = IdealGasThermo(
            frequencies, electronic_energy, mass, geometry, symmetry
        )

    def add_surface(self, name: str, electronic_energy: float) -> None:
        """
        Add a surface (no vibrational correction).

        Parameters
        ----------
        name : str
            Surface name
        electronic_energy : float
            Electronic energy in eV
        """
        self.species[name] = {"E": electronic_energy, "type": "surface"}

    def get_gibbs_energy(
        self,
        name: str,
        T: Optional[float] = None,
        p: Optional[float] = None,
    ) -> float:
        """
        Get Gibbs free energy for a species.

        Parameters
        ----------
        name : str
            Species name
        T : float, optional
            Temperature (default: self.T)
        p : float, optional
            Pressure (default: self.p)

        Returns
        -------
        float
            Gibbs free energy in eV
        """
        T = T if T is not None else self.T
        p = p if p is not None else self.p

        species = self.species[name]

        if isinstance(species, HarmonicThermo):
            return species.get_gibbs_energy(T)
        elif isinstance(species, IdealGasThermo):
            return species.get_gibbs_energy(T, p)
        elif isinstance(species, dict):
            return species["E"]
        else:
            raise ValueError(f"Unknown species type: {name}")

    def calculate_reaction_energy(
        self,
        products: Dict[str, int],
        reactants: Dict[str, int],
        T: Optional[float] = None,
        p: Optional[float] = None,
    ) -> float:
        """
        Calculate reaction Gibbs energy.

        ΔG = Σ(n_i * G_i)_products - Σ(n_i * G_i)_reactants

        Parameters
        ----------
        products : dict
            Product species and stoichiometry
        reactants : dict
            Reactant species and stoichiometry
        T : float, optional
            Temperature
        p : float, optional
            Pressure

        Returns
        -------
        float
            Reaction Gibbs energy in eV
        """
        G_prod = sum(n * self.get_gibbs_energy(s, T, p) for s, n in products.items())
        G_react = sum(n * self.get_gibbs_energy(s, T, p) for s, n in reactants.items())

        return G_prod - G_react

    def get_equilibrium_constant(
        self,
        products: Dict[str, int],
        reactants: Dict[str, int],
        T: Optional[float] = None,
    ) -> float:
        """
        Calculate equilibrium constant K.

        K = exp(-ΔG / kT)

        Parameters
        ----------
        products : dict
            Products and stoichiometry
        reactants : dict
            Reactants and stoichiometry
        T : float, optional
            Temperature

        Returns
        -------
        float
            Equilibrium constant (dimensionless)
        """
        T = T if T is not None else self.T
        dG = self.calculate_reaction_energy(products, reactants, T)

        return np.exp(-dG / (KB_EV * T))


def calculate_zpe(frequencies: List[float]) -> float:
    """
    Calculate zero-point energy from frequencies.

    Parameters
    ----------
    frequencies : list
        Frequencies in cm^-1

    Returns
    -------
    float
        ZPE in eV
    """
    thermo = HarmonicThermo(frequencies)
    return thermo.get_zpe()


def calculate_gibbs_correction(
    frequencies: List[float],
    temperature: float,
    is_gas: bool = False,
    mass: Optional[float] = None,
    pressure: float = 1.0,
) -> float:
    """
    Calculate Gibbs free energy correction (G - E_elec).

    Parameters
    ----------
    frequencies : list
        Frequencies in cm^-1
    temperature : float
        Temperature in K
    is_gas : bool
        Whether species is gas phase
    mass : float, optional
        Molecular mass in amu (required for gas)
    pressure : float
        Pressure in bar

    Returns
    -------
    float
        Gibbs correction in eV
    """
    if is_gas:
        if mass is None:
            raise ValueError("Mass required for gas phase species")
        thermo = IdealGasThermo(frequencies, 0.0, mass)
        return thermo.get_gibbs_energy(temperature, pressure)
    else:
        thermo = HarmonicThermo(frequencies, 0.0)
        return thermo.get_gibbs_energy(temperature)


def get_thermal_properties(
    frequencies: List[float],
    temperatures: List[float],
) -> Dict[str, np.ndarray]:
    """
    Calculate thermal properties at multiple temperatures.

    Parameters
    ----------
    frequencies : list
        Frequencies in cm^-1
    temperatures : list
        Temperatures in K

    Returns
    -------
    dict
        Arrays of ZPE, U, S, F, Cv vs T
    """
    thermo = HarmonicThermo(frequencies)
    temps = np.array(temperatures)

    results = {
        "T": temps,
        "ZPE": np.full_like(temps, thermo.get_zpe()),
        "U": np.array([thermo.get_vibrational_energy(T) for T in temps]),
        "S": np.array([thermo.get_vibrational_entropy(T) for T in temps]),
        "F": np.array([thermo.get_vibrational_free_energy(T) for T in temps]),
        "Cv": np.array([thermo.get_heat_capacity(T) for T in temps]),
    }

    return results
