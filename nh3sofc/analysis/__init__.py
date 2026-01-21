"""Analysis modules for energetics and thermochemistry.

Provides tools for:
- Adsorption energy calculations
- Gibbs free energy and thermochemistry
- Rate-determining step identification
- Surface comparison and catalyst screening
"""

from .energetics import (
    AdsorptionEnergyCalculator,
    SurfaceEnergyCalculator,
    DecompositionEnergetics,
    BindingEnergyAnalyzer,
    calculate_adsorption_energy,
    calculate_coverage_dependent_energy,
)
from .thermochemistry import (
    HarmonicThermo,
    IdealGasThermo,
    ReactionThermodynamics,
    calculate_zpe,
    calculate_gibbs_correction,
    get_thermal_properties,
)
from .rds import (
    EnergySpanModel,
    BEPRelation,
    RDSAnalyzer,
    find_thermodynamic_rds,
    estimate_barrier_bep,
    calculate_rate_constant,
)
from .surface_comparison import (
    SurfaceComparator,
    ActivityDescriptor,
    compare_surfaces,
    get_best_catalyst,
)
from .microkinetic import (
    RateConstantCalculator,
    MicroKineticModel,
    NH3DecompositionModel,
    calculate_tof,
    arrhenius_analysis,
)
from .exsolution import (
    ExsolutionEnergetics,
    calculate_exsolution_driving_force,
)

__all__ = [
    # Energetics
    "AdsorptionEnergyCalculator",
    "SurfaceEnergyCalculator",
    "DecompositionEnergetics",
    "BindingEnergyAnalyzer",
    "calculate_adsorption_energy",
    "calculate_coverage_dependent_energy",
    # Thermochemistry
    "HarmonicThermo",
    "IdealGasThermo",
    "ReactionThermodynamics",
    "calculate_zpe",
    "calculate_gibbs_correction",
    "get_thermal_properties",
    # RDS
    "EnergySpanModel",
    "BEPRelation",
    "RDSAnalyzer",
    "find_thermodynamic_rds",
    "estimate_barrier_bep",
    "calculate_rate_constant",
    # Surface comparison
    "SurfaceComparator",
    "ActivityDescriptor",
    "compare_surfaces",
    "get_best_catalyst",
    # Microkinetic
    "RateConstantCalculator",
    "MicroKineticModel",
    "NH3DecompositionModel",
    "calculate_tof",
    "arrhenius_analysis",
    # Exsolution
    "ExsolutionEnergetics",
    "calculate_exsolution_driving_force",
]
