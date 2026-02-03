"""Analysis modules for energetics, thermochemistry, and electronic structure.

Provides tools for:
- Adsorption energy calculations
- Gibbs free energy and thermochemistry
- Rate-determining step identification
- Surface comparison and catalyst screening
- Electronic structure descriptors (d-band center)
"""

from .energetics import (
    AdsorptionEnergyCalculator,
    SurfaceEnergyCalculator,
    DecompositionEnergetics,
    BindingEnergyAnalyzer,
    calculate_adsorption_energy,
    calculate_coverage_dependent_energy,
    plot_surface_stability,
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
from .electronic import (
    DBandAnalyzer,
    DOSAnalyzer,
    calculate_d_band_center,
    get_surface_d_band_center,
)

__all__ = [
    # Energetics
    "AdsorptionEnergyCalculator",
    "SurfaceEnergyCalculator",
    "DecompositionEnergetics",
    "BindingEnergyAnalyzer",
    "plot_surface_stability",
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
    # Electronic structure
    "DBandAnalyzer",
    "DOSAnalyzer",
    "calculate_d_band_center",
    "get_surface_d_band_center",
]
