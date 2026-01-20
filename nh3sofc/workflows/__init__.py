"""Workflow modules for automated calculations.

Provides high-level workflows for common DFT/ML calculations:
- Geometry optimization (relaxation)
- NH3 decomposition pathway
- NEB transition states
- Frequency/thermochemistry
"""

from .relaxation import (
    RelaxationWorkflow,
    BatchRelaxation,
    relax_structure,
)
from .decomposition import (
    DecompositionWorkflow,
    H2FormationWorkflow,
    run_decomposition_study,
)
from .neb import (
    NEBWorkflow,
    AutoNEB,
    setup_neb_calculation,
    create_neb_from_pathway,
)
from .frequency import (
    FrequencyWorkflow,
    GasPhaseThermo,
    ThermochemistryWorkflow,
    run_frequency_calculation,
)
from .screening import (
    ScreeningWorkflow,
    CompositionScreening,
    AdsorbateScreening,
    run_screening,
)

__all__ = [
    # Relaxation
    "RelaxationWorkflow",
    "BatchRelaxation",
    "relax_structure",
    # Decomposition
    "DecompositionWorkflow",
    "H2FormationWorkflow",
    "run_decomposition_study",
    # NEB
    "NEBWorkflow",
    "AutoNEB",
    "setup_neb_calculation",
    "create_neb_from_pathway",
    # Frequency
    "FrequencyWorkflow",
    "GasPhaseThermo",
    "ThermochemistryWorkflow",
    "run_frequency_calculation",
    # Screening
    "ScreeningWorkflow",
    "CompositionScreening",
    "AdsorbateScreening",
    "run_screening",
]
