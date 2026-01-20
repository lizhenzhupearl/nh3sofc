"""Structure building module."""

from .bulk import BulkStructure
from .surface import SurfaceBuilder, SlabStructure, create_slab_from_cif
from .adsorbates import AdsorbatePlacer, filter_unique_configs, save_configs
from .defects import DefectBuilder, OxynitrideStructure, generate_vacancy_series
from .decomposition import (
    DecompositionBuilder,
    generate_decomposition_pathway,
    save_decomposition_configs,
)

__all__ = [
    # Bulk
    "BulkStructure",
    # Surface
    "SurfaceBuilder",
    "SlabStructure",
    "create_slab_from_cif",
    # Adsorbates
    "AdsorbatePlacer",
    "filter_unique_configs",
    "save_configs",
    # Defects
    "DefectBuilder",
    "OxynitrideStructure",
    "generate_vacancy_series",
    # Decomposition
    "DecompositionBuilder",
    "generate_decomposition_pathway",
    "save_decomposition_configs",
]
