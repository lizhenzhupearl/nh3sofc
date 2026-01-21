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
from .exsolution import (
    ExsolutionBuilder,
    ExsolutionStructure,
    generate_exsolution_series,
    create_metallic_cluster,
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
    # Exsolution
    "ExsolutionBuilder",
    "ExsolutionStructure",
    "generate_exsolution_series",
    "create_metallic_cluster",
]
