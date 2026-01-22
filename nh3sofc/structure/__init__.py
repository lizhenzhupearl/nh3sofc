"""Structure building module."""

from .bulk import BulkStructure
from .surface import (
    SurfaceBuilder,
    SlabStructure,
    create_slab_from_cif,
    PerovskiteSurfaceBuilder,
    RocksaltSurfaceBuilder,
    FluoriteSurfaceBuilder,
)
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
    # Surface - General
    "SurfaceBuilder",
    "SlabStructure",
    "create_slab_from_cif",
    # Surface - Specialized builders
    "PerovskiteSurfaceBuilder",
    "RocksaltSurfaceBuilder",
    "FluoriteSurfaceBuilder",
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
