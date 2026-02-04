"""Structure building module."""

from .bulk import BulkStructure
from .surface import (
    SurfaceBuilder,
    SlabStructure,
    create_slab_from_cif,
)
from .adsorbates import AdsorbatePlacer, filter_unique_configs, save_configs
from .defects import (
    DefectBuilder,
    OxynitrideStructure,
    generate_vacancy_series,
    analyze_defect_distribution,
    analyze_oxynitride_pool,
    print_defect_analysis,
    plot_defect_distribution,
)
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
from .dopants import (
    DopantBuilder,
    DopedCeriaStructure,
    analyze_dopant_distribution,
    generate_dopant_series,
    print_dopant_analysis,
)

__all__ = [
    # Bulk
    "BulkStructure",
    # Surface - General
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
    "analyze_defect_distribution",
    "analyze_oxynitride_pool",
    "print_defect_analysis",
    "plot_defect_distribution",
    # Decomposition
    "DecompositionBuilder",
    "generate_decomposition_pathway",
    "save_decomposition_configs",
    # Exsolution
    "ExsolutionBuilder",
    "ExsolutionStructure",
    "generate_exsolution_series",
    "create_metallic_cluster",
    # Dopants
    "DopantBuilder",
    "DopedCeriaStructure",
    "analyze_dopant_distribution",
    "generate_dopant_series",
    "print_dopant_analysis",
]
