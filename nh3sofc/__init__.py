"""
nh3sofc - NH3 Solid Oxide Fuel Cell Simulation Package

A Python package for automating DFT and ML force field calculations
in NH3 solid oxide fuel cell research.
"""

__version__ = "0.1.0"

from .structure.bulk import BulkStructure
from .structure.surface import SurfaceBuilder
from .structure.adsorbates import AdsorbatePlacer
from .core.io import (
    write_poscar,
    write_cif,
    save_structure,
    save_configurations,
    sort_atoms_by_element,
)

__all__ = [
    "BulkStructure",
    "SurfaceBuilder",
    "AdsorbatePlacer",
    "write_poscar",
    "write_cif",
    "save_structure",
    "save_configurations",
    "sort_atoms_by_element",
]
