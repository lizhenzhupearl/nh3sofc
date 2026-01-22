"""Shared test fixtures for nh3sofc tests."""

import pytest
import numpy as np
from ase import Atoms
from ase.build import bulk, fcc111
from ase.spacegroup import crystal


@pytest.fixture
def perovskite_bulk():
    """Create a simple LaVO3-like perovskite bulk structure."""
    a = 3.9
    atoms = crystal(
        ['La', 'V', 'O'],
        basis=[(0, 0, 0), (0.5, 0.5, 0.5), (0.5, 0.5, 0)],
        spacegroup=221,
        cellpar=[a, a, a, 90, 90, 90]
    )
    return atoms


@pytest.fixture
def rocksalt_bulk():
    """Create NiO rocksalt bulk structure."""
    return bulk('NiO', 'rocksalt', a=4.17)


@pytest.fixture
def fluorite_bulk():
    """Create CeO2 fluorite bulk structure."""
    return bulk('CeO2', 'fluorite', a=5.41)


@pytest.fixture
def simple_slab():
    """Create a simple Pt(111) slab for adsorbate testing."""
    slab = fcc111('Pt', size=(2, 2, 4), vacuum=10.0)
    return slab


@pytest.fixture
def mixed_element_slab():
    """Create a slab with mixed elements (La, V, O, N) for POSCAR testing."""
    # Create a simple structure with multiple element types
    atoms = Atoms(
        symbols=['La', 'V', 'O', 'O', 'N', 'La', 'V', 'O', 'N', 'O'],
        positions=[
            [0, 0, 0], [2, 0, 0], [1, 1, 0], [3, 1, 0], [1, 3, 0],
            [0, 0, 2], [2, 0, 2], [1, 1, 2], [3, 1, 2], [1, 3, 2],
        ],
        cell=[4, 4, 15],
        pbc=True
    )
    return atoms
