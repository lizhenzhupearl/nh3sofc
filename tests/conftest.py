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


@pytest.fixture
def ceo2_slab():
    """Create a CeO2 slab-like structure for dopant testing.

    Creates a 2x2x2 fluorite supercell with vacuum to mimic a slab.
    Contains 8 Ce and 16 O atoms.
    """
    # Start with CeO2 bulk
    ceo2 = bulk('CeO2', 'fluorite', a=5.41)
    # Create 2x2x2 supercell
    supercell = ceo2 * (2, 2, 2)
    # Add vacuum in z-direction to make it slab-like
    cell = supercell.cell.copy()
    cell[2, 2] += 15.0  # Add 15 A vacuum
    supercell.set_cell(cell)
    # Shift atoms to center in z
    positions = supercell.get_positions()
    positions[:, 2] += 7.5
    supercell.set_positions(positions)
    return supercell


@pytest.fixture
def zro2_slab():
    """Create a ZrO2 slab-like structure for dopant testing (YSZ, ScSZ).

    Creates a 2x2x2 fluorite supercell with vacuum to mimic a slab.
    Contains 8 Zr and 16 O atoms.
    """
    # Start with ZrO2 bulk (cubic fluorite phase)
    zro2 = bulk('ZrO2', 'fluorite', a=5.07)
    # Create 2x2x2 supercell
    supercell = zro2 * (2, 2, 2)
    # Add vacuum in z-direction to make it slab-like
    cell = supercell.cell.copy()
    cell[2, 2] += 15.0  # Add 15 A vacuum
    supercell.set_cell(cell)
    # Shift atoms to center in z
    positions = supercell.get_positions()
    positions[:, 2] += 7.5
    supercell.set_positions(positions)
    return supercell
