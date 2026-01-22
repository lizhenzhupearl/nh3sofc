"""Regression tests for bugs that have been fixed.

Each test documents a bug that was found and fixed.
If these tests fail, we've reintroduced the bug.
"""

import pytest
import tempfile
import numpy as np


class TestPOSCARAtomOrdering:
    """
    Bug: POSCAR had scattered element ordering like 'La V N O La V O'
    which VESTA couldn't read.

    Fix: VASPInputGenerator.generate_all() now sorts atoms by element
    before writing POSCAR.

    Date fixed: 2026-01-22
    """

    def test_poscar_elements_are_grouped(self, mixed_element_slab):
        """POSCAR should have grouped elements like 'La N O V' not scattered."""
        from nh3sofc.calculators.vasp import VASPInputGenerator

        with tempfile.TemporaryDirectory() as tmpdir:
            vasp = VASPInputGenerator(
                mixed_element_slab,
                calc_type='relax',
                work_dir=tmpdir
            )
            files = vasp.generate_all(encut=400)

            with open(files['POSCAR']) as f:
                lines = f.readlines()

            # Line 6 (index 5) has element symbols
            elements_line = lines[5].strip()
            elements = elements_line.split()

            # Elements should appear only once each (grouped)
            assert len(elements) == len(set(elements)), \
                f"Elements are scattered: {elements_line}. Expected grouped elements."

    def test_poscar_counts_match_atoms(self, mixed_element_slab):
        """POSCAR atom counts should match actual structure."""
        from nh3sofc.calculators.vasp import VASPInputGenerator

        with tempfile.TemporaryDirectory() as tmpdir:
            vasp = VASPInputGenerator(
                mixed_element_slab,
                calc_type='relax',
                work_dir=tmpdir
            )
            files = vasp.generate_all(encut=400)

            with open(files['POSCAR']) as f:
                lines = f.readlines()

            # Line 7 (index 6) has counts
            counts_line = lines[6].strip()
            counts = [int(x) for x in counts_line.split()]

            # Total should match number of atoms
            assert sum(counts) == len(mixed_element_slab), \
                f"POSCAR counts {sum(counts)} != atoms {len(mixed_element_slab)}"


class TestAdsorbatePlacerAddsAtoms:
    """
    Bug: AdsorbatePlacer appeared to not add molecules to surface.

    Root cause: POSCAR ordering issue made VESTA unable to read the file.
    The molecules were actually being added correctly.

    This test verifies molecules are actually added.

    Date fixed: 2026-01-22
    """

    def test_nh3_actually_added_to_surface(self, simple_slab):
        """NH3 should actually be present in the resulting structure."""
        from nh3sofc.structure import AdsorbatePlacer

        original_count = len(simple_slab)
        original_symbols = set(simple_slab.get_chemical_symbols())

        placer = AdsorbatePlacer(simple_slab)
        result = placer.add_simple('NH3', position=(0.5, 0.5), height=2.0)

        # Check atom count increased by 4 (NH3 = N + 3H)
        assert len(result) == original_count + 4, \
            f"Expected {original_count + 4} atoms, got {len(result)}"

        # Check N and H are now present
        new_symbols = set(result.get_chemical_symbols())
        assert 'N' in new_symbols, "N should be present after adding NH3"
        assert 'H' in new_symbols, "H should be present after adding NH3"

    def test_adsorbate_preserved_in_vasp_output(self, simple_slab):
        """Adsorbate should still be present after VASP input generation."""
        from nh3sofc.structure import AdsorbatePlacer
        from nh3sofc.calculators.vasp import VASPInputGenerator

        placer = AdsorbatePlacer(simple_slab)
        with_nh3 = placer.add_simple('NH3', position=(0.5, 0.5), height=2.0)

        with tempfile.TemporaryDirectory() as tmpdir:
            vasp = VASPInputGenerator(with_nh3, calc_type='relax', work_dir=tmpdir)
            files = vasp.generate_all(encut=400)

            with open(files['POSCAR']) as f:
                lines = f.readlines()

            elements_line = lines[5].strip()
            counts_line = lines[6].strip()

            # H and N should be in POSCAR
            assert 'H' in elements_line, "H should be in POSCAR elements"
            assert 'N' in elements_line, "N should be in POSCAR elements"

            # Total count should match
            counts = [int(x) for x in counts_line.split()]
            assert sum(counts) == len(with_nh3)


class TestSupercellPreservesConstraints:
    """
    Bug: repeat_xy() didn't preserve fixed atom constraints.

    Fix: repeat_xy() now maps old fixed indices to new structure.

    Date fixed: 2026-01-22
    """

    def test_repeat_xy_preserves_fixed_atoms(self, rocksalt_bulk):
        """Fixed atoms should be preserved when creating supercell."""
        from nh3sofc.structure import SurfaceBuilder

        builder = SurfaceBuilder(rocksalt_bulk)
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=4,
            fix_bottom=2
        )

        original_fixed = len(slab.fixed_indices)
        assert original_fixed > 0, "Should have fixed atoms before repeat"

        # Repeat 2x2
        repeated = slab.repeat_xy(2, 2)

        # Should have 4x the fixed atoms
        expected_fixed = original_fixed * 4
        assert len(repeated.fixed_indices) == expected_fixed, \
            f"Expected {expected_fixed} fixed atoms, got {len(repeated.fixed_indices)}"


class TestSpecializedBuildersAutoDetect:
    """
    Feature: Specialized builders should auto-detect cation/anion sites.

    This ensures the auto-detection works correctly.
    """

    def test_perovskite_auto_detects_sites(self, perovskite_bulk):
        """PerovskiteSurfaceBuilder should auto-detect A and B sites."""
        from nh3sofc.structure import PerovskiteSurfaceBuilder

        builder = PerovskiteSurfaceBuilder(perovskite_bulk)

        assert builder.A_site == 'La', f"Expected La, got {builder.A_site}"
        assert builder.B_site == 'V', f"Expected V, got {builder.B_site}"

    def test_rocksalt_auto_detects_cation(self, rocksalt_bulk):
        """RocksaltSurfaceBuilder should auto-detect cation."""
        from nh3sofc.structure import RocksaltSurfaceBuilder

        builder = RocksaltSurfaceBuilder(rocksalt_bulk)

        assert builder.cation == 'Ni', f"Expected Ni, got {builder.cation}"
        assert builder.anion == 'O', f"Expected O, got {builder.anion}"

    def test_fluorite_auto_detects_cation(self, fluorite_bulk):
        """FluoriteSurfaceBuilder should auto-detect cation."""
        from nh3sofc.structure import FluoriteSurfaceBuilder

        builder = FluoriteSurfaceBuilder(fluorite_bulk)

        assert builder.cation == 'Ce', f"Expected Ce, got {builder.cation}"
        assert builder.anion == 'O', f"Expected O, got {builder.anion}"
