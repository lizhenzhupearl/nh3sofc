"""Regression tests for bugs that have been fixed.

Each test documents a bug that was found and fixed.
If these tests fail, we've reintroduced the bug.
"""

import pytest
import tempfile
import numpy as np
from pathlib import Path


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


class TestLayerToleranceAutoCalculation:
    """
    Bug: Fixed 0.5 A tolerance in identify_layers() splits V and O atoms
    in LaVO3 VO2 layers into separate layers because they have z-position
    differences of ~0.5-0.6 A.

    Fix: Added estimate_layer_tolerance() method that calculates tolerance
    from covalent radii (0.5 * min_bond_length). Changed default tolerance
    to "auto" which uses this calculation.

    Date fixed: 2026-01-22
    """

    def test_auto_tolerance_is_calculated_from_radii(self, perovskite_bulk):
        """Auto tolerance should be calculated from covalent radii."""
        from nh3sofc.structure import SurfaceBuilder, SlabStructure

        builder = SurfaceBuilder(perovskite_bulk)
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=6,
            vacuum=15.0,
        )

        tol = slab.estimate_layer_tolerance()

        # For LaVO3: V-O covalent radii sum = ~2.19 A
        # Expected tolerance ~1.1 A (0.5 * 2.19)
        assert tol > 0.8, f"Tolerance {tol} too small for perovskite"
        assert tol < 1.5, f"Tolerance {tol} unexpectedly large"

    def test_identify_layers_uses_auto_tolerance_by_default(self, perovskite_bulk):
        """identify_layers() should use auto tolerance by default."""
        from nh3sofc.structure import SurfaceBuilder

        builder = SurfaceBuilder(perovskite_bulk)
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=6,
            vacuum=15.0,
        )

        # With auto tolerance, should get fewer, properly-grouped layers
        layers_auto = slab.identify_layers()  # default is "auto"
        layers_small = slab.identify_layers(tolerance=0.3)  # small tolerance

        # Small tolerance should create more (incorrectly split) layers
        assert len(layers_small) >= len(layers_auto), \
            "Auto tolerance should group atoms better than small tolerance"

    def test_auto_tolerance_groups_perovskite_layers_correctly(self, perovskite_bulk):
        """V and O in VO2 layer should be grouped together with auto tolerance."""
        from nh3sofc.structure import SurfaceBuilder

        builder = SurfaceBuilder(perovskite_bulk)
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=8,
            vacuum=15.0,
        )

        layers = slab.identify_layers()

        # Check that we have layers with both V and O (VO2-like)
        vo2_layers = [
            layer for layer in layers
            if 'V' in layer['composition'] and 'O' in layer['composition']
        ]
        assert len(vo2_layers) > 0, \
            "Should have VO2-like layers with both V and O"

        # Check that we have layers with both La and O (LaO-like)
        lao_layers = [
            layer for layer in layers
            if 'La' in layer['composition'] and 'O' in layer['composition']
            and 'V' not in layer['composition']
        ]
        assert len(lao_layers) > 0, \
            "Should have LaO-like layers with La and O but not V"


class TestSymmetricSlabTrimming:
    """
    Bug: create_symmetric_slab() did not create truly symmetric slabs.
    Top layer had different composition than bottom layer (e.g., VO2 on
    top vs O-only on bottom).

    Fix: Added trim_to_symmetric_termination() method that creates an
    oversized slab and trims to layers matching target termination.
    Updated create_symmetric_slab() to use this approach.

    Date fixed: 2026-01-22
    """

    def test_symmetric_slab_has_matching_terminations(self, perovskite_bulk):
        """Top and bottom layers should have matching element ratios."""
        from nh3sofc.structure import SurfaceBuilder

        builder = SurfaceBuilder(perovskite_bulk)
        slab = builder.create_symmetric_slab(
            miller_index=(0, 0, 1),
            layers=7,
            vacuum=15.0,
            fix_bottom=0,  # Don't fix to simplify testing
        )

        layers = slab.identify_layers()
        top_comp = layers[-1]["composition"]
        bottom_comp = layers[0]["composition"]

        # Get element ratios
        def get_ratios(comp):
            if not comp:
                return {}
            min_val = min(comp.values())
            return {k: v / min_val for k, v in comp.items()}

        top_ratios = get_ratios(top_comp)
        bottom_ratios = get_ratios(bottom_comp)

        # Must have same elements
        assert set(top_ratios.keys()) == set(bottom_ratios.keys()), \
            f"Top {top_comp} and bottom {bottom_comp} have different elements"

        # Ratios should match (within tolerance)
        for elem in top_ratios:
            assert abs(top_ratios[elem] - bottom_ratios[elem]) < 0.1, \
                f"Element {elem} ratio differs: top={top_ratios[elem]}, bottom={bottom_ratios[elem]}"

    def test_trim_to_lao_termination(self, perovskite_bulk):
        """Should be able to create LaO-terminated symmetric slab."""
        from nh3sofc.structure import SurfaceBuilder

        builder = SurfaceBuilder(perovskite_bulk)

        # First create a regular slab to trim
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=10,
            vacuum=15.0,
        )

        # Trim to LaO termination
        symmetric = slab.trim_to_symmetric_termination(
            termination={"La": 1, "O": 1},
            min_layers=5,
        )

        layers = symmetric.identify_layers()
        top_comp = layers[-1]["composition"]
        bottom_comp = layers[0]["composition"]

        # Both should have La and O, with same ratio
        assert "La" in top_comp and "O" in top_comp, \
            f"Top layer {top_comp} should have La and O"
        assert "La" in bottom_comp and "O" in bottom_comp, \
            f"Bottom layer {bottom_comp} should have La and O"

        # V should NOT be in the termination layers
        assert "V" not in top_comp, f"Top layer should not have V: {top_comp}"
        assert "V" not in bottom_comp, f"Bottom layer should not have V: {bottom_comp}"

    def test_trim_to_vo2_termination(self, perovskite_bulk):
        """Should be able to create VO2-terminated symmetric slab."""
        from nh3sofc.structure import SurfaceBuilder

        builder = SurfaceBuilder(perovskite_bulk)

        # First create a regular slab to trim
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=10,
            vacuum=15.0,
        )

        # Trim to VO2 termination
        symmetric = slab.trim_to_symmetric_termination(
            termination={"V": 1, "O": 2},
            min_layers=5,
        )

        layers = symmetric.identify_layers()
        top_comp = layers[-1]["composition"]
        bottom_comp = layers[0]["composition"]

        # Both should have V and O
        assert "V" in top_comp and "O" in top_comp, \
            f"Top layer {top_comp} should have V and O"
        assert "V" in bottom_comp and "O" in bottom_comp, \
            f"Bottom layer {bottom_comp} should have V and O"

        # Check V:O ratio is approximately 1:2
        top_ratio = top_comp.get("O", 0) / top_comp.get("V", 1)
        bottom_ratio = bottom_comp.get("O", 0) / bottom_comp.get("V", 1)
        assert 1.5 < top_ratio < 2.5, f"Top V:O ratio wrong: {top_ratio}"
        assert 1.5 < bottom_ratio < 2.5, f"Bottom V:O ratio wrong: {bottom_ratio}"

    def test_symmetric_slab_auto_finds_best_termination(self, perovskite_bulk):
        """create_symmetric_slab() without termination should auto-find valid one."""
        from nh3sofc.structure import SurfaceBuilder

        builder = SurfaceBuilder(perovskite_bulk)

        # Without specifying termination, should still create symmetric slab
        slab = builder.create_symmetric_slab(
            miller_index=(0, 0, 1),
            layers=7,
            vacuum=15.0,
            fix_bottom=0,
        )

        layers = slab.identify_layers()
        top_comp = layers[-1]["composition"]
        bottom_comp = layers[0]["composition"]

        # Get element ratios
        def get_ratios(comp):
            if not comp:
                return {}
            min_val = min(comp.values())
            return {k: v / min_val for k, v in comp.items()}

        top_ratios = get_ratios(top_comp)
        bottom_ratios = get_ratios(bottom_comp)

        # Must have same elements and ratios
        assert set(top_ratios.keys()) == set(bottom_ratios.keys()), \
            f"Top {top_comp} and bottom {bottom_comp} have different elements"
        for elem in top_ratios:
            assert abs(top_ratios[elem] - bottom_ratios[elem]) < 0.1, \
                f"Element {elem} ratio differs"

    def test_termination_ratio_equivalence(self, perovskite_bulk):
        """{"La": 1, "O": 1} and {"La": 8, "O": 8} should give same results."""
        from nh3sofc.structure import SurfaceBuilder

        builder = SurfaceBuilder(perovskite_bulk)

        # Create slab first
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=10,
            vacuum=15.0,
        )

        # Trim with ratio 1:1
        sym1 = slab.trim_to_symmetric_termination(
            termination={"La": 1, "O": 1},
            min_layers=5,
        )

        # Trim with ratio 8:8 (same as 1:1)
        sym8 = slab.trim_to_symmetric_termination(
            termination={"La": 8, "O": 8},
            min_layers=5,
        )

        # Should have same number of atoms
        assert len(sym1.atoms) == len(sym8.atoms), \
            f"Different atom counts: {len(sym1.atoms)} vs {len(sym8.atoms)}"

        # Should have same layer count
        layers1 = sym1.identify_layers()
        layers8 = sym8.identify_layers()
        assert len(layers1) == len(layers8), \
            f"Different layer counts: {len(layers1)} vs {len(layers8)}"


class TestPOSCARStandardFormat:
    """
    Bug: ASE's default POSCAR writer can produce non-standard format with
    repeated atom types like "La O N La N O O" instead of unique labels.

    Fix: Added write_poscar() function in core/io.py that writes proper
    VASP5 format with unique element labels only.

    Date fixed: 2026-01-22
    """

    def test_write_poscar_has_unique_elements(self, mixed_element_slab):
        """write_poscar should produce only unique element labels."""
        from nh3sofc import write_poscar

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = f"{tmpdir}/test.vasp"
            write_poscar(mixed_element_slab, filepath)

            with open(filepath) as f:
                lines = f.readlines()

            # Line 6 (index 5) has element symbols
            elements_line = lines[5].strip()
            elements = elements_line.split()

            # Elements should appear only once each
            assert len(elements) == len(set(elements)), \
                f"Elements repeated: {elements_line}. Expected unique labels only."

    def test_write_poscar_preserves_atom_count(self, mixed_element_slab):
        """write_poscar should preserve total atom count."""
        from nh3sofc import write_poscar

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = f"{tmpdir}/test.vasp"
            write_poscar(mixed_element_slab, filepath)

            with open(filepath) as f:
                lines = f.readlines()

            # Line 7 (index 6) has counts
            counts_line = lines[6].strip()
            counts = [int(x) for x in counts_line.split()]

            assert sum(counts) == len(mixed_element_slab), \
                f"POSCAR counts {sum(counts)} != atoms {len(mixed_element_slab)}"

    def test_save_configurations_multiple_formats(self, simple_slab):
        """save_configurations should save in multiple formats."""
        from nh3sofc import save_configurations

        configs = [simple_slab, simple_slab.copy()]

        with tempfile.TemporaryDirectory() as tmpdir:
            result = save_configurations(
                configs,
                tmpdir,
                name_prefix="config",
                formats=["poscar", "cif"]
            )

            assert "work_dir" in result
            assert "configs" in result
            assert len(result["configs"]) == 2
            for p in result["configs"]:
                assert "poscar" in p
                assert "cif" in p
                assert p["poscar"].exists()
                assert p["cif"].exists()

    def test_save_structure_creates_work_dir(self, simple_slab):
        """save_structure should create work directory if needed."""
        from nh3sofc import save_structure

        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = f"{tmpdir}/new_subdir/configs"
            paths = save_structure(simple_slab, work_dir, "test_struct")

            assert "work_dir" in paths
            assert paths["poscar"].exists()
            assert paths["poscar"].parent == Path(work_dir)

    def test_save_structure_auto_generates_work_dir(self, simple_slab):
        """save_structure should auto-generate meaningful work_dir when not provided."""
        from nh3sofc import save_structure
        import os

        # Change to temp directory to avoid polluting cwd
        original_cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            try:
                paths = save_structure(simple_slab, name="test_struct")

                # Should have auto-generated work_dir
                assert "work_dir" in paths
                assert paths["work_dir"].exists()
                assert paths["poscar"].exists()

                # Work dir name should contain formula info
                work_dir_name = paths["work_dir"].name
                assert "atoms" in work_dir_name  # Should have atom count
            finally:
                os.chdir(original_cwd)

    def test_save_configurations_auto_generates_work_dir(self, simple_slab):
        """save_configurations should auto-generate work_dir when not provided."""
        from nh3sofc import save_configurations
        import os

        configs = [simple_slab, simple_slab.copy()]

        original_cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            try:
                result = save_configurations(configs, name_prefix="config")

                # Should have auto-generated work_dir
                assert "work_dir" in result
                assert result["work_dir"].exists()
                assert len(result["configs"]) == 2

                # All configs should be saved
                for p in result["configs"]:
                    assert p["poscar"].exists()
            finally:
                os.chdir(original_cwd)


class TestDefectPlacementProbabilistic:
    """
    Bug: DefectBuilder placement strategies were deterministic.
    "surface" placement always selected the topmost atoms, so all
    configurations with same strategy were identical.

    Fix: Changed to probability-weighted selection where surface/near_N
    atoms have higher probability but selection is still stochastic.
    Added surface_n_preference and vacancy_preference parameters.

    Date fixed: 2026-01-22
    """

    def test_surface_placement_produces_different_configs(self, perovskite_bulk):
        """Surface placement should produce different configurations."""
        from nh3sofc.structure import SurfaceBuilder, DefectBuilder

        builder = SurfaceBuilder(perovskite_bulk)
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=6,
            vacuum=15.0,
        )

        defect = DefectBuilder(slab)

        # Generate multiple configs with same strategy
        configs = []
        for i in range(5):
            atoms = defect.create_oxynitride(
                nitrogen_fraction=0.5,
                placement="surface",
                surface_n_preference=0.7,
                random_seed=i,
            )
            configs.append(atoms)

        # Configs should have different N positions
        n_positions_list = []
        for atoms in configs:
            n_indices = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "N"]
            n_positions_list.append(tuple(sorted(n_indices)))

        unique_configs = set(n_positions_list)
        assert len(unique_configs) > 1, \
            "Surface placement should produce different configurations"

    def test_preference_affects_n_distribution(self, perovskite_bulk):
        """Higher surface_n_preference should put more N at surface."""
        from nh3sofc.structure import SurfaceBuilder, DefectBuilder

        builder = SurfaceBuilder(perovskite_bulk)
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=6,
            vacuum=15.0,
        )

        defect = DefectBuilder(slab)
        z_positions = slab.atoms.positions[:, 2]
        z_cutoff = z_positions.max() - 0.3 * (z_positions.max() - z_positions.min())

        # Low preference (nearly random)
        low_pref_surface_n = 0
        for i in range(10):
            atoms = defect.create_oxynitride(
                nitrogen_fraction=0.5,
                placement="surface",
                surface_n_preference=0.55,  # Nearly random
                random_seed=i,
            )
            n_indices = [j for j, s in enumerate(atoms.get_chemical_symbols()) if s == "N"]
            surface_n = sum(1 for idx in n_indices if atoms.positions[idx, 2] > z_cutoff)
            low_pref_surface_n += surface_n

        # High preference
        high_pref_surface_n = 0
        for i in range(10):
            atoms = defect.create_oxynitride(
                nitrogen_fraction=0.5,
                placement="surface",
                surface_n_preference=0.95,  # Strongly prefer surface
                random_seed=i + 100,
            )
            n_indices = [j for j, s in enumerate(atoms.get_chemical_symbols()) if s == "N"]
            surface_n = sum(1 for idx in n_indices if atoms.positions[idx, 2] > z_cutoff)
            high_pref_surface_n += surface_n

        # High preference should have more N at surface on average
        assert high_pref_surface_n > low_pref_surface_n, \
            f"High preference ({high_pref_surface_n}) should have more surface N than low ({low_pref_surface_n})"

    def test_analyze_defect_distribution(self, perovskite_bulk):
        """analyze_defect_distribution should return expected statistics."""
        from nh3sofc.structure import SurfaceBuilder, DefectBuilder
        from nh3sofc.structure import analyze_defect_distribution

        builder = SurfaceBuilder(perovskite_bulk)
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=6,
            vacuum=15.0,
        )

        defect = DefectBuilder(slab)
        oxynitride = defect.create_oxynitride(
            nitrogen_fraction=0.5,
            placement="surface",
            surface_n_preference=0.8,
            random_seed=42,
        )

        stats = analyze_defect_distribution(oxynitride, z_threshold=0.3)

        # Should have required keys
        assert "n_total" in stats
        assert "n_surface" in stats
        assert "surface_n_ratio" in stats
        assert "bulk_n_ratio" in stats

        # Surface N ratio should be higher than bulk for surface placement
        # (with high preference)
        assert stats["surface_n_ratio"] >= stats["bulk_n_ratio"], \
            f"Surface N ratio ({stats['surface_n_ratio']}) should be >= bulk ({stats['bulk_n_ratio']})"

    def test_analyze_oxynitride_pool(self, perovskite_bulk):
        """analyze_oxynitride_pool should return statistics by strategy."""
        from nh3sofc.structure import SurfaceBuilder, DefectBuilder
        from nh3sofc.structure import analyze_oxynitride_pool

        builder = SurfaceBuilder(perovskite_bulk)
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=6,
            vacuum=15.0,
        )

        defect = DefectBuilder(slab)
        pool = defect.create_oxynitride_pool(
            nitrogen_fraction=0.5,
            n_configs_per_strategy=3,
            strategies=["random", "surface"],
            surface_n_preference=0.8,
            random_seed=42,
        )

        stats = analyze_oxynitride_pool(pool, z_threshold=0.3)

        # Should have stats by strategy
        assert "by_strategy" in stats
        assert "random" in stats["by_strategy"]
        assert "surface" in stats["by_strategy"]

        # Surface strategy should have higher surface N ratio than random
        surface_ratio = stats["by_strategy"]["surface"]["surface_n_ratio_mean"]
        random_ratio = stats["by_strategy"]["random"]["surface_n_ratio_mean"]
        assert surface_ratio >= random_ratio, \
            f"Surface strategy ({surface_ratio:.2f}) should have >= N ratio than random ({random_ratio:.2f})"


class TestAdsorbateRotationConsistency:
    """
    Bug: Rotation handling was inconsistent across placement methods:
    - add_with_collision() was missing x-rotation (only z and y)
    - add_catkit() had no rotation support at all
    - Euler angle sampling was biased (not uniform on SO(3))

    Fix: Created _apply_random_rotation() helper with uniform SO(3) sampling
    and applied it consistently across all methods.

    Date fixed: 2026-02-02
    """

    def test_all_methods_produce_varied_orientations(self, simple_slab):
        """All placement methods should produce varied molecular orientations."""
        from nh3sofc.structure import AdsorbatePlacer

        placer = AdsorbatePlacer(simple_slab)

        # Test add_random
        configs_random = placer.add_random("NH3", n_configs=5, random_seed=42)
        assert len(configs_random) == 5

        # Check that orientations are different (not all H atoms at same positions)
        h_positions = []
        for config in configs_random:
            h_atoms = [a for a in config if a.symbol == "H"]
            h_positions.append([a.position.copy() for a in h_atoms])

        # At least some configs should have different H positions
        all_same = all(
            np.allclose(h_positions[0][0], h_positions[i][0], atol=0.1)
            for i in range(1, len(h_positions))
        )
        assert not all_same, "Random placement should produce varied orientations"

    def test_add_with_collision_has_full_rotation(self, simple_slab):
        """add_with_collision should apply full 3D rotation (z, y, x)."""
        from nh3sofc.structure import AdsorbatePlacer

        placer = AdsorbatePlacer(simple_slab)

        # Generate multiple configs and check orientations vary in all dimensions
        configs = placer.add_with_collision(
            "NH3", n_configs=10, min_distance=1.5, random_seed=42
        )

        assert len(configs) > 0, "Should generate at least one config"

        # Get NH3 positions (last 4 atoms: N + 3H)
        nh3_positions = [c.positions[-4:] for c in configs]

        # Check variance in all three dimensions
        all_positions = np.array(nh3_positions)
        variance = np.var(all_positions, axis=0).mean(axis=0)

        # Should have variance in all xyz (not just xy)
        assert variance[2] > 0.01, \
            f"Z variance too low ({variance[2]:.4f}), x-rotation may be missing"

    def test_helper_method_applies_three_rotations(self, simple_slab):
        """_apply_random_rotation should apply rotations around all 3 axes."""
        from nh3sofc.structure import AdsorbatePlacer
        from ase.build import molecule

        placer = AdsorbatePlacer(simple_slab)

        np.random.seed(42)
        orientations = []

        for _ in range(20):
            mol = molecule("NH3")
            original_positions = mol.positions.copy()
            placer._apply_random_rotation(mol)
            orientations.append(mol.positions - mol.positions.mean(axis=0))

        orientations = np.array(orientations)

        # Check variance exists in all dimensions
        variance = np.var(orientations, axis=0).mean(axis=0)
        assert all(v > 0.01 for v in variance), \
            f"Rotation should vary in all dimensions, got variance: {variance}"


class TestAddOnSiteSurfaceValidation:
    """
    Bug: add_on_site() with site_indices would place adsorbates on any
    atom index, including subsurface atoms, bypassing surface detection.

    Fix: site_indices are now validated against surface atoms. Non-surface
    indices are filtered out with a warning.

    Date fixed: 2026-02-02
    """

    def test_site_indices_filtered_to_surface(self, perovskite_bulk):
        """site_indices should only use atoms that are on the surface."""
        from nh3sofc.structure import SurfaceBuilder, AdsorbatePlacer
        from nh3sofc.core.base import get_surface_atoms

        # Create a slab
        builder = SurfaceBuilder(perovskite_bulk)
        slab = builder.create_surface(miller_index=(0, 0, 1), layers=4, vacuum=15.0)

        placer = AdsorbatePlacer(slab.atoms)

        # Get actual surface atoms
        surface_indices, z_cutoff = get_surface_atoms(slab.atoms, z_threshold=0.2)

        # Try to place on a mix of surface and non-surface indices
        # Index 0 is likely in the bulk (bottom of slab)
        mixed_indices = [0, 1] + surface_indices[:2]

        configs = placer.add_on_site(
            "NH3",
            site_indices=mixed_indices,
            n_orientations=1,
            random_seed=42,
        )

        # Should only get configs for the surface atoms in the list
        n_surface_in_mixed = len([i for i in mixed_indices if i in surface_indices])
        assert len(configs) == n_surface_in_mixed, \
            f"Expected {n_surface_in_mixed} configs (surface atoms only), got {len(configs)}"

    def test_add_on_site_only_uses_surface_atoms(self, perovskite_bulk):
        """add_on_site with atom_types should only place on surface atoms of that type."""
        from nh3sofc.structure import SurfaceBuilder, AdsorbatePlacer
        from nh3sofc.core.base import get_surface_atoms

        builder = SurfaceBuilder(perovskite_bulk)
        slab = builder.create_surface(miller_index=(0, 0, 1), layers=4, vacuum=15.0)

        placer = AdsorbatePlacer(slab.atoms)

        # Get surface atoms
        surface_indices, z_cutoff = get_surface_atoms(slab.atoms, z_threshold=0.2)

        # Place on La atoms only
        configs = placer.add_on_site(
            "NH3",
            atom_types=["La"],
            n_orientations=1,
            random_seed=42,
        )

        # Count La atoms on surface
        n_la_surface = sum(1 for i in surface_indices if slab.atoms[i].symbol == "La")

        assert len(configs) == n_la_surface, \
            f"Expected {n_la_surface} configs (La on surface), got {len(configs)}"

        # Verify adsorbate z positions are above surface
        for config in configs:
            nh3_z = config.positions[-4:, 2].min()  # NH3 is last 4 atoms
            assert nh3_z > z_cutoff, \
                f"Adsorbate placed below surface: z={nh3_z:.2f} < cutoff={z_cutoff:.2f}"
