"""Smoke tests - verify basic imports and instantiation work."""

import pytest


class TestImports:
    """Test that all main modules can be imported."""

    def test_import_structure_module(self):
        from nh3sofc.structure import (
            BulkStructure,
            SurfaceBuilder,
            SlabStructure,
            AdsorbatePlacer,
            DefectBuilder,
        )

    def test_import_specialized_builders(self):
        from nh3sofc.structure import PerovskiteSurfaceBuilder

    def test_import_calculators(self):
        from nh3sofc.calculators.vasp import VASPInputGenerator

    def test_import_constants(self):
        from nh3sofc.core.constants import (
            DEFAULT_FORMAL_CHARGES,
            PEROVSKITE_SITES,
            DEFAULT_VASP_PARAMS,
        )


class TestBasicInstantiation:
    """Test that main classes can be instantiated."""

    def test_surface_builder(self, rocksalt_bulk):
        from nh3sofc.structure import SurfaceBuilder
        builder = SurfaceBuilder(rocksalt_bulk)
        assert builder.bulk is not None

    def test_perovskite_builder(self, perovskite_bulk):
        from nh3sofc.structure import PerovskiteSurfaceBuilder
        builder = PerovskiteSurfaceBuilder(perovskite_bulk)
        # Should auto-detect La as A-site, V as B-site
        assert builder.A_site == 'La'
        assert builder.B_site == 'V'

    def test_adsorbate_placer(self, simple_slab):
        from nh3sofc.structure import AdsorbatePlacer
        placer = AdsorbatePlacer(simple_slab)
        assert placer.slab is not None

    def test_vasp_input_generator(self, simple_slab):
        from nh3sofc.calculators.vasp import VASPInputGenerator
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            vasp = VASPInputGenerator(simple_slab, calc_type='relax', work_dir=tmpdir)
            assert vasp.atoms is not None


class TestBasicSurfaceCreation:
    """Test basic surface creation works."""

    def test_create_surface(self, rocksalt_bulk):
        from nh3sofc.structure import SurfaceBuilder
        builder = SurfaceBuilder(rocksalt_bulk)
        slab = builder.create_surface(
            miller_index=(0, 0, 1),
            layers=4,
            vacuum=10.0
        )
        assert len(slab.atoms) > 0
        assert slab.miller_index == (0, 0, 1)

    def test_create_symmetric_slab(self, rocksalt_bulk):
        from nh3sofc.structure import SurfaceBuilder
        builder = SurfaceBuilder(rocksalt_bulk)
        slab = builder.create_symmetric_slab(
            miller_index=(0, 0, 1),
            layers=5,
            vacuum=10.0
        )
        assert len(slab.atoms) > 0
        assert slab.termination == 'symmetric'


class TestPolarityAnalysis:
    """Test polarity analysis methods."""

    def test_check_polarity(self, rocksalt_bulk):
        from nh3sofc.structure import SurfaceBuilder
        builder = SurfaceBuilder(rocksalt_bulk)
        slab = builder.create_surface(miller_index=(0, 0, 1), layers=4)

        polarity = slab.check_polarity()

        assert 'is_polar' in polarity
        assert 'dipole_z' in polarity
        assert 'recommendation' in polarity

    def test_identify_layers(self, rocksalt_bulk):
        from nh3sofc.structure import SurfaceBuilder
        builder = SurfaceBuilder(rocksalt_bulk)
        slab = builder.create_surface(miller_index=(0, 0, 1), layers=4)

        layers = slab.identify_layers()

        assert len(layers) > 0
        assert 'z' in layers[0]
        assert 'composition' in layers[0]


class TestAdsorbatePlacement:
    """Test adsorbate placement methods."""

    def test_add_simple(self, simple_slab):
        from nh3sofc.structure import AdsorbatePlacer
        placer = AdsorbatePlacer(simple_slab)

        original_count = len(simple_slab)
        result = placer.add_simple('NH3', position=(0.5, 0.5), height=2.0)

        # NH3 has 4 atoms
        assert len(result) == original_count + 4

    def test_add_random(self, simple_slab):
        from nh3sofc.structure import AdsorbatePlacer
        placer = AdsorbatePlacer(simple_slab)

        configs = placer.add_random('NH3', n_configs=3, random_seed=42)

        assert len(configs) == 3
        for config in configs:
            assert len(config) > len(simple_slab)
