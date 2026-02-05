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

    def test_import_surface_builders(self):
        from nh3sofc.structure import SurfaceBuilder, SlabStructure

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


class TestDopantModule:
    """Test dopant module for doped ceria (GDC, SDC, PDC)."""

    def test_import_dopant_classes(self):
        from nh3sofc.structure import (
            DopantBuilder,
            DopedCeriaStructure,
            analyze_dopant_distribution,
            generate_dopant_series,
            print_dopant_analysis,
        )

    def test_import_dopant_constants(self):
        from nh3sofc.core.constants import (
            ACCEPTOR_DOPANTS,
            HOST_CATIONS,
            CHARGE_COMPENSATION,
        )
        assert "Gd" in ACCEPTOR_DOPANTS
        assert "Sm" in ACCEPTOR_DOPANTS
        assert "Pr" in ACCEPTOR_DOPANTS
        assert "Ce" in HOST_CATIONS
        assert 3 in CHARGE_COMPENSATION

    def test_dopant_builder_instantiation(self, ceo2_slab):
        from nh3sofc.structure import DopantBuilder
        builder = DopantBuilder(ceo2_slab)
        assert builder.atoms is not None
        assert len(builder.get_element_indices("Ce")) == 8
        assert len(builder.get_element_indices("O")) == 16

    def test_create_gdc(self, ceo2_slab):
        """Test creating Gd-doped ceria (GDC)."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(ceo2_slab)
        original_ce = len(builder.get_element_indices("Ce"))
        original_o = len(builder.get_element_indices("O"))

        # 25% Gd doping (2 out of 8 Ce atoms)
        gdc = builder.create_doped_structure(
            dopant="Gd",
            dopant_fraction=0.25,
            random_seed=42,
        )

        # Check dopant substitution
        gd_count = sum(1 for s in gdc.get_chemical_symbols() if s == "Gd")
        ce_count = sum(1 for s in gdc.get_chemical_symbols() if s == "Ce")
        o_count = sum(1 for s in gdc.get_chemical_symbols() if s == "O")

        assert gd_count == 2  # 25% of 8 Ce
        assert ce_count == 6  # 8 - 2
        # Check charge compensation: 2 Gd -> 1 vacancy
        assert o_count == original_o - 1  # 16 - 1 = 15

    def test_create_sdc(self, ceo2_slab):
        """Test creating Sm-doped ceria (SDC)."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(ceo2_slab)

        sdc = builder.create_doped_structure(
            dopant="Sm",
            dopant_fraction=0.50,  # 50% = 4 dopants
            random_seed=42,
        )

        sm_count = sum(1 for s in sdc.get_chemical_symbols() if s == "Sm")
        o_count = sum(1 for s in sdc.get_chemical_symbols() if s == "O")

        assert sm_count == 4  # 50% of 8 Ce
        # 4 Sm -> 2 vacancies
        assert o_count == 14  # 16 - 2

    def test_create_doped_pool(self, ceo2_slab):
        """Test generating pool of doped configurations."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(ceo2_slab)
        pool = builder.create_doped_pool(
            dopant="Gd",
            dopant_fraction=0.25,
            n_configs=2,
            strategies=["random", "surface"],
            random_seed=42,
        )

        assert len(pool) == 4  # 2 strategies * 2 configs
        assert all("atoms" in config for config in pool)
        assert all("vacancy_placement" in config for config in pool)

    def test_create_doped_pool_with_fraction_list(self, ceo2_slab):
        """Test generating pool with multiple dopant fractions."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(ceo2_slab)
        pool = builder.create_doped_pool(
            dopant="Gd",
            dopant_fraction=[0.125, 0.25, 0.375],  # 1, 2, 3 dopants
            n_configs=2,
            random_seed=42,
        )

        # 3 fractions * 1 strategy (default "random" for list) * 2 configs = 6
        assert len(pool) == 6

        # Check that different fractions are present
        fractions = set(c["dopant_fraction"] for c in pool)
        assert fractions == {0.125, 0.25, 0.375}

    def test_doped_ceria_structure_from_ceria(self, ceo2_slab):
        """Test DopedCeriaStructure.from_ceria factory method."""
        from nh3sofc.structure import DopedCeriaStructure

        doped = DopedCeriaStructure.from_ceria(
            ceo2_slab,
            dopant="Gd",
            dopant_fraction=0.25,
            random_seed=42,
        )

        assert doped.dopant == "Gd"
        assert doped.dopant_fraction == 0.25
        assert doped.n_vacancies == 1  # 2 Gd -> 1 vacancy
        assert doped.get_dopant_name() == "GDC"

    def test_doped_ceria_stoichiometry(self, ceo2_slab):
        """Test stoichiometry calculation for doped ceria."""
        from nh3sofc.structure import DopedCeriaStructure

        doped = DopedCeriaStructure.from_ceria(
            ceo2_slab,
            dopant="Sm",
            dopant_fraction=0.50,
            random_seed=42,
        )

        stoich = doped.get_stoichiometry()

        # 4 Sm, 4 Ce, 14 O -> total cations = 8
        assert abs(stoich["Sm"] - 0.5) < 0.01  # 4/8
        assert abs(stoich["Ce"] - 0.5) < 0.01  # 4/8
        assert abs(stoich["O"] - 1.75) < 0.01  # 14/8

    def test_analyze_dopant_distribution(self, ceo2_slab):
        """Test dopant distribution analysis."""
        from nh3sofc.structure import DopantBuilder, analyze_dopant_distribution

        builder = DopantBuilder(ceo2_slab)
        gdc = builder.create_doped_structure(
            dopant="Gd",
            dopant_fraction=0.25,
            random_seed=42,
        )

        stats = analyze_dopant_distribution(
            gdc,
            dopant="Gd",
            reference_atoms=ceo2_slab,
        )

        assert stats["dopant"] == "Gd"
        assert stats["dopant_total"] == 2
        assert stats["vacancy_total"] == 1
        assert "dopant_surface_fraction" in stats
        assert "vacancy_near_dopant_fraction" in stats

    def test_generate_dopant_series_deprecated(self, ceo2_slab):
        """Test deprecated generate_dopant_series (backwards compatibility)."""
        import warnings
        from nh3sofc.structure import generate_dopant_series

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            series = generate_dopant_series(
                ceo2_slab,
                dopant="Gd",
                dopant_fractions=[0.125, 0.25],  # 1 and 2 dopants
                n_configs=2,
                random_seed=42,
            )
            # Check deprecation warning was issued
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "deprecated" in str(w[0].message).lower()

        assert len(series) == 4  # 2 fractions * 2 configs
        assert all("atoms" in config for config in series)
        assert all("dopant_fraction" in config for config in series)

    def test_pr_mixed_valence(self, ceo2_slab):
        """Test Pr doping with mixed Pr³⁺/Pr⁴⁺."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(ceo2_slab)
        original_o = len(builder.get_element_indices("O"))

        # 50% Pr with only half being Pr³⁺
        pdc = builder.create_doped_structure(
            dopant="Pr",
            dopant_fraction=0.50,  # 4 Pr
            pr_trivalent_fraction=0.5,  # Only 2 effective for vacancies
            random_seed=42,
        )

        pr_count = sum(1 for s in pdc.get_chemical_symbols() if s == "Pr")
        o_count = sum(1 for s in pdc.get_chemical_symbols() if s == "O")

        assert pr_count == 4
        # Only 2 effective Pr³⁺ -> 1 vacancy
        assert o_count == original_o - 1

    def test_tb_mixed_valence(self, ceo2_slab):
        """Test Tb doping with mixed Tb³⁺/Tb⁴⁺."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(ceo2_slab)
        original_o = len(builder.get_element_indices("O"))

        # 50% Tb with only half being Tb³⁺
        tdc = builder.create_doped_structure(
            dopant="Tb",
            dopant_fraction=0.50,  # 4 Tb
            pr_trivalent_fraction=0.5,  # Only 2 effective for vacancies
            random_seed=42,
        )

        tb_count = sum(1 for s in tdc.get_chemical_symbols() if s == "Tb")
        o_count = sum(1 for s in tdc.get_chemical_symbols() if s == "O")

        assert tb_count == 4
        # Only 2 effective Tb³⁺ -> 1 vacancy
        assert o_count == original_o - 1

    def test_vacancy_placement_strategies(self, ceo2_slab):
        """Test different vacancy placement strategies."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(ceo2_slab)

        for strategy in ["random", "surface", "bulk", "near_dopant"]:
            gdc = builder.create_doped_structure(
                dopant="Gd",
                dopant_fraction=0.50,  # Need enough for meaningful test
                vacancy_placement=strategy,
                vacancy_preference=0.8,
                random_seed=42,
            )
            # Should complete without error
            assert gdc is not None
            gd_count = sum(1 for s in gdc.get_chemical_symbols() if s == "Gd")
            assert gd_count == 4


class TestDopedFluoriteGeneralization:
    """Test generalized dopant module for multiple fluorite materials (YSZ, ScSZ, etc.)."""

    def test_import_doped_fluorite_structure(self):
        """Test that DopedFluoriteStructure can be imported."""
        from nh3sofc.structure import DopedFluoriteStructure
        assert DopedFluoriteStructure is not None

    def test_import_new_constants(self):
        """Test that new constants are available."""
        from nh3sofc.core.constants import (
            ACCEPTOR_DOPANTS,
            HOST_CATIONS,
            FLUORITE_MATERIAL_NAMES,
        )
        # New trivalent dopant
        assert "Sc" in ACCEPTOR_DOPANTS
        assert ACCEPTOR_DOPANTS["Sc"]["charge"] == +3
        # New divalent dopants
        assert "Ca" in ACCEPTOR_DOPANTS
        assert ACCEPTOR_DOPANTS["Ca"]["charge"] == +2
        assert "Mg" in ACCEPTOR_DOPANTS
        assert ACCEPTOR_DOPANTS["Mg"]["charge"] == +2
        # New host cations
        assert "Hf" in HOST_CATIONS
        assert "Th" in HOST_CATIONS
        # Material names mapping
        assert ("Zr", "Y") in FLUORITE_MATERIAL_NAMES
        assert FLUORITE_MATERIAL_NAMES[("Zr", "Y")] == "YSZ"
        assert FLUORITE_MATERIAL_NAMES[("Zr", "Sc")] == "ScSZ"

    def test_create_ysz(self, zro2_slab):
        """Test creating Y-doped zirconia (YSZ)."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(zro2_slab)
        original_zr = len(builder.get_element_indices("Zr"))
        original_o = len(builder.get_element_indices("O"))

        # 25% Y doping (2 out of 8 Zr atoms)
        ysz = builder.create_doped_structure(
            dopant="Y",
            dopant_fraction=0.25,
            host_cation="Zr",
            random_seed=42,
        )

        # Check dopant substitution
        y_count = sum(1 for s in ysz.get_chemical_symbols() if s == "Y")
        zr_count = sum(1 for s in ysz.get_chemical_symbols() if s == "Zr")
        o_count = sum(1 for s in ysz.get_chemical_symbols() if s == "O")

        assert y_count == 2  # 25% of 8 Zr
        assert zr_count == 6  # 8 - 2
        # Check charge compensation: 2 Y -> 1 vacancy
        assert o_count == original_o - 1  # 16 - 1 = 15

    def test_create_scsz(self, zro2_slab):
        """Test creating Sc-doped zirconia (ScSZ)."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(zro2_slab)

        scsz = builder.create_doped_structure(
            dopant="Sc",
            dopant_fraction=0.50,  # 50% = 4 dopants
            host_cation="Zr",
            random_seed=42,
        )

        sc_count = sum(1 for s in scsz.get_chemical_symbols() if s == "Sc")
        o_count = sum(1 for s in scsz.get_chemical_symbols() if s == "O")

        assert sc_count == 4  # 50% of 8 Zr
        # 4 Sc -> 2 vacancies (trivalent)
        assert o_count == 14  # 16 - 2

    def test_divalent_dopant_ca(self, zro2_slab):
        """Test divalent Ca dopant (CSZ - 1:1 dopant:vacancy ratio)."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(zro2_slab)
        original_o = len(builder.get_element_indices("O"))

        # 25% Ca doping (2 out of 8 Zr atoms)
        csz = builder.create_doped_structure(
            dopant="Ca",
            dopant_fraction=0.25,
            host_cation="Zr",
            random_seed=42,
        )

        ca_count = sum(1 for s in csz.get_chemical_symbols() if s == "Ca")
        o_count = sum(1 for s in csz.get_chemical_symbols() if s == "O")

        assert ca_count == 2
        # Divalent: 1 Ca -> 1 vacancy (1:1 ratio)
        assert o_count == original_o - 2  # 16 - 2 = 14

    def test_divalent_dopant_mg(self, zro2_slab):
        """Test divalent Mg dopant (MSZ)."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(zro2_slab)
        original_o = len(builder.get_element_indices("O"))

        msz = builder.create_doped_structure(
            dopant="Mg",
            dopant_fraction=0.125,  # 1 out of 8 Zr atoms
            host_cation="Zr",
            random_seed=42,
        )

        mg_count = sum(1 for s in msz.get_chemical_symbols() if s == "Mg")
        o_count = sum(1 for s in msz.get_chemical_symbols() if s == "O")

        assert mg_count == 1
        # Divalent: 1 Mg -> 1 vacancy
        assert o_count == original_o - 1  # 16 - 1 = 15

    def test_from_zirconia_factory(self, zro2_slab):
        """Test DopedFluoriteStructure.from_zirconia factory method."""
        from nh3sofc.structure import DopedFluoriteStructure

        ysz = DopedFluoriteStructure.from_zirconia(
            zro2_slab,
            dopant="Y",
            dopant_fraction=0.25,
            random_seed=42,
        )

        assert ysz.dopant == "Y"
        assert ysz.dopant_fraction == 0.25
        assert ysz.host_cation == "Zr"
        assert ysz.n_vacancies == 1  # 2 Y -> 1 vacancy
        assert ysz.get_material_name() == "YSZ"

    def test_from_parent_factory(self, zro2_slab):
        """Test DopedFluoriteStructure.from_parent general factory."""
        from nh3sofc.structure import DopedFluoriteStructure

        scsz = DopedFluoriteStructure.from_parent(
            zro2_slab,
            dopant="Sc",
            dopant_fraction=0.25,
            host_cation="Zr",
            random_seed=42,
        )

        assert scsz.dopant == "Sc"
        assert scsz.host_cation == "Zr"
        assert scsz.get_material_name() == "ScSZ"

    def test_doped_fluorite_from_ceria(self, ceo2_slab):
        """Test DopedFluoriteStructure.from_ceria convenience method."""
        from nh3sofc.structure import DopedFluoriteStructure

        gdc = DopedFluoriteStructure.from_ceria(
            ceo2_slab,
            dopant="Gd",
            dopant_fraction=0.25,
            random_seed=42,
        )

        assert gdc.dopant == "Gd"
        assert gdc.host_cation == "Ce"
        assert gdc.get_material_name() == "GDC"

    def test_get_material_name_fallback(self, ceo2_slab):
        """Test material name fallback for unknown combinations."""
        from nh3sofc.structure import DopedFluoriteStructure

        # Sc-doped ceria is uncommon but should work
        doped = DopedFluoriteStructure.from_parent(
            ceo2_slab,
            dopant="Sc",
            dopant_fraction=0.125,
            host_cation="Ce",
            random_seed=42,
        )

        # Should fall back to descriptive name
        name = doped.get_material_name()
        assert "Sc" in name and "CeO2" in name

    def test_doped_ceria_structure_deprecated(self, ceo2_slab):
        """Test that DopedCeriaStructure shows deprecation warning."""
        import warnings
        from nh3sofc.structure import DopedCeriaStructure
        from ase.build import bulk

        atoms = ceo2_slab.copy()

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            # This should trigger deprecation warning
            doped = DopedCeriaStructure(
                atoms=atoms,
                dopant="Gd",
                dopant_fraction=0.25,
                n_vacancies=1,
            )
            # Check deprecation warning was issued
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "deprecated" in str(w[0].message).lower()

    def test_get_dopant_name_deprecated(self, ceo2_slab):
        """Test that get_dopant_name() shows deprecation warning."""
        import warnings
        from nh3sofc.structure import DopedFluoriteStructure

        gdc = DopedFluoriteStructure.from_ceria(
            ceo2_slab,
            dopant="Gd",
            dopant_fraction=0.25,
            random_seed=42,
        )

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            name = gdc.get_dopant_name()
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert name == "GDC"

    def test_backwards_compat_from_ceria(self, ceo2_slab):
        """Test backwards compatibility: DopedCeriaStructure.from_ceria still works."""
        import warnings
        from nh3sofc.structure import DopedCeriaStructure

        # from_ceria should work without deprecation warning (it's on DopedFluoriteStructure)
        doped = DopedCeriaStructure.from_ceria(
            ceo2_slab,
            dopant="Sm",
            dopant_fraction=0.25,
            random_seed=42,
        )

        assert doped.dopant == "Sm"
        assert doped.host_cation == "Ce"
        assert doped.get_material_name() == "SDC"

    def test_create_doped_pool_zirconia(self, zro2_slab):
        """Test creating pool of YSZ configurations."""
        from nh3sofc.structure import DopantBuilder

        builder = DopantBuilder(zro2_slab)
        pool = builder.create_doped_pool(
            dopant="Y",
            dopant_fraction=0.25,
            n_configs=2,
            strategies=["random", "surface"],
            host_cation="Zr",
            random_seed=42,
        )

        assert len(pool) == 4  # 2 strategies * 2 configs
        for config in pool:
            y_count = sum(1 for s in config["atoms"].get_chemical_symbols() if s == "Y")
            assert y_count == 2  # 25% of 8 Zr
