# Key Directories

```
nh3sofc/
├── structure/      # Bulk, surface, defects, adsorbates
├── calculators/    # VASP, MACE interfaces
├── analysis/       # Energetics, thermochemistry, kinetics
├── workflows/      # Automated calculation workflows
└── database/       # Naming conventions, ASE DB

tests/
├── conftest.py     # Shared fixtures
├── test_smoke.py   # Import and basic tests
└── test_regressions.py  # Tests for fixed bugs

docs/
├── plans/          # Implementation plans
├── tutorials/      # Step-by-step guides
├── api/            # API reference
└── claude/         # Claude Code reference docs

.github/workflows/
├── ci.yml          # Run tests on push/PR
└── docs.yml        # Deploy documentation
```

## Verification Commands

```bash
# Full verification (recommended — same as init.sh)
./init.sh

# Individual checks
pytest tests/ -v
pytest tests/test_smoke.py
pytest tests/test_regressions.py
python -c "from nh3sofc.structure import SurfaceBuilder, AdsorbatePlacer"
python -c "from nh3sofc.calculators.vasp import VASPInputGenerator"
```
