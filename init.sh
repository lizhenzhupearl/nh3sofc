#!/bin/bash
set -e

echo "=== NH3SOFC Harness Initialization ==="
echo ""

# 1. Verify we're in the right directory
if [ ! -f pyproject.toml ]; then
  echo "ERROR: pyproject.toml not found. Are you in the NH3SOFC repo root?"
  exit 1
fi

# 2. Install package in editable mode
echo "=== Installing dependencies ==="
pip install -e . --quiet 2>&1 | tail -1
echo "Install complete."
echo ""

# 3. Run tests
echo "=== Running tests ==="
pytest tests/ -v --tb=short
echo ""

# 4. Import smoke check (mirrors CI lint job)
echo "=== Import smoke check ==="
python -c "from nh3sofc.structure import SurfaceBuilder, AdsorbatePlacer; print('  SurfaceBuilder, AdsorbatePlacer: OK')"
python -c "from nh3sofc.calculators.vasp import VASPInputGenerator; print('  VASPInputGenerator: OK')"
echo ""

echo "=== Verification Complete ==="
echo ""
echo "Next steps:"
echo "1. Read feature-list.json to see current feature state"
echo "2. Read progress.md for session context"
echo "3. Pick ONE unfinished feature to work on"
echo "4. Implement only that feature"
echo "5. Re-run ./init.sh before claiming done"
