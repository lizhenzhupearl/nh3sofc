# Claude Code Instructions for NH3SOFC

## Project Overview

NH3SOFC is a Python package for computational catalysis research focused on ammonia decomposition for solid oxide fuel cells. It provides structure generation, VASP/MACE calculations, and analysis tools.

## Plan Management

### IMPORTANT: Save All Plans

Every implementation plan MUST be saved to `docs/plans/` with proper naming:

1. **Filename format:** `YYYY-MM-DD-descriptive-title.md`
   - Example: `2026-01-22-fix-adsorbate-placer-api.md`

2. **When to save:**
   - After a plan is approved and before starting implementation
   - Update the plan with "Completed" status and commit references after finishing

3. **Plan template:**
   ```markdown
   # Plan: [Descriptive Title]

   **Date:** YYYY-MM-DD
   **Status:** Draft | Approved | Completed

   ## Summary
   [Brief description]

   ## Problem / Motivation
   [Why this change is needed]

   ## Proposed Solution
   [Implementation approach]

   ## Files to Modify
   - `path/to/file.py` - Description

   ## Verification
   [How to test]

   ## Commits
   [Add after completion]
   ```

4. **Update the index:** Add entry to `docs/plans/README.md` index table

### Why This Matters

Plans represent significant intellectual work and decision-making. They serve as:
- Decision records
- Methodology documentation
- Reproducibility artifacts
- Project history

## Code Style

- Follow existing patterns in the codebase
- Use type hints
- Include docstrings with Examples section
- Prefer editing existing files over creating new ones

## Testing

### Run Tests Before and After Changes

The `tests/` directory contains automated tests. **Always run tests** to verify changes work:

```bash
# Run all tests
pytest tests/

# Run with verbose output
pytest tests/ -v

# Run specific test file
pytest tests/test_regressions.py
```

### Test Structure

| File | Purpose |
|------|---------|
| `tests/conftest.py` | Shared fixtures (mock structures) |
| `tests/test_smoke.py` | Basic imports and instantiation |
| `tests/test_regressions.py` | Tests for bugs that were fixed |

### When Developing

1. **After making changes:** Run `pytest tests/` to catch regressions
2. **When fixing a bug:** Add a regression test to `tests/test_regressions.py`
3. **When adding features:** Add smoke tests for new classes/methods

### Adding Regression Tests

When fixing a bug, add a test to `tests/test_regressions.py`:

```python
class TestYourBugFix:
    """
    Bug: [Description of the bug]
    Fix: [How it was fixed]
    Date fixed: YYYY-MM-DD
    """

    def test_bug_is_fixed(self, fixture):
        # Test that verifies the fix works
        assert result == expected
```

### CI/CD

GitHub Actions runs tests automatically on push/PR to `main`. See `.github/workflows/ci.yml`.

## Documentation

- Update relevant docs when changing APIs
- Keep `docs/index.md` feature lists current
- Add tutorials for significant new features

## Key Directories

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
├── plans/          # Implementation plans (SAVE PLANS HERE!)
├── tutorials/      # Step-by-step guides
└── api/            # API reference

.github/workflows/
├── ci.yml          # Run tests on push/PR
└── docs.yml        # Deploy documentation
```
