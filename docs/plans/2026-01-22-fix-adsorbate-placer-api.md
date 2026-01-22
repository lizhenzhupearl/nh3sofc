# Plan: Fix AdsorbatePlacer API Documentation Bugs and Enhance Implementation

**Date:** 2026-01-22
**Status:** Completed

## Summary

Fixed documentation bugs where API examples didn't match actual implementation, then enhanced the implementation to support the documented (and more intuitive) API.

## Problems Found

### Bug 1: Wrong attribute access on Atoms objects
- `oxynitride.atoms.get_chemical_formula()` should be `oxynitride.get_chemical_formula()`
- `AdsorbatePlacer(oxynitride.atoms)` should be `AdsorbatePlacer(oxynitride)`

### Bug 2: Non-existent `site_type` parameter in `add_on_site()`
Documentation used `site_type="ontop"` but parameter didn't exist.

### Bug 3: Non-existent `n_configs` parameter in `add_with_collision()`
Documentation used `n_configs=10` but method only returned single config.

### Bug 4: Parameter name mismatches
- `orientation` vs `rotation`
- `seed` vs `random_seed`
- `nx, ny` vs `grid_size`

## Solution

### Phase 1: Fix Documentation
Updated `docs/quickstart.md` and `docs/tutorials/adsorbate_placement.md` to use correct parameter names.

### Phase 2: Enhance Implementation
Modified `nh3sofc/structure/adsorbates.py`:

1. **Added `n_configs` to `add_with_collision()`**
   - If `n_configs=1`: returns single `Atoms` or `None` (backward compatible)
   - If `n_configs>1`: returns `List[Atoms]`

2. **Added `site_type` to `add_on_site()`**
   - `"ontop"`: above surface atoms (existing behavior)
   - `"bridge"`: midpoints between atom pairs
   - `"hollow"`: centroids via Delaunay triangulation
   - Falls back to computed sites when CatKit unavailable

3. **Added backward compatibility aliases**
   - `add_simple()`: `orientation` → `rotation`
   - `add_random()`: `seed` → `random_seed`
   - `add_grid()`: `nx`, `ny` → `grid_size`

## Files Modified

- `docs/quickstart.md` - Fixed API examples
- `docs/tutorials/adsorbate_placement.md` - Fixed parameter names
- `nh3sofc/structure/adsorbates.py` - Added features and aliases

## Verification

All tests passed:
```python
# Parameter signatures verified
# Functional tests: add_simple, add_random, add_grid, add_with_collision, add_on_site
# Bridge and hollow site fallbacks work without CatKit
```

## Commits

- `54499a9` - Fix AdsorbatePlacer API docs and enhance implementation
