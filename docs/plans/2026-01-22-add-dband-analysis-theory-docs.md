# Plan: Add D-Band Center Analysis and Document Theoretical Framework

**Date:** 2026-01-22
**Status:** Completed

## Summary

Added d-band center analysis capability from VASP DOSCAR files and created comprehensive documentation of the theoretical framework implemented in NH3SOFC.

## Motivation

1. **D-band center** is a key descriptor for predicting adsorption energies on transition metals (Hammer-Norskov model)
2. Package had strong theory implementation but it wasn't well documented
3. Users need to understand the theoretical capabilities to use them effectively

## Solution

### Part 1: D-Band Analysis Module

Created `nh3sofc/analysis/electronic.py` with:

**DBandAnalyzer class:**
- Parse VASP DOSCAR files
- Calculate d-band center (1st moment): ε_d = ∫ε·ρ_d(ε)dε / ∫ρ_d(ε)dε
- Calculate d-band width (2nd moment)
- Calculate d-band filling
- Surface-averaged descriptors

**Convenience functions:**
- `calculate_d_band_center()`
- `get_surface_d_band_center()`

### Part 2: Theory Documentation

Created `docs/tutorials/theory.md` covering:

| Topic | Description |
|-------|-------------|
| Electronic Descriptors | D-band center, width, filling |
| Adsorption Energetics | E_ads calculations, coverage dependence |
| Thermochemistry | Harmonic (adsorbates), Ideal gas |
| BEP Relations | Built-in parameters, barrier estimation |
| Energy Span Model | Kozuch-Shaik TDTS/TDI identification |
| Microkinetics | Steady-state solver, TOF calculation |
| Volcano Plots | Descriptor-activity relationships |
| Exsolution | Vacancy formation, segregation energetics |

### Part 3: Website Updates

- Added "Theoretical Framework" card to main index
- Added feature table summarizing theory capabilities
- Updated tutorials index with "Analysis & Theory" section
- Added DBandAnalyzer to API documentation

## Files Modified

- `nh3sofc/analysis/electronic.py` - NEW: D-band analysis
- `nh3sofc/analysis/__init__.py` - Export new module
- `docs/tutorials/theory.md` - NEW: Theory documentation
- `docs/tutorials/index.md` - Added theory section
- `docs/index.md` - Added theory highlights
- `docs/api/analysis.md` - Added electronic structure docs

## Verification

```python
from nh3sofc.analysis import DBandAnalyzer, calculate_d_band_center
# Imports successful, module functional
```

## Commits

- `ca74141` - Add d-band center analysis and document theoretical framework

## References

1. Hammer, B. & Norskov, J.K. Adv. Catal. 45, 71-129 (2000) - D-band model
2. Kozuch, S. & Shaik, S. Acc. Chem. Res. 44, 101 (2011) - Energy span model
3. Norskov, J.K. et al. Nat. Chem. 1, 37 (2009) - Computational catalysis
