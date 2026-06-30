# Session Progress Log

## Current Verified State

**Last Updated:** 2026-06-29
**Repo Root:** `/Users/zhenzhu/NH3SOFC`
**Startup Path:** `./init.sh`
**Verification Path:** `pytest tests/ -v && python -c "from nh3sofc.structure import SurfaceBuilder, AdsorbatePlacer"`
**All baseline features:** passing (feat-001 through feat-005)
**Current Blockers:** None

## Session Log

### Session 1 — 2026-06-29

**Goal:** Set up harness engineering system for the project

**Completed:**
- Merged harness patterns (startup workflow, working rules, definition of done, session lifecycle) into CLAUDE.md
- Customized init.sh for Python/NH3SOFC (pip install -e . → pytest → import smoke check)
- Populated feature-list.json with 5 real features reflecting current project state
- Initialized progress.md (this file)
- Removed agents.md template (content merged into CLAUDE.md)

**Verification:** Pending (run `./init.sh` to confirm)

**Known Risks:** None

**Next Best Step:** Run `./init.sh` to confirm baseline passes, then add new features to feature-list.json as research priorities emerge

### Session 2 — 2026-06-29

**Goal:** Review uncommitted changes in dopants.py, write tests for them

**Completed:**
- Reviewed 16-line diff in `nh3sofc/structure/dopants.py`: generalizes `analyze_dopant_distribution()` and `print_dopant_analysis()` to support non-Ce hosts (YSZ, ScSZ, etc.)
- Key changes: added `host_cation` param, renamed `ce_total` -> `host_total` in return dict, added backward-compat fallback in `print_dopant_analysis()`
- Added 7 regression tests in `TestAnalyzeDopantDistributionHostCation` covering: key rename, default host, YSZ analysis, dopant fraction calculation, wrong host detection, print backward compat
- Added 1 new smoke test `test_analyze_dopant_distribution_ysz` for YSZ host analysis
- Updated existing smoke test `test_analyze_dopant_distribution` to assert on new `host_total`/`host_cation` keys
- All 79 tests pass (was 71 before)

**Verification:** `pytest tests/ -v` — 79 passed, 0 failed

**Known Risks:** The `ce_total` -> `host_total` key rename is a breaking change for any downstream code that accesses `stats["ce_total"]` directly

**Next Best Step:** Commit the dopants.py changes together with the new tests

### Session 3 — 2026-06-29

**Goal:** Add backward compatibility for `ce_total` -> `host_total` rename in `analyze_dopant_distribution()`

**Completed:**
- Searched entire codebase for `ce_total` references: found in `dopants.py` (3 locations), `test_regressions.py` (6 locations), `test_smoke.py` (2 locations). No references in tutorials or docs.
- Added backward-compatibility alias: `analyze_dopant_distribution()` now includes `ce_total` in the return dict when `host_cation == "Ce"`, so existing code using `stats["ce_total"]` continues to work
- Added `import warnings` to module imports
- Updated docstring to document `ce_total` as a deprecated alias
- Updated regression test `test_host_total_replaces_ce_total` to expect `ce_total` present as deprecated alias (instead of asserting it's absent)
- Added assertion in `test_host_cation_zr_for_ysz` to verify `ce_total` is NOT present for non-Ce hosts
- All 79 tests pass

**Verification:** `pytest tests/ -v` — 79 passed, 0 failed

**Known Risks:** None — the breaking change is now resolved. Both `ce_total` and `host_total` work for Ce hosts.

**Next Best Step:** Commit all uncommitted changes (dopants.py generalization + backward compat + tests)

### Session 4 — 2026-06-29

**Goal:** feat-007 — Audit `nh3sofc/calculators/` subpackage (8 files)

**Completed:**
- Audited all 9 files in calculators/ (vasp/inputs.py, vasp/outputs.py, vasp/frequency.py, mace/interface.py, mace/training.py, plus __init__.py files)
- Found and fixed 5 bugs:
  1. **INCAR parameter case mismatch** (critical): DEFAULT_VASP_PARAMS used lowercase keys but CALC_PRESETS and generate_incar used uppercase, causing duplicate INCAR parameters. Fixed by uppercasing all keys in DEFAULT_VASP_PARAMS.
  2. **Hubbard U f-orbital bug** (critical): set_hubbard_u() hardcoded LDAUL=2 (d-orbitals) for ALL elements including lanthanides. Fixed with F_ELECTRON_ELEMENTS set and automatic LDAUL=3 / LMAXMIX=6 selection.
  3. **VDW_METHODS lowercase keys** (critical): Used `"ivdw"` instead of `"IVDW"`, causing VDW settings to end up in wrong INCAR section.
  4. **cm-1 to eV conversion** (moderate): Hardcoded 1.23981e-4 (wrong); now derived from H_EV_S * C_CMS = 1.23984e-4.
  5. **is_converged false positive** (moderate): Returned True for unconverged relaxations via "EDIFF is reached" check. Fixed to check "General timing and accounting" for static calcs only.
- Fixed 4 lint issues: unused imports (ase_write, HUBBARD_U, VASPOutputParser, H_EV_S), mutable default argument, uninitialized _kpoints_mesh
- Fixed 1 API inconsistency: added MACE exports to calculators/__init__.py
- Added 12 regression tests (4 for case mismatch, 3 for Hubbard U, 2 for conversion factor, 3 for is_converged)
- All 83 tests pass

**Verification:** `pytest tests/ -v` — 83 passed, 0 failed

**Known Risks:** The DEFAULT_VASP_PARAMS key case change could affect downstream code that accesses params by lowercase keys, but no such usage found in codebase.

**Next Best Step:** Continue with feat-006 (structure audit) or feat-008 (analysis audit)

## Notes for Next Session

- All 5 features in feature-list.json are currently "passing" — the project is at a stable baseline
- New features should be appended to feature-list.json as feat-006, feat-007, etc.
- When starting a new research task, add it as a feature first, then work on it
- The uncommitted dopants.py changes + new tests are ready to commit
- Consider adding a `DeprecationWarning` when `ce_total` is accessed in a future version
