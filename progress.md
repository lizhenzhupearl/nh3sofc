# Session Progress Log

## Current Verified State

**Last Updated:** 2026-07-08
**Repo Root:** `/Users/zhenzhu/NH3SOFC`
**Startup Path:** `./init.sh`
**Verification Path:** `pytest tests/ -v` (114 tests passing)
**All baseline features:** passing (feat-001 through feat-010, T0.1–T0.4)
**Current Blockers:** None

## Session Log

### Session 1 — 2026-06-29

**Goal:** Set up harness engineering system for the project

**Completed:**
- Merged harness patterns (startup workflow, working rules, definition of done, session lifecycle) into CLAUDE.md
- Customized init.sh for Python/NH3SOFC
- Populated feature-list.json with 5 real features
- Wired research discipline files into CLAUDE.md
- Integrated Karpathy's 10 rules as condensed Code Discipline section
- Moved plan.md to docs/plans/

**Verification:** init.sh passes, 71 tests passing

**Commits:** 750f808

### Session 2 — 2026-06-29

**Goal:** Review uncommitted dopants.py changes, write tests, add backward compat

**Completed:**
- Reviewed 16-line diff generalizing analyze_dopant_distribution() for non-Ce hosts
- Added 7 regression tests + 1 smoke test for host_cation support
- Added ce_total backward-compat alias for Ce hosts

**Verification:** 79 tests passing

### Session 3 — 2026-06-29/30

**Goal:** Harnessed codebase audit (feat-006 through feat-008)

**Completed:**
- **Structure audit (feat-006):** Fixed 2 bugs (cluster atom count, molecule fallback), 4 regression tests
- **Calculators audit (feat-007):** Fixed 5 bugs (INCAR case, Hubbard U f-orbitals, VDW case, cm-1 conversion, is_converged false positive), 12 regression tests
- **Analysis audit (feat-008):** Fixed 2 bugs (DOSCAR LORBIT=11, set_bulk_energy), verified all conversion factors correct, 6 regression tests
- Merged all fixes to main

**Verification:** 101 tests passing

**Commits:** 61880cf

**Known Risks:** None

**Next Best Step:** Complete remaining audits (feat-009: workflows, feat-010: database+core)

### Session 4 — 2026-06-29

**Goal:** Audit nh3sofc/database/ + nh3sofc/core/ (feat-010)

**Completed:**
- Audited 8 files: core/__init__.py, core/constants.py, core/base.py, core/io.py, database/__init__.py, database/naming.py, database/asedb.py, cli/__init__.py
- **BUG 1 (critical):** `get_surface_atoms()` ignored `z_threshold` parameter entirely; when `layer_tolerance=None`, it hardcoded 2.0 instead of using `z_threshold`. Fixed to fall back to fractional threshold.
- **BUG 2 (moderate):** `asedb.py` `add_structure()` stored Miller indices without `abs()`, inconsistent with `naming.py` which uses `abs()`. Fixed.
- **BUG 3 (moderate):** `get_decomposition_energies()` passed `miller="001"` but DB stores `"m001"`, so queries silently returned empty. Fixed to add "m" prefix.
- **LINT 1:** Redundant `isinstance` branch in `query()` kwargs loop (both branches did same thing). Fixed.
- **Noted but not fixed (out of scope):** `SURFACE_POLARITY` imported but unused in `surface.py`; `STRUCTURE_TYPES`, `TOLERANCE`, `KB_J`, `H_J_S` defined but never used in constants.py (intentional placeholders)
- Verified all physical constants correct (KB_EV, H_EV_S, EV_TO_KJ_MOL, etc.)
- 6 regression tests added, 77 tests passing

**Verification:** 77 tests passing

**Known Risks:** None

**Next Best Step:** Complete feat-009 (workflows audit)

### Session 5 — 2026-06-30

**Goal:** Complete final audit (feat-009: workflows) and merge all remaining fixes

**Completed:**
- **Workflows audit (feat-009):** Fixed 3 bugs (get_gibbs_energy TypeError, couple_with_decomposition wrong kwarg, NEB dead code), 2 lint fixes, 5 regression tests
- Merged feat-009 and feat-010 fixes to main
- Cleaned up all worktree branches
- Updated all harness artifacts

**Verification:** 110 tests passing

**Commits:** 71db775

**Known Risks:** None

### Session 6 — 2026-07-08

**Goal:** Execute Phase 0 of nh3sofc_code_review_and_plan_v3.md (make harness own the plan)

**Completed:**
- **T0.1:** Imported all 20 tasks (T0.1–T4.3) into feature-list.json with track labels (S/E), dependencies, spec references
- **T0.2:** Added breaking-change deprecation shim policy to CLAUDE.md Working Rules (closes H3)
- **T0.3:** Added metrics ratchet line to CLAUDE.md Definition of Done (closes H4 partially)
- **T0.4:** Created `scripts/metrics.py` (AST+tokenize-based, not regex), `tests/test_metrics_ratchet.py` (4 tests), `metrics_baseline.json`
- Copied v3 plan to `docs/plans/nh3sofc_code_review_and_plan_v3.md`
- Added `docs/claude/` reference files (key-directories, plan-management, testing-guide)

**Authoritative baseline (AST-based, corrects v2's regex estimates):**
- total_functions: 530, annotated: 448 (84.5%), print_calls: 226, broad_excepts: 21

**Verification:** 114 tests passing (110 existing + 4 ratchet)

**Known Risks:** None

## Notes for Next Session

- Phase 0 complete — harness now owns the v3 plan
- Ratchets are active: new print() calls or broad excepts will fail tests
- Next: Phase 1 (T1.1 oracle + T1.2 Ni benchmark) — science track starts
