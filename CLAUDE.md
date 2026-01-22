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

- Verify imports work after adding new modules
- Run functional tests for new features
- Check backward compatibility when modifying APIs

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

docs/
├── plans/          # Implementation plans (SAVE PLANS HERE!)
├── tutorials/      # Step-by-step guides
└── api/            # API reference
```
