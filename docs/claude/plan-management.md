# Plan Management

Every implementation plan MUST be saved to `docs/plans/` with proper naming.

## Filename Format

`YYYY-MM-DD-descriptive-title.md` (e.g., `2026-01-22-fix-adsorbate-placer-api.md`)

## When to Save

- After a plan is approved and before starting implementation
- Update the plan with "Completed" status and commit references after finishing

## Plan Template

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

## Update the Index

Add entry to `docs/plans/README.md` index table.

## Why This Matters

Plans represent significant intellectual work and decision-making. They serve as:
- Decision records
- Methodology documentation
- Reproducibility artifacts
- Project history
