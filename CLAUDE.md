# Claude Code Instructions for NH3SOFC

## Project Overview

NH3SOFC is a Python package for computational catalysis research focused on ammonia decomposition for solid oxide fuel cells. It provides structure generation, VASP/MACE calculations, and analysis tools.

## Startup Workflow

Before writing any code, every session must follow this sequence:

1. **Confirm working directory** with `pwd` (must be the NH3SOFC repo root)
2. **Read this file** completely
3. **Run `./init.sh`** to verify the environment is healthy
4. **Read `feature-list.json`** to see current feature state
5. **Read `progress.md`** for cross-session context
6. **Review recent commits** with `git log --oneline -5`

If baseline verification is failing, **repair that first** before adding any new scope.

### Research Session Startup (additional steps)

When the session involves computation, experiments, or scientific analysis, also do:

7. **Read `RESEARCH.md`** — current problem, hypothesis, approach
8. **Read last 3 entries in `EXPERIMENTS.md`** — where are we?
9. **Read `ORACLE.md`** — what metric and threshold define success?
10. **Declare the budget** — time + compute for this session

## Working Rules

- **One feature at a time**: Pick exactly one unfinished feature from `feature-list.json`
- **Verification required**: Don't claim done without running `./init.sh` or `pytest tests/`
- **Update artifacts**: Before ending session, update `progress.md` and `feature-list.json`
- **Stay in scope**: Don't modify files unrelated to the current feature
- **Leave clean state**: Next session must be able to run `./init.sh` immediately
- **Deprecation shims for breaking changes**: Renaming or removing any public API, function, or returned dict key requires a deprecation shim for one minor version (follow the existing `get_dopant_name` → `get_material_name` pattern). Ship the shim in the same commit as the rename; remove it one minor version later.

## Code Discipline (condensed from [karpathy-rules.md](karpathy-rules.md))

- **Think before you code**: State assumptions. Name tradeoffs. If confused, stop and ask.
- **Simplicity**: Minimum code that solves *this* problem. No premature abstractions.
- **Surgical changes**: Small diffs. Don't touch unrelated code. Match existing style.
- **Dependencies**: Prefer what's already in the project (ASE, NumPy, SciPy). Justify additions.
- **Communication**: Say what you did and why. Flag concerns proactively. Be precise about uncertainty.
- **Failure modes**: Watch for Kitchen Sink, Wrong Abstraction, Invisible Decision, Knowledge Hallucination, Runaway Refactor.

## Definition of Done

- [ ] Target behavior is implemented
- [ ] Verification actually ran (`pytest tests/` + import smoke check)
- [ ] Metrics ratchet tests pass via `./init.sh` (no regression in print count, broad excepts, or annotation fraction)
- [ ] Evidence recorded in `feature-list.json` and `progress.md`
- [ ] Repository remains restartable via `./init.sh`
- [ ] For non-trivial computation code: `BUG_CHECKLIST.md` 3-round review completed

## End of Session

1. Update `progress.md` with current state
2. Update `feature-list.json` with new feature status and evidence
3. Record any unresolved risks or blockers
4. Park tangents in `IDEAS.md` (Rule 13 — don't chase)
5. If experiments ran: update `EXPERIMENTS.md` with result and KEEP/DISCARD decision
6. If hypothesis changed: update `RESEARCH.md`
7. Commit with descriptive message once work is in safe state

## Escalation

- **Architecture decisions**: Check `docs/plans/` and project docs, otherwise ask user
- **Unclear requirements**: Check `feature-list.json` for definition of done, otherwise ask user
- **Repeated test failures**: Update `progress.md`, flag for human review
- **Scope ambiguity**: Re-read `feature-list.json` for the active feature's description

## Code Style

- Follow existing patterns in the codebase
- Use type hints
- Include docstrings with Examples section
- Prefer editing existing files over creating new ones

## Reference (read on demand)

| Topic | File |
|-------|------|
| Plan management & templates | [`docs/claude/plan-management.md`](docs/claude/plan-management.md) |
| Testing guide & regression template | [`docs/claude/testing-guide.md`](docs/claude/testing-guide.md) |
| Directory map & verification commands | [`docs/claude/key-directories.md`](docs/claude/key-directories.md) |

## Harness Artifacts

| File | Purpose |
|------|---------|
| `init.sh` | Reproducible boot sequence — install, verify, report |
| `feature-list.json` | Feature state tracker (source of truth for scope) |
| `progress.md` | Cross-session continuity log |

## Research Artifacts

| File | Purpose | When to read |
|------|---------|--------------|
| `RESEARCH.md` | Strategy — problem, hypothesis, approach | Every research session start |
| `ORACLE.md` | Success metrics — metric + threshold | Before running any experiment |
| `EXPERIMENTS.md` | Live experiment log — KEEP/DISCARD | Session start + after each experiment |
| `IDEAS.md` | Parked tangents (Rule 13) | End of session + retrospectives |
| `LITERATURE.md` | Reference values, benchmarks, constants | When validating results |
| `BUG_CHECKLIST.md` | 3-round scientific code review | Non-trivial computation code |
| `RETROSPECTIVE.md` | Phase-level strategy review | At milestones or when stuck |
| `docs/plans/` | Dated implementation plans | When starting planned work |
