# ORACLE.md — Objective Success Metrics

> Rule 1: Oracle First. These metrics are defined BEFORE experiments run.
> If you can't define the oracle before running, you're not ready to run.

---

## Primary Oracle

- **Metric**: [e.g., validation bits-per-byte, adsorption energy MAE, force MAE]
- **Baseline**: [current best / reference implementation / literature value]
- **Target**: [what counts as success — be specific]
- **How to measure**:
  ```bash
  # Exact command or script
  [command]
  ```

---

## Secondary Metrics

| Metric | Baseline | Target | Why it matters |
|--------|----------|--------|----------------|
| [metric] | [value] | [value] | [reason] |

---

## Regression Guard

*Must NOT get worse than these thresholds (regression = automatic DISCARD):*

| Metric | Floor | Consequence of breach |
|--------|-------|----------------------|
| [metric] | [value] | [e.g., discard run, investigate] |

---

## Sanity Checks

*Physical / mathematical limits the oracle must respect:*

- [ ] [e.g., Energy must be negative for stable adsorption]
- [ ] [e.g., Loss must decrease monotonically during training (or explain why not)]
- [ ] [e.g., Forces must be zero at equilibrium geometry]

---

## Limiting Case Test

*The most powerful sanity check: does the code reproduce a known analytical answer?*

- **Test case**: [describe the simple solvable case]
- **Known answer**: [exact value]
- **How to run**:
  ```bash
  [command]
  ```
- **Acceptable tolerance**: [e.g., < 1e-4 relative error]

---

*Last updated: 2026-03-25*
