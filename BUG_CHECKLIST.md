# BUG_CHECKLIST.md — Scientific Code Review Protocol

> Use for every non-trivial piece of computational code.
> Minimum 3 rounds. Do not skip Round 3.

---

## Pre-Code (Plan Before You Write)

- [ ] **Structure sketch**: write down data flow — what goes in, what comes out at each step
- [ ] **Connection map**: how does this connect to existing modules? What interfaces must it respect?
- [ ] **Units declaration**: declare upfront what units every input/output is in (write as comment block at top of function)
- [ ] **Double-check the plan**: read once more before writing a single line

---

## Round 1 — Self-Review (syntax, logic, science)

**Units and conversions:**
- [ ] Units consistent throughout?
- [ ] Conversion factors correct? (eV↔Hartree↔kJ/mol; Å↔Bohr↔nm; K↔eV via k_B)
- [ ] Prefactors right? (k_B = 8.617e-5 eV/K; ħ = 6.582e-16 eV·s; e = 1.602e-19 C; N_A = 6.022e23)

**Signs and directions:**
- [ ] Sign correct? (forces = −∇E; energy differences; adsorption energies negative for binding)
- [ ] Reference frame correct? (relative vs absolute; slab energy references)

**Numerics:**
- [ ] Normalization correct? (per-atom? per-cell? per-sample? per-batch?)
- [ ] Indexing/axis/dimension correct? (wrong axis, off-by-one, transposed shape)
- [ ] Numerical precision adequate? (float32 vs float64; catastrophic cancellation risk)
- [ ] Exponents right? (e.g., 1e-3 vs 1e3 is a factor of 1e6)

**Domain-specific:**
- [ ] Periodic boundary conditions / cell wrapping handled correctly?
- [ ] Symmetry violations? (result should respect symmetry of the system)
- [ ] Convergence thresholds adequate? (EDIFF, EDIFFG in VASP; loss plateaus in ML)
- [ ] Statistical bugs? (wrong denominator, sampling bias, incorrect aggregation)

---

## Round 2 — AI-Assisted Review

**Prompt to use:**
> "Review this code specifically for: unit errors, sign errors, conversion factor errors, normalization bugs, indexing bugs, and numerical precision issues. Show your reasoning for each."

After first response, always say:
> "Check again — find more."

Repeat until no new issues surface.

**Cross-validate (Rule 6):** For critical derivations or numerical outputs, have a second model check the first.

- [ ] AI review done
- [ ] "Check again" cycle done (at least once)
- [ ] Cross-validated with second model (if critical)

---

## Round 3 — Oracle / Sanity Check

- [ ] **Limiting case test**: does code reproduce known analytical answer in simple solvable case?
  - Test case: [describe]
  - Known answer: [value]
  - Result: [passed/failed]

- [ ] Physical reasonableness: energy in expected range? geometry stable? loss decreasing?
- [ ] Reference implementation match: [e.g., CLASS for Boltzmann; known DFT benchmark]
- [ ] For iterative methods: did it actually *converge*, or just stop?
- [ ] For periodic systems: PBC/wrapping effects handled correctly?
- [ ] Normalization explicitly documented in code comments?
- [ ] Sign convention documented in code comments?

---

## Science-Specific Bug Categories (always check these)

| Category | What to look for |
|----------|-----------------|
| PBC/wrapping | Cell translation, minimum image convention, fractional coordinates |
| Conversion factors | Every unit transition in the code — no mental math |
| Reference frames | Slab energy reference, relative vs absolute energies |
| Symmetry | Result should respect physical symmetry of the system |
| Statistical | N vs N-1 denominators, sampling bias, incorrect aggregation |
| Convergence | Threshold vs tolerance; "converged" vs "stopped" |

---

*File: [filename] | Date: 2026-03-25 | Reviewer: [name/model]*
