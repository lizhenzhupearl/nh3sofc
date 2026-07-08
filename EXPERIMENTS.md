# EXPERIMENTS.md — Live Research Log

> Rule 3: Write hypothesis and expected outcome BEFORE executing.
> Rule 9: End every entry with KEEP or DISCARD + one-sentence reason.
> Rule 12: DISCARD entries must say *why* and *what this rules out*.
> Rule 2: One variable changed per experiment.

---

<!-- TEMPLATE — copy for each new experiment

## EXP-XXX | 2026-03-25 | [One-line description]

**Variable changed:** [What exactly is different from the last experiment?]

**Hypothesis:** [What do we expect to happen and why?]

**Method:**
- [Step 1]
- [Step 2]

**Budget:** [e.g., 2h on Hyades GPU 4 / 4h on Aspire2a A100 / 30min local]

**Oracle:** [Metric + threshold — e.g., val_loss < 0.3 after 100 epochs]

**Expected result:** [Quantitative prediction]

**Actual result:** [Fill in after running]

**Decision:** KEEP — [why this advances the work]
          or DISCARD — [what failed; what this rules out for future attempts]

---
-->

## EXP-001 | 2026-07-08 | MKM with literature DFT barriers (Ni(111), 0K energies)

**Variable changed:** Baseline — first MKM run on Ni(111)

**Hypothesis:** Literature DFT barriers (Duan et al. 2012) with MACE-computed
thermodynamics should reproduce Ea_app in the 133–227 kJ/mol range and give
the correct sign of reaction orders (n_NH3 > 0, n_H2 < 0).

**Method:**
- Barriers: NH3_ads=0, NH3→NH2=1.0, NH2→NH=0.9, NH→N=1.2, 2N→N₂=1.7, 2H→H₂=1.0 eV
- Reaction energies: from MACE NEB + profile data (0K, no thermal corrections)
- Gas pressures: P_NH3=0.1, P_N2=0.01, P_H2=0.01 bar
- Solve steady-state with robust solver, Arrhenius over 500–900K

**Budget:** 5 min local

**Oracle:** Ea_app ∈ [133, 227] kJ/mol; n_NH3 > 0; n_H2 < 0

**Expected result:** Ea_app ~180 kJ/mol, N* dominant (N₂ desorption is RDS on Ni)

**Actual result:**
- Ea_app = 337.5 kJ/mol — **FAIL** (too high by ~2x)
- n_NH3 = NaN (TOF too small for log differentiation)
- n_H2 = −1.8 — **PASS** (correct sign)
- TOF(673K) = 1.5e-11 s⁻¹ — **FAIL** (10+ orders of magnitude too low)
- Dominant: NH3* (64%), H* (36%) — **FAIL** (should be N* on Ni)
- Site balance: conserved (sum = 1.0)

**Root cause:** Using raw 0K DFT energies without thermodynamic corrections.
Key issue: dG_H2 = +1.061 eV (very endothermic) → H* blocks 36% of sites.
At 673K, TΔS(H₂ gas) ≈ 0.7 eV would reduce effective dG to ~+0.36 eV,
dramatically lowering H* coverage and increasing available sites.

**Decision:** DISCARD — raw DFT energies without ZPE + vibrational entropy +
gas-phase entropy corrections cannot reproduce experimental kinetics. The
MKM framework works mechanically (solver converges, site balance conserved)
but thermochemistry corrections are mandatory. Rules out: using 0K energies
directly in MKM at elevated temperatures.

---

## EXP-002 | 2026-07-08 | MKM with MACE NEB barriers (Ni(111), 0K energies)

**Variable changed:** Barrier source (MACE NEB instead of literature DFT)

**Hypothesis:** MACE barriers should give qualitatively different kinetics from
literature, given known MAE of 0.64 eV for MACE NEB barriers vs DFT.

**Method:**
- Same as EXP-001 but barriers from MACE NEB:
  NH3→NH2=0.804, NH2→NH=1.967, NH→N=1.406, 2N→N₂=0.606, 2H→H₂=1.0 eV
- Same MACE thermodynamics (0K)

**Budget:** 5 min local

**Oracle:** Same as EXP-001

**Expected result:** Worse than EXP-001 due to known MACE barrier errors

**Actual result:**
- Ea_app = −275.6 kJ/mol — **FAIL** (nonsensical, negative)
- TOF(673K) = −2.0e-7 s⁻¹ — **FAIL** (negative = running backwards)
- Dominant: NH3* (64%), H* (36%) — same as EXP-001
- MACE gives wrong RDS: NH₂→NH barrier (1.97 eV) is 2x literature (0.9 eV)

**Decision:** DISCARD — MACE NEB barriers are unreliable for Ni NH3
decomposition (MAE 0.64 eV, wrong RDS identification, produces negative TOF).
Rules out: using MACE barriers directly for production MKM on this system.
DFT NEB validation is required for any barrier-sensitive analysis.

---

## Scientific Gaps Found (T1.3)

1. **Thermochemistry corrections are missing** — the MKM needs Gibbs free
   energies at operating temperature, not 0K DFT energies. The existing
   `analysis/thermochemistry.py` module has ZPE/thermal correction tools
   but they are not integrated into the MKM workflow. This is the
   #1 blocker for meaningful MKM results.

2. **MACE barriers are unreliable** for this system (MAE 0.64 eV). DFT NEB
   is needed for production barriers. This means the Ni benchmark must use
   literature DFT barriers until we run our own DFT NEBs.

3. **H₂ desorption thermodynamics dominate the kinetics** — the model is
   extremely sensitive to the H binding energy and H₂ desorption free
   energy. This parameter must be validated carefully.

4. **`sensitivity_analysis()` has a pre-existing API bug** — the `parameter`
   argument is required but unused. Should be made optional.
