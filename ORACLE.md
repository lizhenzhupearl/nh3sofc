# ORACLE.md — Objective Success Metrics

> Rule 1: Oracle First. These metrics are defined BEFORE experiments run.
> If you can't define the oracle before running, you're not ready to run.

---

## Primary Oracle — Ni Benchmark

NH3 decomposition microkinetics on Ni, validated against experimental data.

- **Metric**: Apparent activation energy (Ea_app) and reaction orders (n_NH3, n_H2)
  from Arrhenius analysis of MKM-predicted TOFs over 623–873 K.
- **Baseline**: Literature experimental values:
  - Ea_app = 163–197 kJ/mol (Choudhary et al. 2001: 163 kJ/mol on Ni/SiO2;
    Ganley et al. 2004: 197 kJ/mol on Ni wire)
  - n_NH3 ≈ 0.5–1.0 (positive, sub-first-order at high conversion)
  - n_H2 ≈ −0.5 to −1.5 (negative, H2 inhibits by site blocking)
- **Target**: MKM with DFT barriers reproduces Ea_app within ±30 kJ/mol of
  the literature range (i.e., 133–227 kJ/mol) and correct signs of reaction
  orders (n_NH3 > 0, n_H2 < 0).
- **How to measure**:
  ```bash
  # Once T1.2 populates the energetics:
  python -c "
  from nh3sofc.analysis.microkinetic import arrhenius_analysis
  result = arrhenius_analysis(barriers, reaction_energies,
                              temperatures=[623,673,723,773,823,873])
  print(f'Ea_app = {result[\"apparent_Ea\"] * 96.485:.1f} kJ/mol')
  "
  ```

---

## Secondary Metrics

| Metric | Baseline | Target | Why it matters |
|--------|----------|--------|----------------|
| TOF at 673 K | ~0.1–10 s⁻¹ (typical Ni) | Within 2 orders of magnitude of experiment | Validates absolute rate scale |
| N* coverage at 673 K | Dominant on Ni (literature) | θ_N > 0.3 | N₂ desorption is RDS on Ni — N must accumulate |
| Campbell degree of rate control (N₂ step) | Close to 1.0 (RDS) | X_RC(N₂) > 0.5 | Confirms N₂ desorption is rate-determining |

---

## Regression Guard

*Must NOT get worse than these thresholds (regression = automatic DISCARD):*

| Metric | Floor | Consequence of breach |
|--------|-------|----------------------|
| MKM solver convergence | fsolve converges (info=1) for all T in [400,900] K | Discard run, debug solver |
| Site balance | |Σθ − 1| < 1e-8 | Discard, conservation violated |
| All coverages ≥ 0 | θ_i ≥ −1e-10 for all i | Discard, unphysical solution |
| Analytic two-step test | Relative error < 1e-8 | Block merge (pytest regression) |

---

## Sanity Checks

*Physical / mathematical limits the oracle must respect:*

- [ ] Rate constants k > 0 for all reactions at all temperatures
- [ ] All coverages in [0, 1] and sum to total_sites
- [ ] TOF > 0 (net forward direction for decomposition)
- [ ] Ea_app > 0 (reaction is thermally activated)
- [ ] n_H2 < 0 on Ni (H₂ inhibits — well-established experimentally)
- [ ] ΔG closed cycle = 0 (thermodynamic consistency: sum of reaction energies
      around any closed cycle must be zero)
- [ ] Higher T → higher TOF (Arrhenius behavior in the kinetic regime)

---

## Limiting Case Test

*The most powerful sanity check: does the solver reproduce a known analytical answer?*

- **Test case**: Two-step unimolecular surface reaction
  - Step 1: A(g) + * → A*,  rate r₁ = k₁ · P_A · θ_*
  - Step 2: A* → B(g) + *,  rate r₂ = k₂ · θ_A
  - Site balance: θ_A + θ_* = 1

- **Known answer** (steady-state):
  - θ_A = k₁·P_A / (k₁·P_A + k₂)
  - θ_* = k₂ / (k₁·P_A + k₂)
  - TOF = k₁·k₂·P_A / (k₁·P_A + k₂)

- **How to run**:
  ```bash
  pytest tests/test_mkm_oracle.py -v
  ```

- **Acceptable tolerance**: relative error < 1e-8

---

## Key References

| Ref | Value | Source |
|-----|-------|--------|
| Ea_app Ni wire | 197 kJ/mol | Ganley et al., J. Catal. 227 (2004) 26 |
| Ea_app Ni/SiO₂ | 163 kJ/mol | Choudhary et al., Catal. Lett. 72 (2001) 197 |
| n_NH3 (Ni) | ~0.75 | Yin et al., Appl. Catal. A 277 (2004) 1 |
| n_H2 (Ni) | ~−1.0 | Yin et al., Appl. Catal. A 277 (2004) 1 |
| DFT barriers (Ni) | N₂ RDS, Ea~1.3–1.9 eV | Hansgen et al., J. Chem. Phys. 134 (2010) 184701 |

---

*Last updated: 2026-07-08*
