# Tutorials

Step-by-step guides for common workflows in NH3-SOFC research.

## Getting Started

| Tutorial | Description | Time |
|----------|-------------|------|
| [Surface Building](surface_building.md) | Create surfaces and oxynitride defects | 15 min |
| [Adsorbate Placement](adsorbate_placement.md) | 6 methods for placing molecules | 20 min |
| [VASP Calculations](vasp_calculations.md) | Generate inputs, parse outputs | 25 min |

## NH3 Decomposition

| Tutorial | Description | Time |
|----------|-------------|------|
| [Decomposition Pathway](decomposition.md) | Complete NH3→N2+H2 workflow | 30 min |
| [NEB Transition States](neb.md) | Find reaction barriers | 25 min |
| [Frequency & Thermochemistry](thermochemistry.md) | ZPE, entropy, Gibbs energy | 20 min |

## Analysis & Theory

| Tutorial | Description | Time |
|----------|-------------|------|
| [Theoretical Framework](theory.md) | D-band analysis, BEP relations, energy span model | 30 min |
| [Microkinetic Modeling](microkinetics.md) | TOF and steady-state analysis | 35 min |
| [Surface Comparison](surface_comparison.md) | Ranking catalysts with volcano plots | 20 min |

## Advanced Topics

| Tutorial | Description | Time |
|----------|-------------|------|
| [Doped Ceria (GDC/SDC/PDC)](doped_ceria.md) | Acceptor-doped ceria with O vacancies | 25 min |
| [Exsolution Simulation](exsolution.md) | Metal nanoparticle exsolution from perovskites | 35 min |
| [High-Throughput Screening](screening.md) | Systematic parameter scans | 30 min |
| [MACE ML Force Fields](mace.md) | Training and active learning | 40 min |

## Workflow Diagrams

### Complete NH3 Decomposition Study

```
┌─────────────────┐
│  Build Surface  │
│  from Bulk CIF  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Create Defects  │
│ (O→N, vacancies)│
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Place NH3      │
│  (6 methods)    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Relax Structure │
│  (VASP/MACE)    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Decomposition   │
│  Intermediates  │
│  NH2*, NH*, N*  │
└────────┬────────┘
         │
    ┌────┴────┐
    │         │
    ▼         ▼
┌───────┐ ┌────────┐
│  NEB  │ │  Freq  │
│  TS   │ │  ZPE   │
└───┬───┘ └───┬────┘
    │         │
    └────┬────┘
         │
         ▼
┌─────────────────┐
│    Analysis     │
│ • Energy profile│
│ • RDS/barriers  │
│ • Gibbs ΔG      │
│ • Microkinetics │
└─────────────────┘
```

### Calculation Types

| Calc Type | VASP INCAR | Purpose |
|-----------|------------|---------|
| `relax` | IBRION=2, NSW=300 | Geometry optimization |
| `static` | NSW=0 | Single-point energy |
| `frequency` | IBRION=5, NFREE=2 | Vibrational modes |
| `neb` | IMAGES=5, LCLIMB=T | Transition states |
| `md` | IBRION=0, SMASS=0 | Molecular dynamics |

## Best Practices

### 1. Start Simple
- Begin with a small test system (2x2 surface)
- Use lower ENCUT (400 eV) for initial tests
- Verify convergence before production runs

### 2. Systematic Studies
- Use consistent parameters across comparisons
- Save all structures to database
- Document calculation settings

### 3. Validation
- Compare with literature values for benchmark systems
- Check for imaginary frequencies after relaxation
- Verify energy convergence with k-points and ENCUT

## Example Scripts

All tutorials include complete, runnable code. Download example scripts:

```bash
# Clone with examples
git clone https://github.com/lizhenzhupearl/nh3sofc.git
cd NH3SOFC/examples
```

## Video Tutorials

*(Coming soon)*
