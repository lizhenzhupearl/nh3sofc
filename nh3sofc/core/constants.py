"""Physical constants and default parameters for NH3-SOFC calculations."""

# Physical constants
KB_EV = 8.617333262e-5  # Boltzmann constant in eV/K
KB_J = 1.380649e-23     # Boltzmann constant in J/K
H_EV_S = 4.135667696e-15  # Planck constant in eV·s
H_J_S = 6.62607015e-34    # Planck constant in J·s
AVOGADRO = 6.02214076e23  # Avogadro's number
EV_TO_KJ_MOL = 96.485     # eV to kJ/mol conversion
EV_TO_J = 1.602176634e-19  # eV to J conversion
C_CMS = 2.99792458e10     # Speed of light in cm/s

# Default VASP parameters
DEFAULT_VASP_PARAMS = {
    "encut": 520,           # Plane-wave cutoff (eV)
    "ediff": 1e-6,          # Electronic convergence (eV)
    "ediffg": -0.02,        # Force convergence (eV/A)
    "ismear": 0,            # Gaussian smearing for insulators
    "sigma": 0.05,          # Smearing width (eV)
    "kspacing": 0.03,       # K-point spacing (1/A)
    "ispin": 2,             # Spin-polarized
    "lorbit": 11,           # Write DOSCAR and lm-decomposed PROCAR
    "lwave": False,         # Don't write WAVECAR
    "lcharg": False,        # Don't write CHGCAR
    "ncore": 4,             # Parallelization
}

# Hubbard U values for transition metals (eV)
# From Materials Project / literature
HUBBARD_U = {
    "V": 3.25,   # Vanadium
    "Ti": 3.00,  # Titanium
    "Mn": 3.90,  # Manganese
    "Fe": 5.30,  # Iron
    "Co": 3.32,  # Cobalt
    "Ni": 6.20,  # Nickel
    "Cu": 4.00,  # Copper
    "Zr": 0.00,  # Zirconium (usually no U)
    "Ce": 5.00,  # Cerium
    "La": 0.00,  # Lanthanum (usually no U)
}

# Default surface parameters
DEFAULT_SURFACE_PARAMS = {
    "vacuum": 15.0,         # Vacuum thickness (A)
    "layers": 6,            # Number of layers
    "fix_bottom": 2,        # Fixed bottom layers
}

# Van der Waals correction methods
VDW_METHODS = {
    "D3": {"ivdw": 11},           # DFT-D3 (zero damping)
    "D3BJ": {"ivdw": 12},         # DFT-D3 with BJ damping
    "D4": {"ivdw": 13},           # DFT-D4 (if available)
    "TS": {"ivdw": 20},           # Tkatchenko-Scheffler
    "TSHIRSH": {"ivdw": 21},      # TS with Hirshfeld partitioning
}

# Common adsorbate molecules and their properties
ADSORBATES = {
    "NH3": {
        "formula": "NH3",
        "atoms": ["N", "H", "H", "H"],
        "n_atoms": 4,
        "anchor_atom": 0,  # N is anchor
        "height": 2.0,     # Default height above surface (A)
        "symmetry": 3,     # C3v symmetry
    },
    "NH2": {
        "formula": "NH2",
        "atoms": ["N", "H", "H"],
        "n_atoms": 3,
        "anchor_atom": 0,
        "height": 1.8,
        "symmetry": 2,
    },
    "NH": {
        "formula": "NH",
        "atoms": ["N", "H"],
        "n_atoms": 2,
        "anchor_atom": 0,
        "height": 1.5,
        "symmetry": 1,
    },
    "N": {
        "formula": "N",
        "atoms": ["N"],
        "n_atoms": 1,
        "anchor_atom": 0,
        "height": 1.2,
        "symmetry": 1,
    },
    "H": {
        "formula": "H",
        "atoms": ["H"],
        "n_atoms": 1,
        "anchor_atom": 0,
        "height": 1.0,
        "symmetry": 1,
    },
    "H2": {
        "formula": "H2",
        "atoms": ["H", "H"],
        "n_atoms": 2,
        "anchor_atom": 0,
        "height": 2.5,
        "symmetry": 2,
    },
    "N2": {
        "formula": "N2",
        "atoms": ["N", "N"],
        "n_atoms": 2,
        "anchor_atom": 0,
        "height": 2.5,
        "symmetry": 2,
    },
    "H2O": {
        "formula": "H2O",
        "atoms": ["O", "H", "H"],
        "n_atoms": 3,
        "anchor_atom": 0,
        "height": 2.0,
        "symmetry": 2,
    },
    "O": {
        "formula": "O",
        "atoms": ["O"],
        "n_atoms": 1,
        "anchor_atom": 0,
        "height": 1.2,
        "symmetry": 1,
    },
    "OH": {
        "formula": "OH",
        "atoms": ["O", "H"],
        "n_atoms": 2,
        "anchor_atom": 0,
        "height": 1.5,
        "symmetry": 1,
    },
}

# Reference gas phase energies (placeholder - should be calculated)
# These are approximate DFT-PBE values and should be updated
GAS_PHASE_ENERGIES = {
    "H2": -6.77,    # eV (VASP PBE)
    "N2": -16.64,   # eV
    "NH3": -19.54,  # eV
    "H2O": -14.22,  # eV
    "O2": -9.86,    # eV (triplet)
}

# Typical operating conditions for NH3-SOFC
OPERATING_CONDITIONS = {
    "temperature_low": 400,   # K (~127 C)
    "temperature_mid": 673,   # K (~400 C)
    "temperature_high": 873,  # K (~600 C)
    "pressure_atm": 1.0,      # atm
}

# Tolerance values
TOLERANCE = {
    "distance": 1e-6,         # A
    "energy": 1e-8,           # eV
    "angle": 1e-4,            # degrees
    "rmsd": 0.5,              # A (for structure comparison)
}
