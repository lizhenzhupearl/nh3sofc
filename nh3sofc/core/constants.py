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

# Exsolution-related constants
# Bulk metal energies per atom (VASP PBE) and atomic radii (A)
EXSOLUTION_METALS = {
    "Ni": {"bulk_energy": -5.51, "radius": 1.24, "atomic_number": 28},
    "Co": {"bulk_energy": -7.11, "radius": 1.25, "atomic_number": 27},
    "Fe": {"bulk_energy": -8.31, "radius": 1.26, "atomic_number": 26},
}

# Perovskite ABO3 site classification
PEROVSKITE_SITES = {
    "A_site": ["La", "Sr", "Ba", "Ca", "K", "Na", "Pr", "Nd", "Sm", "Gd"],
    "B_site": ["Ti", "V", "Mn", "Fe", "Co", "Ni", "Cu", "Zr", "Nb", "Cr", "Mo", "W"],
}

# Formal oxidation states for polarity calculations
# Used for estimating layer charges and dipole moments
DEFAULT_FORMAL_CHARGES = {
    # A-site cations (perovskites, typically +2 or +3)
    "La": +3, "Sr": +2, "Ba": +2, "Ca": +2, "K": +1, "Na": +1,
    "Pr": +3, "Nd": +3, "Sm": +3, "Gd": +3, "Y": +3,
    # B-site cations (perovskites, variable oxidation state)
    "Ti": +4, "V": +5, "Mn": +4, "Fe": +3, "Co": +3, "Ni": +2,
    "Cu": +2, "Zr": +4, "Nb": +5, "Cr": +3, "Mo": +6, "W": +6,
    # Reducible cations
    "Ce": +4,  # Can be +3 near vacancies
    # Anions
    "O": -2, "N": -3, "F": -1, "S": -2, "Cl": -1,
    # Adsorbate atoms (approximate)
    "H": +1, "C": 0,
}

# Structure type definitions for automatic detection
STRUCTURE_TYPES = {
    "perovskite": {
        "formula_pattern": "ABO3",  # General pattern
        "coordination": {"A": 12, "B": 6, "O": 2},  # Typical coordination numbers
        "examples": ["LaVO3", "SrTiO3", "BaTiO3", "LaMnO3"],
    },
    "rocksalt": {
        "formula_pattern": "MX",
        "coordination": {"M": 6, "X": 6},
        "examples": ["NiO", "MgO", "CoO", "FeO"],
    },
    "fluorite": {
        "formula_pattern": "MX2",
        "coordination": {"M": 8, "X": 4},
        "examples": ["CeO2", "ZrO2", "UO2"],
    },
}

# Surface polarity information for common facets
SURFACE_POLARITY = {
    "perovskite": {
        (0, 0, 1): {"polar": True, "layers": ["AO", "BO2"], "charges": [+1, -1]},
        (1, 1, 0): {"polar": True, "layers": ["ABO", "O2"], "charges": [+3, -4]},
        (1, 1, 1): {"polar": True, "layers": ["AO3", "B"], "charges": [-3, +5]},
    },
    "rocksalt": {
        (0, 0, 1): {"polar": False, "layers": ["MO"], "charges": [0]},
        (1, 1, 0): {"polar": False, "layers": ["MO"], "charges": [0]},
        (1, 1, 1): {"polar": True, "layers": ["M", "O"], "charges": [+2, -2]},
    },
    "fluorite": {
        (1, 1, 1): {"polar": False, "layers": ["O-M-O"], "charges": [0]},
        (1, 1, 0): {"polar": True, "layers": ["MO", "O"], "charges": [+2, -2]},
        (0, 0, 1): {"polar": True, "layers": ["M", "O2"], "charges": [+4, -4]},
    },
}

# Magic cluster sizes for stable nanoparticles
MAGIC_CLUSTER_SIZES = {
    "hemispherical": [1, 4, 7, 10, 13, 19],
    "icosahedral": [1, 13, 55, 147],
    "cuboctahedral": [1, 13, 55, 147],
}

# Default exsolution parameters
DEFAULT_EXSOLUTION_PARAMS = {
    "interface_distance": 2.0,  # A, metal-oxide interface distance
    "vacancy_coupling": 2,      # O vacancies per reduced B-site
    "socket_depth": 1.0,        # A, particle embedding depth
}

# Acceptor dopant properties for fluorite oxides (GDC, SDC, YSZ, etc.)
# Ionic radii in Angstroms for 8-coordinate sites (Shannon radii)
ACCEPTOR_DOPANTS = {
    # Trivalent dopants (2 dopants : 1 vacancy for charge compensation)
    "Sm": {"charge": +3, "ionic_radius": 1.079, "name": "Samarium"},
    "Gd": {"charge": +3, "ionic_radius": 1.053, "name": "Gadolinium"},
    "Pr": {"charge": +3, "ionic_radius": 1.126, "name": "Praseodymium"},
    "Y":  {"charge": +3, "ionic_radius": 1.019, "name": "Yttrium"},
    "La": {"charge": +3, "ionic_radius": 1.160, "name": "Lanthanum"},
    "Nd": {"charge": +3, "ionic_radius": 1.109, "name": "Neodymium"},
    "Tb": {"charge": +3, "ionic_radius": 1.040, "name": "Terbium"},  # Can also be +4 (0.88 Å)
    "Sc": {"charge": +3, "ionic_radius": 0.870, "name": "Scandium"},
    # Divalent dopants (1 dopant : 1 vacancy for charge compensation)
    "Ca": {"charge": +2, "ionic_radius": 1.120, "name": "Calcium"},
    "Mg": {"charge": +2, "ionic_radius": 0.890, "name": "Magnesium"},
}

# Host cation properties (8-coordinate ionic radii)
HOST_CATIONS = {
    "Ce": {"charge": +4, "ionic_radius": 0.970, "name": "Cerium", "oxide": "CeO2"},
    "Zr": {"charge": +4, "ionic_radius": 0.840, "name": "Zirconium", "oxide": "ZrO2"},
    "Hf": {"charge": +4, "ionic_radius": 0.830, "name": "Hafnium", "oxide": "HfO2"},
    "Th": {"charge": +4, "ionic_radius": 1.050, "name": "Thorium", "oxide": "ThO2"},
}

# Conventional material names for doped fluorites
# Maps (host_cation, dopant) tuples to common abbreviations
FLUORITE_MATERIAL_NAMES = {
    # Ceria-based
    ("Ce", "Gd"): "GDC",   # Gadolinium-doped Ceria
    ("Ce", "Sm"): "SDC",   # Samarium-doped Ceria
    ("Ce", "Pr"): "PDC",   # Praseodymium-doped Ceria
    ("Ce", "Y"):  "YDC",   # Yttrium-doped Ceria
    ("Ce", "La"): "LDC",   # Lanthanum-doped Ceria
    ("Ce", "Nd"): "NDC",   # Neodymium-doped Ceria
    ("Ce", "Tb"): "TDC",   # Terbium-doped Ceria
    # Zirconia-based
    ("Zr", "Y"):  "YSZ",   # Yttria-Stabilized Zirconia
    ("Zr", "Sc"): "ScSZ",  # Scandia-Stabilized Zirconia
    ("Zr", "Ca"): "CSZ",   # Calcia-Stabilized Zirconia
    ("Zr", "Mg"): "MSZ",   # Magnesia-Stabilized Zirconia
    # Hafnia-based
    ("Hf", "Y"):  "YSH",   # Yttria-Stabilized Hafnia
    # Thoria-based
    ("Th", "Ca"): "CTh",   # Calcia-doped Thoria
    ("Th", "Y"):  "YTh",   # Yttria-doped Thoria
}

# Charge compensation ratios for acceptor doping
# Maps dopant charge to (n_dopants, n_vacancies) for charge neutrality
# For trivalent dopants in tetravalent host: 2 M³⁺ + V_O^{2+} → net neutral
CHARGE_COMPENSATION = {
    3: (2, 1),  # 2 M³⁺ + 1 V_O^{2+} → net neutral (e.g., Gd³⁺ in Ce⁴⁺)
    2: (1, 1),  # 1 M²⁺ + 1 V_O^{2+} → net neutral (for divalent dopants)
}
