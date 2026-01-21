"""Exsolution structure generation for perovskite materials.

Creates structures for modeling exsolution processes where transition metal
cations (Ni, Co, Fe) migrate from perovskite bulk to surface under reducing
conditions, forming metallic nanoparticles.

Key features:
1. Perovskite site identification (A-site, B-site, O-site)
2. Defective perovskite with cation and anion vacancies
3. Surface segregation modeling
4. Metallic nanoparticle placement with socketing
5. Adsorption site identification for catalysis studies
"""

from pathlib import Path
from typing import Optional, List, Union, Dict, Any, Tuple
import numpy as np
from ase import Atoms
from ase.io import write
from ase.build import bulk as ase_bulk

from ..core.base import BaseStructure, get_surface_atoms
from ..core.constants import (
    EXSOLUTION_METALS,
    PEROVSKITE_SITES,
    MAGIC_CLUSTER_SIZES,
    DEFAULT_EXSOLUTION_PARAMS,
)
from .surface import SlabStructure


class ExsolutionBuilder:
    """
    Builder class for creating exsolution structures.

    Generates structures representing different stages of exsolution:
    1. Pristine perovskite surface
    2. Defective perovskite (with cation/anion vacancies)
    3. Surface-segregated metal
    4. Exsolved nanoparticle on surface

    Examples
    --------
    >>> slab = SlabStructure.from_file("LaSrTiNiO3_001.cif")
    >>> builder = ExsolutionBuilder(slab)
    >>>
    >>> # Create defective perovskite with vacancies
    >>> defective = builder.create_defective_perovskite(
    ...     a_site_vacancy_fraction=0.1,
    ...     b_site_vacancy_fraction=0.05,
    ...     oxygen_vacancy_fraction=0.08,
    ... )
    >>>
    >>> # Create Ni nanoparticle on surface
    >>> with_particle = builder.create_nanoparticle(
    ...     metal="Ni",
    ...     n_atoms=13,
    ...     shape="hemispherical",
    ... )
    """

    def __init__(
        self,
        structure: Union[SlabStructure, Atoms],
        a_site_elements: Optional[List[str]] = None,
        b_site_elements: Optional[List[str]] = None,
    ):
        """
        Initialize ExsolutionBuilder.

        Parameters
        ----------
        structure : SlabStructure or Atoms
            Surface slab structure (perovskite)
        a_site_elements : list, optional
            Elements occupying A-site (default from PEROVSKITE_SITES)
        b_site_elements : list, optional
            Elements occupying B-site (default from PEROVSKITE_SITES)
        """
        if isinstance(structure, SlabStructure):
            self.atoms = structure.atoms.copy()
            self._structure_obj = structure
        elif isinstance(structure, BaseStructure):
            self.atoms = structure.atoms.copy()
            self._structure_obj = structure
        else:
            self.atoms = structure.copy()
            self._structure_obj = None

        self.a_site_elements = a_site_elements or PEROVSKITE_SITES["A_site"]
        self.b_site_elements = b_site_elements or PEROVSKITE_SITES["B_site"]

        # Cache site identification
        self._site_indices = None

    def identify_perovskite_sites(self) -> Dict[str, List[int]]:
        """
        Identify A-site, B-site, and O-site indices in the perovskite structure.

        Uses element-based identification (A-site: La/Sr/Ba/Ca, B-site: Ti/Ni/Co/Fe).

        Returns
        -------
        dict
            Dictionary with keys 'A_site', 'B_site', 'O_site' containing atom indices
        """
        if self._site_indices is not None:
            return self._site_indices

        symbols = self.atoms.get_chemical_symbols()

        a_indices = []
        b_indices = []
        o_indices = []

        for i, symbol in enumerate(symbols):
            if symbol in self.a_site_elements:
                a_indices.append(i)
            elif symbol in self.b_site_elements:
                b_indices.append(i)
            elif symbol == "O":
                o_indices.append(i)

        self._site_indices = {
            "A_site": a_indices,
            "B_site": b_indices,
            "O_site": o_indices,
        }

        return self._site_indices

    def get_element_indices(self, element: str) -> List[int]:
        """Get indices of atoms of a specific element."""
        symbols = self.atoms.get_chemical_symbols()
        return [i for i, s in enumerate(symbols) if s == element]

    def create_defective_perovskite(
        self,
        a_site_vacancy_fraction: float = 0.0,
        b_site_vacancy_fraction: float = 0.0,
        oxygen_vacancy_fraction: float = 0.0,
        b_site_element: Optional[str] = None,
        random_seed: Optional[int] = None,
    ) -> Atoms:
        """
        Create defective perovskite with cation and anion vacancies.

        Parameters
        ----------
        a_site_vacancy_fraction : float
            Fraction of A-site atoms to remove (0.0 to 1.0)
        b_site_vacancy_fraction : float
            Fraction of B-site atoms to remove (0.0 to 1.0)
        oxygen_vacancy_fraction : float
            Fraction of O atoms to remove (0.0 to 1.0)
        b_site_element : str, optional
            Specific B-site element to create vacancies in (e.g., "Ni")
        random_seed : int, optional
            Random seed for reproducibility

        Returns
        -------
        Atoms
            Defective perovskite structure
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        sites = self.identify_perovskite_sites()
        vacancy_indices = []

        # A-site vacancies
        if a_site_vacancy_fraction > 0 and sites["A_site"]:
            n_remove = max(1, int(len(sites["A_site"]) * a_site_vacancy_fraction))
            selected = np.random.choice(
                sites["A_site"], size=min(n_remove, len(sites["A_site"])), replace=False
            )
            vacancy_indices.extend(selected.tolist())

        # B-site vacancies
        if b_site_vacancy_fraction > 0:
            if b_site_element:
                b_candidates = [
                    i for i in sites["B_site"]
                    if self.atoms.get_chemical_symbols()[i] == b_site_element
                ]
            else:
                b_candidates = sites["B_site"]

            if b_candidates:
                n_remove = max(1, int(len(b_candidates) * b_site_vacancy_fraction))
                selected = np.random.choice(
                    b_candidates, size=min(n_remove, len(b_candidates)), replace=False
                )
                vacancy_indices.extend(selected.tolist())

        # Oxygen vacancies
        if oxygen_vacancy_fraction > 0 and sites["O_site"]:
            n_remove = max(1, int(len(sites["O_site"]) * oxygen_vacancy_fraction))
            selected = np.random.choice(
                sites["O_site"], size=min(n_remove, len(sites["O_site"])), replace=False
            )
            vacancy_indices.extend(selected.tolist())

        # Remove atoms
        vacancy_indices = list(set(vacancy_indices))
        keep_indices = [i for i in range(len(self.atoms)) if i not in vacancy_indices]
        new_atoms = self.atoms[keep_indices]

        return new_atoms

    def create_surface_segregation(
        self,
        metal: str,
        n_atoms: int = 1,
        random_seed: Optional[int] = None,
    ) -> Atoms:
        """
        Move B-site metal atoms from bulk layer to surface layer.

        Models the initial segregation step of exsolution.

        Parameters
        ----------
        metal : str
            Element to segregate (e.g., "Ni", "Co", "Fe")
        n_atoms : int
            Number of atoms to move to surface
        random_seed : int, optional
            Random seed for reproducibility

        Returns
        -------
        Atoms
            Structure with segregated metal at surface
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        new_atoms = self.atoms.copy()
        symbols = list(new_atoms.get_chemical_symbols())
        positions = new_atoms.get_positions().copy()

        # Find metal atoms
        metal_indices = self.get_element_indices(metal)
        if not metal_indices:
            raise ValueError(f"No {metal} atoms found in structure")

        # Get z positions
        z_positions = positions[:, 2]
        z_max = z_positions.max()
        z_min = z_positions.min()
        z_mid = (z_max + z_min) / 2

        # Find bulk metal atoms (below mid-height)
        bulk_metal_indices = [i for i in metal_indices if positions[i, 2] < z_mid]
        if not bulk_metal_indices:
            bulk_metal_indices = metal_indices

        # Select atoms to move
        n_to_move = min(n_atoms, len(bulk_metal_indices))
        selected = np.random.choice(bulk_metal_indices, size=n_to_move, replace=False)

        # Move selected atoms to surface region
        surface_z = z_max + 1.5  # Place slightly above current surface
        for idx in selected:
            # Keep x, y position, move z to surface
            positions[idx, 2] = surface_z
            surface_z += 0.5  # Stagger multiple atoms slightly

        new_atoms.set_positions(positions)
        return new_atoms

    def create_nanoparticle(
        self,
        metal: str,
        n_atoms: int,
        shape: str = "hemispherical",
        position: str = "hollow",
        interface_distance: float = 2.0,
        socketed: bool = True,
        random_seed: Optional[int] = None,
    ) -> Atoms:
        """
        Place metallic nanoparticle on surface.

        Parameters
        ----------
        metal : str
            Metal element ("Ni", "Co", "Fe")
        n_atoms : int
            Number of atoms in cluster (1, 4, 7, 13, 19 for hemispherical)
        shape : str
            Cluster shape ("hemispherical" or "icosahedral")
        position : str
            Placement position ("hollow", "ontop", "bridge", "random")
        interface_distance : float
            Distance between particle bottom and surface (A)
        socketed : bool
            If True, remove surface atoms under particle (models real exsolution)
        random_seed : int, optional
            Random seed for reproducibility

        Returns
        -------
        Atoms
            Structure with nanoparticle
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        # Create metal cluster
        cluster = create_metallic_cluster(metal, n_atoms, shape)

        # Get surface information
        surface_indices, z_cutoff = get_surface_atoms(self.atoms, z_threshold=0.15)
        z_positions = self.atoms.positions[:, 2]
        surface_z = z_positions[surface_indices].max() if surface_indices else z_positions.max()

        # Get cell dimensions
        cell = self.atoms.get_cell()
        cell_x, cell_y = cell[0, 0], cell[1, 1]

        # Determine placement position
        if position == "random":
            x = np.random.uniform(cell_x * 0.2, cell_x * 0.8)
            y = np.random.uniform(cell_y * 0.2, cell_y * 0.8)
        elif position == "hollow":
            # Place at cell center
            x = cell_x / 2
            y = cell_y / 2
        elif position == "ontop":
            # Place above a B-site atom
            sites = self.identify_perovskite_sites()
            if sites["B_site"]:
                b_surface = [i for i in sites["B_site"] if z_positions[i] > z_cutoff]
                if b_surface:
                    idx = np.random.choice(b_surface)
                    x, y = self.atoms.positions[idx, :2]
                else:
                    x, y = cell_x / 2, cell_y / 2
            else:
                x, y = cell_x / 2, cell_y / 2
        else:  # bridge
            x = cell_x / 2 + cell_x * 0.1
            y = cell_y / 2

        # Position cluster
        cluster_positions = cluster.get_positions()
        cluster_center = cluster_positions.mean(axis=0)
        cluster_bottom = cluster_positions[:, 2].min()

        # Shift cluster to target position
        shift = np.array([
            x - cluster_center[0],
            y - cluster_center[1],
            surface_z + interface_distance - cluster_bottom
        ])
        cluster.set_positions(cluster_positions + shift)

        # Handle socketing (remove surface atoms under particle)
        new_slab = self.atoms.copy()
        if socketed and n_atoms > 1:
            # Find surface atoms near particle center
            particle_xy = np.array([x, y])
            socket_radius = EXSOLUTION_METALS.get(metal, {}).get("radius", 1.25) * 1.5

            atoms_to_remove = []
            for i in surface_indices:
                atom_xy = self.atoms.positions[i, :2]
                dist = np.linalg.norm(atom_xy - particle_xy)
                if dist < socket_radius:
                    atoms_to_remove.append(i)

            if atoms_to_remove:
                keep_indices = [i for i in range(len(new_slab)) if i not in atoms_to_remove]
                new_slab = new_slab[keep_indices]

        # Combine slab and cluster
        combined = new_slab + cluster

        return combined

    def create_exsolution_pathway(
        self,
        metal: str,
        particle_size: int,
        vacancy_fraction: float = 0.1,
        random_seed: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        """
        Generate complete exsolution pathway structures.

        Creates 4 stages:
        1. Pristine perovskite surface
        2. Defective perovskite (with vacancies)
        3. Surface-segregated metal
        4. Exsolved nanoparticle

        Parameters
        ----------
        metal : str
            Exsolution metal (e.g., "Ni")
        particle_size : int
            Number of atoms in final nanoparticle
        vacancy_fraction : float
            Vacancy concentration (0.0 to 1.0)
        random_seed : int, optional
            Random seed for reproducibility

        Returns
        -------
        list
            List of dicts with 'atoms', 'stage', 'description'
        """
        pathway = []
        seed = random_seed if random_seed is not None else 42

        # Stage 1: Pristine
        pathway.append({
            "atoms": self.atoms.copy(),
            "stage": "pristine",
            "description": "Pristine perovskite surface",
        })

        # Stage 2: Defective
        # Calculate B-site vacancy fraction based on particle size
        sites = self.identify_perovskite_sites()
        metal_indices = self.get_element_indices(metal)
        n_metal = len(metal_indices)
        b_site_vac_frac = min(particle_size / max(n_metal, 1), 0.5)

        defective = self.create_defective_perovskite(
            a_site_vacancy_fraction=vacancy_fraction * 0.5,
            b_site_vacancy_fraction=b_site_vac_frac,
            oxygen_vacancy_fraction=vacancy_fraction,
            b_site_element=metal,
            random_seed=seed,
        )
        pathway.append({
            "atoms": defective,
            "stage": "defective",
            "description": f"Defective perovskite (O vac: {vacancy_fraction:.0%})",
        })

        # Stage 3: Segregated
        segregated = self.create_surface_segregation(
            metal=metal,
            n_atoms=min(particle_size, 3),
            random_seed=seed + 1,
        )
        pathway.append({
            "atoms": segregated,
            "stage": "segregated",
            "description": f"{metal} segregated to surface",
        })

        # Stage 4: Exsolved nanoparticle
        exsolved = self.create_nanoparticle(
            metal=metal,
            n_atoms=particle_size,
            shape="hemispherical",
            socketed=True,
            random_seed=seed + 2,
        )
        pathway.append({
            "atoms": exsolved,
            "stage": "exsolved",
            "description": f"Exsolved {metal}{particle_size} nanoparticle",
        })

        return pathway

    def get_adsorption_sites(
        self,
        structure: Optional[Atoms] = None,
        particle_indices: Optional[List[int]] = None,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """
        Identify adsorption sites on exsolved particle system.

        Returns different site types:
        - metal_top: Top of metal particle atoms
        - interface_edge: Metal-oxide boundary
        - vacancy_site: Near oxygen vacancies
        - oxide_surface: Remaining perovskite surface

        Parameters
        ----------
        structure : Atoms, optional
            Structure with exsolved particle (uses self.atoms if None)
        particle_indices : list, optional
            Indices of particle atoms (auto-detected if None)

        Returns
        -------
        dict
            Site types with positions and indices
        """
        atoms = structure if structure is not None else self.atoms
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()

        # Auto-detect particle (metal atoms at top of structure)
        if particle_indices is None:
            metal_elements = list(EXSOLUTION_METALS.keys())
            z_positions = positions[:, 2]
            z_threshold = z_positions.max() - 3.0  # Top 3 A

            particle_indices = [
                i for i, (s, z) in enumerate(zip(symbols, z_positions))
                if s in metal_elements and z > z_threshold
            ]

        sites = {
            "metal_top": [],
            "interface_edge": [],
            "vacancy_site": [],
            "oxide_surface": [],
        }

        # Metal top sites
        for i in particle_indices:
            sites["metal_top"].append({
                "index": i,
                "position": positions[i].tolist(),
                "element": symbols[i],
            })

        # Identify surface oxide atoms
        surface_indices, z_cutoff = get_surface_atoms(atoms, z_threshold=0.2)

        # Oxide surface sites (non-particle surface atoms)
        for i in surface_indices:
            if i not in particle_indices:
                sites["oxide_surface"].append({
                    "index": i,
                    "position": positions[i].tolist(),
                    "element": symbols[i],
                })

        # Interface edge sites (between particle and oxide)
        if particle_indices:
            particle_center = positions[particle_indices].mean(axis=0)[:2]
            particle_radius = 2.5  # A

            for i in surface_indices:
                if i not in particle_indices:
                    atom_xy = positions[i, :2]
                    dist = np.linalg.norm(atom_xy - particle_center)
                    if dist < particle_radius * 1.5:
                        sites["interface_edge"].append({
                            "index": i,
                            "position": positions[i].tolist(),
                            "element": symbols[i],
                            "distance_to_particle": float(dist),
                        })

        return sites


class ExsolutionStructure(BaseStructure):
    """
    Structure class for exsolution systems with metadata tracking.

    Tracks vacancy concentrations, nanoparticle info, and exsolution state.
    """

    def __init__(
        self,
        atoms: Atoms,
        exsolution_metal: Optional[str] = None,
        nanoparticle_size: Optional[int] = None,
        a_site_vacancy_concentration: float = 0.0,
        b_site_vacancy_concentration: float = 0.0,
        oxygen_vacancy_concentration: float = 0.0,
        exsolution_state: str = "pristine",
        parent_formula: Optional[str] = None,
    ):
        """
        Initialize ExsolutionStructure.

        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object
        exsolution_metal : str, optional
            Metal being exsolved (Ni, Co, Fe)
        nanoparticle_size : int, optional
            Number of atoms in nanoparticle
        a_site_vacancy_concentration : float
            A-site vacancy fraction
        b_site_vacancy_concentration : float
            B-site vacancy fraction
        oxygen_vacancy_concentration : float
            Oxygen vacancy fraction
        exsolution_state : str
            State: "pristine", "defective", "segregated", "exsolved"
        parent_formula : str, optional
            Original perovskite formula
        """
        super().__init__(atoms)
        self.exsolution_metal = exsolution_metal
        self.nanoparticle_size = nanoparticle_size
        self.a_site_vacancy_concentration = a_site_vacancy_concentration
        self.b_site_vacancy_concentration = b_site_vacancy_concentration
        self.oxygen_vacancy_concentration = oxygen_vacancy_concentration
        self.exsolution_state = exsolution_state
        self.parent_formula = parent_formula

    @classmethod
    def from_pathway_step(
        cls,
        pathway_step: Dict[str, Any],
        metal: str,
        vacancy_fraction: float = 0.0,
    ) -> "ExsolutionStructure":
        """
        Create ExsolutionStructure from pathway step dict.

        Parameters
        ----------
        pathway_step : dict
            Dict with 'atoms', 'stage', 'description' keys
        metal : str
            Exsolution metal
        vacancy_fraction : float
            Vacancy concentration

        Returns
        -------
        ExsolutionStructure
            Structure with metadata
        """
        atoms = pathway_step["atoms"]
        stage = pathway_step["stage"]

        return cls(
            atoms=atoms,
            exsolution_metal=metal,
            exsolution_state=stage,
            oxygen_vacancy_concentration=vacancy_fraction if stage != "pristine" else 0.0,
        )

    def get_defect_info(self) -> Dict[str, Any]:
        """Get information about defects and composition."""
        counts = self.get_element_count()

        return {
            "formula": self.formula,
            "n_atoms": self.n_atoms,
            "element_counts": counts,
            "exsolution_metal": self.exsolution_metal,
            "nanoparticle_size": self.nanoparticle_size,
            "exsolution_state": self.exsolution_state,
            "vacancy_concentrations": {
                "A_site": self.a_site_vacancy_concentration,
                "B_site": self.b_site_vacancy_concentration,
                "O": self.oxygen_vacancy_concentration,
            },
        }

    def __repr__(self) -> str:
        return (
            f"ExsolutionStructure({self.formula}, "
            f"metal={self.exsolution_metal}, "
            f"state={self.exsolution_state})"
        )


def create_metallic_cluster(
    metal: str,
    n_atoms: int,
    shape: str = "hemispherical",
) -> Atoms:
    """
    Create isolated metallic cluster of specified size and shape.

    Parameters
    ----------
    metal : str
        Metal element (Ni, Co, Fe)
    n_atoms : int
        Number of atoms in cluster
    shape : str
        Cluster shape ("hemispherical", "icosahedral")

    Returns
    -------
    Atoms
        Metal cluster
    """
    if n_atoms == 1:
        return Atoms(metal, positions=[[0, 0, 0]])

    # Get metal properties
    metal_info = EXSOLUTION_METALS.get(metal, {"radius": 1.25})
    bond_length = metal_info["radius"] * 2 * 0.9  # Slightly compressed for cluster

    if shape == "icosahedral" and n_atoms in [13, 55, 147]:
        # Build icosahedral cluster
        from ase.cluster import Icosahedron
        try:
            # Determine shell number
            if n_atoms == 13:
                nshells = 1
            elif n_atoms == 55:
                nshells = 2
            else:
                nshells = 3
            cluster = Icosahedron(metal, nshells, latticeconstant=bond_length * 1.4)
        except Exception:
            # Fallback to simple cluster
            cluster = _build_simple_cluster(metal, n_atoms, bond_length)
    else:
        # Build hemispherical cluster
        cluster = _build_hemispherical_cluster(metal, n_atoms, bond_length)

    return cluster


def _build_hemispherical_cluster(metal: str, n_atoms: int, bond_length: float) -> Atoms:
    """Build hemispherical cluster for surface deposition."""
    positions = []

    if n_atoms >= 1:
        positions.append([0, 0, 0])  # Center atom

    if n_atoms >= 4:
        # First shell (triangle around center)
        for i in range(3):
            angle = i * 2 * np.pi / 3
            x = bond_length * np.cos(angle)
            y = bond_length * np.sin(angle)
            positions.append([x, y, 0])

    if n_atoms >= 7:
        # Second layer above
        z_up = bond_length * 0.8
        positions.append([0, 0, z_up])
        for i in range(2):
            angle = i * np.pi + np.pi / 6
            x = bond_length * 0.7 * np.cos(angle)
            y = bond_length * 0.7 * np.sin(angle)
            positions.append([x, y, z_up])

    if n_atoms >= 10:
        # Outer shell
        for i in range(3):
            angle = i * 2 * np.pi / 3 + np.pi / 3
            x = bond_length * 1.5 * np.cos(angle)
            y = bond_length * 1.5 * np.sin(angle)
            positions.append([x, y, 0])

    if n_atoms >= 13:
        # Top atoms
        z_top = bond_length * 1.5
        positions.append([0, 0, z_top])
        for i in range(2):
            angle = i * np.pi
            x = bond_length * 0.5 * np.cos(angle)
            y = bond_length * 0.5 * np.sin(angle)
            positions.append([x, y, z_top * 0.8])

    if n_atoms >= 19:
        # Additional outer atoms
        for i in range(6):
            angle = i * np.pi / 3
            x = bond_length * 2.0 * np.cos(angle)
            y = bond_length * 2.0 * np.sin(angle)
            positions.append([x, y, 0])

    positions = positions[:n_atoms]
    return Atoms(metal * len(positions), positions=positions)


def _build_simple_cluster(metal: str, n_atoms: int, bond_length: float) -> Atoms:
    """Build simple compact cluster."""
    positions = [[0, 0, 0]]

    # Add atoms in shells
    r = bond_length
    idx = 1
    layer = 0

    while idx < n_atoms:
        n_in_ring = max(3, layer * 3 + 3)
        for i in range(n_in_ring):
            if idx >= n_atoms:
                break
            angle = i * 2 * np.pi / n_in_ring
            x = r * np.cos(angle)
            y = r * np.sin(angle)
            z = (layer % 2) * bond_length * 0.5
            positions.append([x, y, z])
            idx += 1
        r += bond_length * 0.8
        layer += 1

    return Atoms(metal * len(positions), positions=positions[:n_atoms])


def generate_exsolution_series(
    base_structure: Union[Atoms, BaseStructure],
    metal: str,
    vacancy_range: List[float],
    particle_sizes: List[int],
    n_configs: int = 3,
    random_seed: Optional[int] = None,
) -> List[Dict[str, Any]]:
    """
    Generate series of exsolution structures for screening.

    Parameters
    ----------
    base_structure : Atoms or BaseStructure
        Base perovskite surface structure
    metal : str
        Exsolution metal (Ni, Co, Fe)
    vacancy_range : list
        List of vacancy concentrations to test
    particle_sizes : list
        List of nanoparticle sizes to test
    n_configs : int
        Number of random configurations per combination
    random_seed : int, optional
        Base random seed

    Returns
    -------
    list
        List of dicts with structure and metadata
    """
    builder = ExsolutionBuilder(base_structure)
    results = []
    config_id = 0

    for vac_conc in vacancy_range:
        for size in particle_sizes:
            for i in range(n_configs):
                seed = None
                if random_seed is not None:
                    seed = random_seed + config_id

                # Create exsolved structure
                atoms = builder.create_nanoparticle(
                    metal=metal,
                    n_atoms=size,
                    socketed=True,
                    random_seed=seed,
                )

                # Also create defective substrate version
                defective_builder = ExsolutionBuilder(base_structure)
                defective = defective_builder.create_defective_perovskite(
                    oxygen_vacancy_fraction=vac_conc,
                    random_seed=seed,
                )

                results.append({
                    "atoms": atoms,
                    "defective_substrate": defective,
                    "metal": metal,
                    "particle_size": size,
                    "vacancy_concentration": vac_conc,
                    "config_id": config_id,
                })
                config_id += 1

    return results
