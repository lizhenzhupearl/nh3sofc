"""Adsorbate placement on surfaces.

Provides 6 methods for placing adsorbates:
1. add_simple() - Manual (x,y) + height using add_adsorbate()
2. add_random() - Random placement with Euler angle rotation
3. add_grid() - Grid-based systematic placement
4. add_on_site() - Site-specific (above specific atom types)
5. add_with_collision() - Random with minimum distance check
6. add_catkit() - CatKit integration (top/bridge/hollow sites)
"""

from pathlib import Path
from typing import Optional, List, Union, Tuple, Dict, Any
import numpy as np
from ase import Atoms
from ase.build import add_adsorbate, molecule
from ase.io import write

from ..core.base import BaseStructure, get_surface_atoms
from ..core.constants import ADSORBATES
from .surface import SlabStructure


class AdsorbatePlacer:
    """
    Class for placing adsorbates on surfaces.

    Provides multiple methods for adsorbate placement with different
    levels of control and automation.

    Examples
    --------
    >>> slab = SlabStructure.from_file("slab.cif")
    >>> placer = AdsorbatePlacer(slab)
    >>>
    >>> # Simple placement
    >>> config = placer.add_simple("NH3", position=(5.0, 5.0), height=2.0)
    >>>
    >>> # Site-specific placement
    >>> configs = placer.add_on_site("NH3", atom_types=["La", "V"])
    >>>
    >>> # Random placement with collision detection
    >>> config = placer.add_with_collision("NH3", min_distance=2.0)
    """

    def __init__(self, slab: Union[SlabStructure, Atoms]):
        """
        Initialize AdsorbatePlacer.

        Parameters
        ----------
        slab : SlabStructure or Atoms
            Surface slab to place adsorbates on
        """
        if isinstance(slab, SlabStructure):
            self.slab = slab.atoms.copy()
            self._slab_obj = slab
        else:
            self.slab = slab.copy()
            self._slab_obj = None

        self.n_slab_atoms = len(self.slab)

    def _get_molecule(self, name: str) -> Atoms:
        """
        Get adsorbate molecule.

        Parameters
        ----------
        name : str
            Molecule name (e.g., "NH3", "H2O")

        Returns
        -------
        Atoms
            Molecule as ASE Atoms object
        """
        try:
            mol = molecule(name)
        except KeyError:
            # Try common variations
            name_map = {"NH2": "NH2", "NH": "NH", "OH": "OH"}
            if name in name_map:
                mol = molecule(name_map[name])
            else:
                raise ValueError(f"Unknown molecule: {name}")
        return mol

    def _get_surface_z(self) -> float:
        """Get the z-coordinate of the surface (highest z)."""
        return self.slab.positions[:, 2].max()

    def _apply_random_rotation(self, mol: Atoms, uniform_sampling: bool = True) -> None:
        """
        Apply random 3D rotation to a molecule.

        Parameters
        ----------
        mol : Atoms
            Molecule to rotate (modified in place)
        uniform_sampling : bool
            If True, use uniform sampling on SO(3) for unbiased orientations.
            If False, use simple uniform Euler angles (biased toward poles).
        """
        alpha = np.random.uniform(0, 2 * np.pi)  # z rotation
        gamma = np.random.uniform(0, 2 * np.pi)  # x rotation

        if uniform_sampling:
            # Uniform sampling on SO(3): beta = arccos(1 - 2u) for u ~ U(0,1)
            # This gives uniform distribution on the sphere
            beta = np.arccos(1 - 2 * np.random.uniform(0, 1))
        else:
            beta = np.random.uniform(0, np.pi)  # y rotation (biased)

        mol.rotate(np.degrees(alpha), "z")
        mol.rotate(np.degrees(beta), "y")
        mol.rotate(np.degrees(gamma), "x")

    def _get_cell_dimensions(self) -> Tuple[float, float]:
        """Get x and y dimensions of the cell."""
        cell = self.slab.get_cell()
        return cell[0, 0], cell[1, 1]

    def add_simple(
        self,
        adsorbate: str,
        position: Tuple[float, float],
        height: float = 2.0,
        mol_index: int = 0,
        rotation: Optional[Tuple[float, float, float]] = None,
        orientation: Optional[Tuple[float, float, float]] = None,
    ) -> Atoms:
        """
        Add adsorbate at a specific position (Method 1).

        Uses ASE's add_adsorbate function for simple manual placement.

        Parameters
        ----------
        adsorbate : str
            Adsorbate molecule name (e.g., "NH3")
        position : tuple
            (x, y) position on the surface
        height : float
            Height above the surface (Angstrom)
        mol_index : int
            Index of atom in molecule to use as anchor
        rotation : tuple, optional
            Euler angles (alpha, beta, gamma) in radians
        orientation : tuple, optional
            Alias for rotation (deprecated, use rotation instead)

        Returns
        -------
        Atoms
            Slab with adsorbate
        """
        # Handle orientation alias (deprecated)
        if orientation is not None and rotation is None:
            rotation = orientation

        slab_copy = self.slab.copy()
        mol = self._get_molecule(adsorbate)

        # Apply rotation if specified
        if rotation is not None:
            alpha, beta, gamma = rotation
            mol.rotate(np.degrees(alpha), "z")
            mol.rotate(np.degrees(beta), "y")
            mol.rotate(np.degrees(gamma), "x")

        add_adsorbate(
            slab_copy,
            mol,
            height=height,
            position=position,
            mol_index=mol_index
        )

        return slab_copy

    def add_random(
        self,
        adsorbate: str,
        n_configs: int = 10,
        height: float = 2.0,
        random_seed: Optional[int] = None,
        seed: Optional[int] = None,
    ) -> List[Atoms]:
        """
        Add adsorbate at random positions with random orientations (Method 2).

        Parameters
        ----------
        adsorbate : str
            Adsorbate molecule name
        n_configs : int
            Number of random configurations to generate
        height : float
            Height above surface
        random_seed : int, optional
            Random seed for reproducibility
        seed : int, optional
            Alias for random_seed (deprecated, use random_seed instead)

        Returns
        -------
        list
            List of Atoms objects with different configurations
        """
        # Handle seed alias (deprecated)
        if seed is not None and random_seed is None:
            random_seed = seed

        if random_seed is not None:
            np.random.seed(random_seed)

        a, b = self._get_cell_dimensions()
        z_top = self._get_surface_z()

        configs = []

        for _ in range(n_configs):
            slab_copy = self.slab.copy()
            mol = self._get_molecule(adsorbate)

            # Random position
            x = np.random.uniform(0, a)
            y = np.random.uniform(0, b)
            z = z_top + height

            # Apply uniform random rotation
            self._apply_random_rotation(mol)

            # Translate to position
            mol.translate([x, y, z])

            # Combine
            combined = slab_copy + mol
            configs.append(combined)

        return configs

    def add_grid(
        self,
        adsorbate: str,
        grid_size: Tuple[int, int] = (3, 3),
        orientations: int = 4,
        height: float = 2.0,
        random_seed: Optional[int] = None,
        nx: Optional[int] = None,
        ny: Optional[int] = None,
    ) -> List[Atoms]:
        """
        Add adsorbate on a grid of positions (Method 3).

        Parameters
        ----------
        adsorbate : str
            Adsorbate molecule name
        grid_size : tuple
            Grid dimensions (nx, ny)
        orientations : int
            Number of random orientations per grid point
        height : float
            Height above surface
        random_seed : int, optional
            Random seed
        nx : int, optional
            Grid x dimension (deprecated, use grid_size instead)
        ny : int, optional
            Grid y dimension (deprecated, use grid_size instead)

        Returns
        -------
        list
            List of Atoms objects for all grid positions
        """
        # Handle nx/ny aliases (deprecated)
        if nx is not None or ny is not None:
            grid_size = (nx if nx is not None else grid_size[0],
                        ny if ny is not None else grid_size[1])

        if random_seed is not None:
            np.random.seed(random_seed)

        a, b = self._get_cell_dimensions()
        z_top = self._get_surface_z()

        nx, ny = grid_size
        x_grid = np.linspace(0, a, nx, endpoint=False) + a / (2 * nx)
        y_grid = np.linspace(0, b, ny, endpoint=False) + b / (2 * ny)

        configs = []

        for x in x_grid:
            for y in y_grid:
                for _ in range(orientations):
                    slab_copy = self.slab.copy()
                    mol = self._get_molecule(adsorbate)

                    # Apply uniform random rotation
                    self._apply_random_rotation(mol)

                    # Position
                    mol.translate([x, y, z_top + height])

                    combined = slab_copy + mol
                    configs.append(combined)

        return configs

    def add_on_site(
        self,
        adsorbate: str,
        site_type: str = "ontop",
        atom_types: Optional[List[str]] = None,
        site_indices: Optional[List[int]] = None,
        n_orientations: int = 5,
        height: float = 2.0,
        layer_tolerance: float = 2.0,
        random_seed: Optional[int] = None,
    ) -> List[Atoms]:
        """
        Add adsorbate above specific surface atoms (Method 4).

        This is the most physically meaningful placement method,
        targeting specific adsorption sites. The selection process is:

        1. First, identify surface atoms (within layer_tolerance of topmost atom)
        2. Then, filter by atom_types if specified (e.g., only La, V on surface)
        3. Or, filter by site_indices if specified (validated against surface)

        This ensures adsorbates are only placed on atoms that are actually
        exposed at the surface, not buried in the bulk.

        Parameters
        ----------
        adsorbate : str
            Adsorbate molecule name
        site_type : str
            Site type: "ontop" (above atoms), "bridge" (between pairs),
            or "hollow" (between triplets). Bridge and hollow require
            CatKit for accurate site detection.
        atom_types : list, optional
            List of element symbols to target (e.g., ["La", "V"]).
            Only atoms of these types that are ON THE SURFACE will be used.
        site_indices : list, optional
            Specific atom indices to use as sites. These are validated
            against surface atoms; non-surface indices are filtered out
            with a warning.
        n_orientations : int
            Number of orientations per site
        height : float
            Height above surface
        layer_tolerance : float
            Distance (Angstrom) from topmost atom to consider as surface.
            Default 2.0 Å captures the top bilayer in perovskites.
        random_seed : int, optional
            Random seed

        Returns
        -------
        list
            List of Atoms objects for each site/orientation

        Examples
        --------
        >>> # Place NH3 on La atoms in the top bilayer
        >>> configs = placer.add_on_site("NH3", atom_types=["La"])
        >>> # La atoms in the bulk are automatically excluded
        >>>
        >>> # Use smaller tolerance for only the topmost layer
        >>> configs = placer.add_on_site("NH3", atom_types=["V"], layer_tolerance=1.0)
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        # Handle bridge and hollow sites
        site_type_lower = site_type.lower()
        if site_type_lower in ("bridge", "hollow", "4fold"):
            # Delegate to CatKit for accurate site detection
            try:
                return self.add_catkit(adsorbate, site_type=site_type, height=height)
            except ImportError:
                # Fallback: compute simple bridge/hollow sites
                return self._add_on_computed_site(
                    adsorbate, site_type_lower, atom_types,
                    n_orientations, height, layer_tolerance
                )

        # Ontop placement: directly above surface atoms
        # Get surface atoms (within layer_tolerance of topmost atom)
        surface_indices, z_cutoff = get_surface_atoms(
            self.slab, layer_tolerance=layer_tolerance
        )
        surface_set = set(surface_indices)

        # Filter by atom type if specified
        if atom_types is not None:
            surface_indices = [
                i for i in surface_indices
                if self.slab[i].symbol in atom_types
            ]

        # Filter by specific indices if provided
        if site_indices is not None:
            # Validate that specified indices are on the surface
            non_surface = [i for i in site_indices if i not in surface_set]
            if non_surface:
                z_positions = self.slab.positions[non_surface, 2]
                print(f"Warning: site_indices {non_surface} are not on the surface "
                      f"(z < {z_cutoff:.2f}). Their z positions: {z_positions}. "
                      f"Filtering to surface atoms only.")
            # Intersect with surface atoms
            surface_indices = [i for i in site_indices if i in surface_set]

        if not surface_indices:
            raise ValueError("No suitable adsorption sites found. "
                           "Check that atom_types or site_indices are on the surface.")

        configs = []

        for idx in surface_indices:
            site_pos = self.slab.positions[idx]

            for _ in range(n_orientations):
                slab_copy = self.slab.copy()
                mol = self._get_molecule(adsorbate)

                # Apply uniform random rotation
                self._apply_random_rotation(mol)

                # Place above the site
                mol.translate([site_pos[0], site_pos[1], site_pos[2] + height])

                combined = slab_copy + mol
                configs.append(combined)

        return configs

    def _add_on_computed_site(
        self,
        adsorbate: str,
        site_type: str,
        atom_types: Optional[List[str]],
        n_orientations: int,
        height: float,
        layer_tolerance: float,
    ) -> List[Atoms]:
        """
        Fallback for bridge/hollow sites when CatKit is not available.

        Computes approximate site positions from surface atom geometry.
        """
        from scipy.spatial import Delaunay

        surface_indices, _ = get_surface_atoms(self.slab, layer_tolerance=layer_tolerance)

        # Filter by atom type if specified
        if atom_types is not None:
            surface_indices = [
                i for i in surface_indices
                if self.slab[i].symbol in atom_types
            ]

        if len(surface_indices) < 2:
            raise ValueError("Not enough surface atoms for bridge/hollow sites")

        surface_positions = self.slab.positions[surface_indices]
        z_top = surface_positions[:, 2].max()

        # Get xy positions for triangulation
        xy_positions = surface_positions[:, :2]

        sites = []

        if site_type == "bridge":
            # Bridge sites: midpoints between neighboring atoms
            # Use distance cutoff to find neighbors
            cutoff = 4.0  # Angstrom
            for i, pos_i in enumerate(surface_positions):
                for j, pos_j in enumerate(surface_positions):
                    if j <= i:
                        continue
                    dist = np.linalg.norm(pos_i[:2] - pos_j[:2])
                    if dist < cutoff:
                        midpoint = (pos_i + pos_j) / 2
                        sites.append(midpoint[:2])

        elif site_type == "hollow":
            # Hollow sites: centroids of triangles from Delaunay triangulation
            if len(xy_positions) >= 3:
                try:
                    tri = Delaunay(xy_positions)
                    for simplex in tri.simplices:
                        centroid = xy_positions[simplex].mean(axis=0)
                        sites.append(centroid)
                except Exception:
                    # Fallback: use centroids of close triplets
                    pass

        if not sites:
            raise ValueError(f"Could not find {site_type} sites. "
                           "Consider installing CatKit for better site detection.")

        configs = []

        for site_xy in sites:
            for _ in range(n_orientations):
                slab_copy = self.slab.copy()
                mol = self._get_molecule(adsorbate)

                # Apply uniform random rotation
                self._apply_random_rotation(mol)

                # Place above the site
                mol.translate([site_xy[0], site_xy[1], z_top + height])

                combined = slab_copy + mol
                configs.append(combined)

        return configs

    def add_with_collision(
        self,
        adsorbate: str,
        n_configs: int = 1,
        min_distance: float = 2.0,
        height: float = 2.0,
        max_attempts: int = 100,
        random_seed: Optional[int] = None,
    ) -> Union[Optional[Atoms], List[Atoms]]:
        """
        Add adsorbate with collision detection (Method 5).

        Ensures no atoms are closer than min_distance.

        Parameters
        ----------
        adsorbate : str
            Adsorbate molecule name
        n_configs : int
            Number of configurations to generate. If 1, returns single
            Atoms or None (backward compatible). If >1, returns list.
        min_distance : float
            Minimum allowed distance between atoms (Angstrom)
        height : float
            Height above surface
        max_attempts : int
            Maximum placement attempts per configuration
        random_seed : int, optional
            Random seed

        Returns
        -------
        Atoms, None, or list
            If n_configs=1: Slab with adsorbate, or None if placement failed
            If n_configs>1: List of Atoms objects (may be shorter than n_configs)
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        a, b = self._get_cell_dimensions()
        z_top = self._get_surface_z()

        configs = []
        total_attempts = 0
        max_total_attempts = max_attempts * n_configs

        while len(configs) < n_configs and total_attempts < max_total_attempts:
            total_attempts += 1

            # Random position
            x = np.random.uniform(0, a)
            y = np.random.uniform(0, b)
            z = z_top + height

            mol = self._get_molecule(adsorbate)

            # Apply uniform random rotation
            self._apply_random_rotation(mol)

            mol.translate([x, y, z])

            # Check distances
            min_dist = float("inf")
            for slab_atom in self.slab:
                for mol_atom in mol:
                    dist = np.linalg.norm(
                        slab_atom.position - mol_atom.position
                    )
                    min_dist = min(min_dist, dist)

            if min_dist > min_distance:
                combined = self.slab.copy() + mol
                configs.append(combined)

        if len(configs) < n_configs:
            print(f"Warning: Only generated {len(configs)}/{n_configs} configs "
                  f"after {total_attempts} attempts")

        # Backward compatibility: return single Atoms or None if n_configs=1
        if n_configs == 1:
            return configs[0] if configs else None

        return configs

    def add_catkit(
        self,
        adsorbate: str,
        site_type: str = "ontop",
        height: float = 2.0,
        n_sites: Optional[int] = None,
        n_orientations: int = 1,
        random_seed: Optional[int] = None,
    ) -> List[Atoms]:
        """
        Add adsorbate using CatKit for site detection (Method 6).

        Uses CatKit's sophisticated site-finding algorithms for
        top, bridge, hollow, and 4-fold sites.

        Parameters
        ----------
        adsorbate : str
            Adsorbate molecule name
        site_type : str
            Site type: "ontop", "bridge", "hollow", or "4fold"
        height : float
            Height above surface
        n_sites : int, optional
            Maximum number of sites to use
        n_orientations : int
            Number of random orientations per site
        random_seed : int, optional
            Random seed for reproducibility

        Returns
        -------
        list
            List of Atoms objects for each site

        Raises
        ------
        ImportError
            If CatKit is not installed
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        try:
            from catkit.gen.adsorption import AdsorptionSites
        except ImportError:
            raise ImportError(
                "CatKit is required for this method. "
                "Install with: pip install catkit"
            )

        # Get adsorption sites
        sites = AdsorptionSites(self.slab)

        site_map = {
            "ontop": "top",
            "top": "top",
            "bridge": "bridge",
            "hollow": "hollow",
            "4fold": "4fold",
        }

        site_key = site_map.get(site_type.lower(), site_type)
        coordinates = sites.get_coordinates(site_type=site_key)

        if n_sites is not None:
            coordinates = coordinates[:n_sites]

        configs = []
        z_top = self._get_surface_z()

        for coord in coordinates:
            for _ in range(n_orientations):
                slab_copy = self.slab.copy()
                mol = self._get_molecule(adsorbate)

                # Apply uniform random rotation
                self._apply_random_rotation(mol)

                # Place at site
                mol.translate([coord[0], coord[1], z_top + height])

                combined = slab_copy + mol
                configs.append(combined)

        return configs

    def add_at_all_sites(
        self,
        adsorbate: str,
        site_type: str = "ontop",
        atom_types: Optional[List[str]] = None,
        height: float = 2.0,
        layer_tolerance: float = 2.0,
    ) -> List[Atoms]:
        """
        Generate configurations for all sites of a given type.

        Convenience method that combines add_on_site with single orientation.

        Parameters
        ----------
        adsorbate : str
            Adsorbate molecule name
        site_type : str
            Site type ("ontop" uses surface atoms)
        atom_types : list, optional
            Element types to use as sites
        height : float
            Height above surface
        layer_tolerance : float
            Distance (Angstrom) from topmost atom to consider as surface

        Returns
        -------
        list
            List of Atoms objects for each site
        """
        if site_type == "ontop":
            return self.add_on_site(
                adsorbate,
                atom_types=atom_types,
                n_orientations=1,
                height=height,
                layer_tolerance=layer_tolerance,
            )
        else:
            # Try CatKit for other site types
            try:
                return self.add_catkit(adsorbate, site_type=site_type, height=height)
            except ImportError:
                raise ValueError(
                    f"Site type '{site_type}' requires CatKit. "
                    "Install with: pip install catkit"
                )

    def get_site_info(self, layer_tolerance: float = 2.0) -> Dict[str, Any]:
        """
        Get information about available adsorption sites.

        Parameters
        ----------
        layer_tolerance : float
            Distance (Angstrom) from topmost atom to consider as surface

        Returns
        -------
        dict
            Site information including counts by element
        """
        surface_indices, z_cutoff = get_surface_atoms(
            self.slab, layer_tolerance=layer_tolerance
        )

        # Count by element
        element_counts = {}
        for idx in surface_indices:
            symbol = self.slab[idx].symbol
            element_counts[symbol] = element_counts.get(symbol, 0) + 1

        # Get positions
        positions = {
            symbol: [
                self.slab.positions[i].tolist()
                for i in surface_indices
                if self.slab[i].symbol == symbol
            ]
            for symbol in element_counts
        }

        return {
            "n_surface_atoms": len(surface_indices),
            "z_cutoff": z_cutoff,
            "element_counts": element_counts,
            "positions": positions,
            "surface_indices": surface_indices,
        }


def filter_unique_configs(
    configs: List[Atoms],
    threshold: float = 0.5,
    n_slab_atoms: Optional[int] = None,
) -> List[Atoms]:
    """
    Filter configurations to remove duplicates based on adsorbate RMSD.

    Compares only the adsorbate atoms (not the slab) to determine uniqueness.
    This prevents the slab atoms from dominating the RMSD calculation.

    Parameters
    ----------
    configs : list
        List of Atoms objects (slab + adsorbate)
    threshold : float
        RMSD threshold in Angstroms for considering adsorbates as duplicates.
        Default 0.5 Å works well for most adsorbates.
    n_slab_atoms : int, optional
        Number of atoms in the original slab (before adsorbate was added).
        Adsorbate atoms are assumed to be at indices n_slab_atoms onwards.
        If not provided, it's inferred by finding atoms with identical
        positions across all configs.

    Returns
    -------
    list
        Filtered list of unique configurations

    Examples
    --------
    >>> # Filter configs generated from a 100-atom slab
    >>> unique = filter_unique_configs(configs, threshold=0.5, n_slab_atoms=100)
    """
    if not configs:
        return []

    if len(configs) == 1:
        return configs

    # Infer n_slab_atoms if not provided
    if n_slab_atoms is None:
        # Assume all configs have same slab, find minimum atom count
        # and check if positions match for those atoms
        min_atoms = min(len(c) for c in configs)
        # Use the first config's non-adsorbate atoms as reference
        # Heuristic: slab atoms are all atoms except the last few
        # Try to detect by checking which atoms have identical positions
        ref_positions = configs[0].positions
        for n in range(min_atoms, 0, -1):
            all_match = True
            for config in configs[1:]:
                if not np.allclose(ref_positions[:n], config.positions[:n], atol=0.01):
                    all_match = False
                    break
            if all_match:
                n_slab_atoms = n
                break
        else:
            # Fallback: assume last 4-10 atoms are adsorbate
            n_slab_atoms = min_atoms - 4
            print(f"Warning: Could not infer n_slab_atoms, assuming {n_slab_atoms}")

    def get_adsorbate_positions(config: Atoms) -> np.ndarray:
        """Extract adsorbate positions (atoms after slab)."""
        return config.positions[n_slab_atoms:]

    def adsorbate_rmsd(config1: Atoms, config2: Atoms) -> float:
        """Calculate RMSD of adsorbate atoms only."""
        pos1 = get_adsorbate_positions(config1)
        pos2 = get_adsorbate_positions(config2)

        if len(pos1) != len(pos2):
            return float('inf')

        diff = pos1 - pos2
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    unique = [configs[0]]

    for config in configs[1:]:
        is_unique = True
        for existing in unique:
            rmsd = adsorbate_rmsd(config, existing)
            if rmsd < threshold:
                is_unique = False
                break

        if is_unique:
            unique.append(config)

    n_adsorbate = len(configs[0]) - n_slab_atoms
    print(f"Filtered {len(configs)} configs to {len(unique)} unique "
          f"(compared {n_adsorbate} adsorbate atoms, threshold={threshold} Å)")
    return unique


def save_configs(
    configs: List[Atoms],
    output_dir: Union[str, Path],
    prefix: str = "config",
    format: str = "cif",
) -> List[Path]:
    """
    Save configurations to files.

    Parameters
    ----------
    configs : list
        List of Atoms objects
    output_dir : str or Path
        Output directory
    prefix : str
        Filename prefix
    format : str
        Output format

    Returns
    -------
    list
        List of output file paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = []
    for i, config in enumerate(configs):
        filename = output_dir / f"{prefix}_{i:03d}.{format}"
        write(str(filename), config)
        paths.append(filename)

    print(f"Saved {len(configs)} configurations to {output_dir}")
    return paths
