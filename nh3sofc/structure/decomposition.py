"""NH3 decomposition intermediate configuration generation.

Generates structures for the NH3 decomposition pathway:
NH3* → NH2* + H* → NH* + 2H* → N* + 3H*

Uses co-adsorption approach where H atoms remain on the surface
during decomposition (physically realistic for catalytic conditions).
"""

from pathlib import Path
from typing import Optional, List, Union, Tuple, Dict, Any
import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.io import write

from ..core.base import BaseStructure, get_surface_atoms, calculate_rmsd
from .surface import SlabStructure


class DecompositionBuilder:
    """
    Generator for NH3 decomposition intermediate configurations.

    Creates multiple configurations for each decomposition step
    with different H atom placements.

    Examples
    --------
    >>> # Start from optimized NH3 on surface
    >>> decomp = DecompositionBuilder(nh3_on_slab)
    >>>
    >>> # Generate NH2 + H configurations
    >>> nh2_h_configs = decomp.create_NH2_H_configs(n_configs=10)
    >>>
    >>> # Generate NH + 2H configurations
    >>> nh_2h_configs = decomp.create_NH_2H_configs(n_configs=10)
    >>>
    >>> # Generate N + 3H configurations
    >>> n_3h_configs = decomp.create_N_3H_configs(n_configs=5)
    """

    def __init__(
        self,
        nh3_on_slab: Atoms,
        n_slab_atoms: Optional[int] = None,
    ):
        """
        Initialize DecompositionBuilder.

        Parameters
        ----------
        nh3_on_slab : Atoms
            Structure with NH3 adsorbed on surface
        n_slab_atoms : int, optional
            Number of atoms in the bare slab. If None, auto-detect
            by assuming last 4 atoms are NH3.
        """
        self.atoms = nh3_on_slab.copy()

        # Identify slab vs adsorbate atoms
        if n_slab_atoms is None:
            # Assume NH3 was added last (4 atoms: N + 3H)
            n_slab_atoms = len(self.atoms) - 4

        self.n_slab_atoms = n_slab_atoms
        self.slab = self.atoms[:n_slab_atoms]

        # Identify NH3 atoms
        adsorbate_atoms = self.atoms[n_slab_atoms:]
        self.nh3_indices = list(range(n_slab_atoms, len(self.atoms)))

        # Find N and H in NH3
        self.n_index = None
        self.h_indices = []
        for i, idx in enumerate(self.nh3_indices):
            if self.atoms[idx].symbol == "N":
                self.n_index = idx
            elif self.atoms[idx].symbol == "H":
                self.h_indices.append(idx)

        if self.n_index is None:
            raise ValueError("Could not find N atom in adsorbate")
        if len(self.h_indices) != 3:
            raise ValueError(f"Expected 3 H atoms, found {len(self.h_indices)}")

    def _get_surface_sites(
        self,
        z_threshold: float = 0.2,
        atom_types: Optional[List[str]] = None,
    ) -> List[np.ndarray]:
        """Get positions of surface sites for H placement."""
        surface_indices, _ = get_surface_atoms(self.slab, z_threshold)

        if atom_types is not None:
            surface_indices = [
                i for i in surface_indices
                if self.slab[i].symbol in atom_types
            ]

        return [self.slab.positions[i] for i in surface_indices]

    def _get_random_surface_position(
        self,
        z_offset: float = 1.0,
        exclude_region: Optional[Tuple[np.ndarray, float]] = None,
    ) -> np.ndarray:
        """
        Get a random position on the surface.

        Parameters
        ----------
        z_offset : float
            Height above surface
        exclude_region : tuple, optional
            (center, radius) to avoid placing H too close to NH2/NH/N

        Returns
        -------
        ndarray
            (x, y, z) position
        """
        cell = self.slab.get_cell()
        a, b = cell[0, 0], cell[1, 1]
        z_top = self.slab.positions[:, 2].max()

        max_attempts = 100
        for _ in range(max_attempts):
            x = np.random.uniform(0, a)
            y = np.random.uniform(0, b)
            z = z_top + z_offset

            if exclude_region is not None:
                center, radius = exclude_region
                dist = np.sqrt((x - center[0])**2 + (y - center[1])**2)
                if dist < radius:
                    continue

            return np.array([x, y, z])

        # If we couldn't find a good spot, return anyway
        return np.array([np.random.uniform(0, a),
                        np.random.uniform(0, b),
                        z_top + z_offset])

    def create_NH2_H_configs(
        self,
        n_configs: int = 10,
        h_height: float = 1.0,
        min_h_distance: float = 2.0,
        site_types: Optional[List[str]] = None,
        random_seed: Optional[int] = None,
    ) -> List[Atoms]:
        """
        Create NH2* + H* configurations (Step 1: NH3* → NH2* + H*).

        Removes one H from NH3 and places it on different surface sites.

        Parameters
        ----------
        n_configs : int
            Number of configurations to generate
        h_height : float
            Height of dissociated H above surface
        min_h_distance : float
            Minimum distance between H and NH2
        site_types : list, optional
            Atom types to use as H adsorption sites
        random_seed : int, optional
            Random seed

        Returns
        -------
        list
            List of Atoms objects
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        configs = []
        n_pos = self.atoms.positions[self.n_index]

        # Get surface sites for H placement
        surface_sites = self._get_surface_sites(atom_types=site_types)
        z_top = self.slab.positions[:, 2].max()

        for config_num in range(n_configs):
            # Choose which H to remove (rotate through the 3 H atoms)
            h_to_remove = self.h_indices[config_num % 3]

            # Create new structure
            new_atoms = self.atoms.copy()

            # Get position for the dissociated H
            if config_num < len(surface_sites):
                # Use a surface site
                site_pos = surface_sites[config_num % len(surface_sites)]
                h_new_pos = np.array([site_pos[0], site_pos[1], z_top + h_height])
            else:
                # Random position avoiding NH2
                h_new_pos = self._get_random_surface_position(
                    z_offset=h_height,
                    exclude_region=(n_pos[:2], min_h_distance),
                )

            # Move the H to its new position
            new_atoms.positions[h_to_remove] = h_new_pos

            configs.append(new_atoms)

        return configs

    def create_NH_2H_configs(
        self,
        n_configs: int = 10,
        h_height: float = 1.0,
        min_h_h_distance: float = 1.5,
        random_seed: Optional[int] = None,
    ) -> List[Atoms]:
        """
        Create NH* + 2H* configurations (Step 2: NH2* → NH* + 2H*).

        Parameters
        ----------
        n_configs : int
            Number of configurations
        h_height : float
            Height of H atoms above surface
        min_h_h_distance : float
            Minimum distance between H atoms
        random_seed : int, optional
            Random seed

        Returns
        -------
        list
            List of Atoms objects
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        configs = []
        n_pos = self.atoms.positions[self.n_index]
        z_top = self.slab.positions[:, 2].max()
        cell = self.slab.get_cell()
        a, b = cell[0, 0], cell[1, 1]

        for config_num in range(n_configs):
            new_atoms = self.atoms.copy()

            # Keep one H on N (the first one), move the other two
            h_to_move = self.h_indices[1:]  # Move 2nd and 3rd H

            h_positions = []
            for i, h_idx in enumerate(h_to_move):
                max_attempts = 50
                for _ in range(max_attempts):
                    # Random position
                    x = np.random.uniform(0, a)
                    y = np.random.uniform(0, b)
                    z = z_top + h_height

                    pos = np.array([x, y, z])

                    # Check distance from N
                    dist_to_n = np.linalg.norm(pos[:2] - n_pos[:2])
                    if dist_to_n < 2.0:
                        continue

                    # Check distance from other H
                    too_close = False
                    for prev_pos in h_positions:
                        if np.linalg.norm(pos - prev_pos) < min_h_h_distance:
                            too_close = True
                            break

                    if not too_close:
                        h_positions.append(pos)
                        break
                else:
                    # Couldn't find good position, use random
                    h_positions.append(np.array([
                        np.random.uniform(0, a),
                        np.random.uniform(0, b),
                        z_top + h_height
                    ]))

            # Update positions
            for i, h_idx in enumerate(h_to_move):
                new_atoms.positions[h_idx] = h_positions[i]

            configs.append(new_atoms)

        return configs

    def create_N_3H_configs(
        self,
        n_configs: int = 5,
        h_height: float = 1.0,
        arrangement: str = "mixed",
        random_seed: Optional[int] = None,
    ) -> List[Atoms]:
        """
        Create N* + 3H* configurations (Step 3: NH* → N* + 3H*).

        Parameters
        ----------
        n_configs : int
            Number of configurations
        h_height : float
            Height of H atoms
        arrangement : str
            H arrangement: "spread", "clustered", or "mixed"
        random_seed : int, optional
            Random seed

        Returns
        -------
        list
            List of Atoms objects
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        configs = []
        n_pos = self.atoms.positions[self.n_index]
        z_top = self.slab.positions[:, 2].max()
        cell = self.slab.get_cell()
        a, b = cell[0, 0], cell[1, 1]

        for config_num in range(n_configs):
            new_atoms = self.atoms.copy()

            # Determine arrangement for this config
            if arrangement == "spread":
                use_spread = True
            elif arrangement == "clustered":
                use_spread = False
            else:  # mixed
                use_spread = config_num % 2 == 0

            if use_spread:
                # Spread H atoms across the surface
                angles = np.linspace(0, 2*np.pi, 3, endpoint=False)
                angles += np.random.uniform(0, 2*np.pi)  # Random rotation
                radius = min(a, b) / 4

                for i, h_idx in enumerate(self.h_indices):
                    x = n_pos[0] + radius * np.cos(angles[i])
                    y = n_pos[1] + radius * np.sin(angles[i])
                    # Wrap to cell
                    x = x % a
                    y = y % b
                    new_atoms.positions[h_idx] = [x, y, z_top + h_height]
            else:
                # Clustered H atoms (may form H2)
                cluster_center = self._get_random_surface_position(
                    z_offset=h_height,
                    exclude_region=(n_pos[:2], 2.0),
                )
                for i, h_idx in enumerate(self.h_indices):
                    offset = np.random.uniform(-1.0, 1.0, size=2)
                    x = cluster_center[0] + offset[0]
                    y = cluster_center[1] + offset[1]
                    x = x % a
                    y = y % b
                    new_atoms.positions[h_idx] = [x, y, z_top + h_height]

            configs.append(new_atoms)

        return configs

    def create_H2_formation_configs(
        self,
        n_configs: int = 5,
        h_h_distances: Optional[List[float]] = None,
        random_seed: Optional[int] = None,
    ) -> List[Atoms]:
        """
        Create 2H* → H2* configurations for H2 formation study.

        Parameters
        ----------
        n_configs : int
            Number of configurations
        h_h_distances : list, optional
            H-H distances to sample (default: [0.74, 1.0, 1.5, 2.0, 2.5])
        random_seed : int, optional
            Random seed

        Returns
        -------
        list
            List of Atoms objects with 2H on surface at various distances
        """
        if random_seed is not None:
            np.random.seed(random_seed)

        if h_h_distances is None:
            h_h_distances = [0.74, 1.0, 1.5, 2.0, 2.5]  # 0.74 is H2 bond length

        configs = []
        z_top = self.slab.positions[:, 2].max()
        cell = self.slab.get_cell()
        a, b = cell[0, 0], cell[1, 1]

        for d in h_h_distances:
            for _ in range(max(1, n_configs // len(h_h_distances))):
                # Create slab with just 2H
                new_atoms = self.slab.copy()

                # Random center position
                cx = np.random.uniform(0, a)
                cy = np.random.uniform(0, b)

                # Random orientation
                angle = np.random.uniform(0, 2*np.pi)

                # H positions
                h1_pos = [cx + d/2 * np.cos(angle),
                         cy + d/2 * np.sin(angle),
                         z_top + 1.0]
                h2_pos = [cx - d/2 * np.cos(angle),
                         cy - d/2 * np.sin(angle),
                         z_top + 1.0]

                # Add H atoms
                h1 = Atoms("H", positions=[h1_pos])
                h2 = Atoms("H", positions=[h2_pos])
                new_atoms = new_atoms + h1 + h2

                configs.append(new_atoms)

        return configs

    def filter_unique_configs(
        self,
        configs: List[Atoms],
        threshold: float = 0.5,
    ) -> List[Atoms]:
        """
        Filter configurations to remove duplicates based on RMSD.

        Parameters
        ----------
        configs : list
            List of Atoms objects
        threshold : float
            RMSD threshold for duplicate detection

        Returns
        -------
        list
            Filtered list of unique configurations
        """
        if not configs:
            return []

        unique = [configs[0]]

        for config in configs[1:]:
            is_unique = True
            for existing in unique:
                if len(config) != len(existing):
                    continue
                try:
                    rmsd = calculate_rmsd(config, existing)
                    if rmsd < threshold:
                        is_unique = False
                        break
                except Exception:
                    continue

            if is_unique:
                unique.append(config)

        return unique


def generate_decomposition_pathway(
    nh3_on_slab: Atoms,
    n_configs_per_step: int = 5,
    random_seed: Optional[int] = None,
) -> Dict[str, List[Atoms]]:
    """
    Generate complete NH3 decomposition pathway configurations.

    Parameters
    ----------
    nh3_on_slab : Atoms
        NH3 adsorbed on surface
    n_configs_per_step : int
        Number of configurations per step
    random_seed : int, optional
        Random seed

    Returns
    -------
    dict
        Dictionary with keys 'NH3', 'NH2_H', 'NH_2H', 'N_3H'
        and lists of Atoms as values
    """
    builder = DecompositionBuilder(nh3_on_slab)

    pathway = {
        "NH3": [nh3_on_slab.copy()],
        "NH2_H": builder.create_NH2_H_configs(n_configs_per_step, random_seed=random_seed),
        "NH_2H": builder.create_NH_2H_configs(n_configs_per_step, random_seed=random_seed),
        "N_3H": builder.create_N_3H_configs(n_configs_per_step, random_seed=random_seed),
    }

    return pathway


def save_decomposition_configs(
    pathway: Dict[str, List[Atoms]],
    output_dir: Union[str, Path],
    format: str = "cif",
) -> Dict[str, List[Path]]:
    """
    Save all decomposition pathway configurations to files.

    Parameters
    ----------
    pathway : dict
        Dictionary from generate_decomposition_pathway
    output_dir : str or Path
        Output directory
    format : str
        Output file format

    Returns
    -------
    dict
        Dictionary mapping step names to list of file paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    file_paths = {}

    for step_name, configs in pathway.items():
        step_dir = output_dir / step_name
        step_dir.mkdir(exist_ok=True)

        paths = []
        for i, config in enumerate(configs):
            filename = step_dir / f"{step_name}_{i:03d}.{format}"
            write(str(filename), config)
            paths.append(filename)

        file_paths[step_name] = paths

    return file_paths
