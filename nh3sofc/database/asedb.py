"""ASE database integration for NH3-SOFC calculations.

Provides wrapper around ASE's database functionality with
specialized methods for catalysis calculations.
"""

from pathlib import Path
from typing import Optional, List, Dict, Any, Union, Iterator
import json
from ase import Atoms
from ase.db import connect
from ase.db.core import Database

from .naming import NamingConvention, generate_calc_id


class NH3SOFCDatabase:
    """
    ASE database wrapper for NH3-SOFC calculations.

    Provides convenient methods for storing and querying
    catalysis calculation results.

    Examples
    --------
    >>> db = NH3SOFCDatabase("calculations.db")
    >>> db.add_structure(atoms, material="LaVO3", calc_type="relax")
    >>> results = db.query(material="LaVO3", converged=True)
    """

    # Standard key-value pairs for NH3-SOFC calculations
    STANDARD_KEYS = [
        "material",
        "miller",
        "termination",
        "calc_type",
        "adsorbate",
        "converged",
        "vacancy_concentration",
        "nitrogen_fraction",
        "config_id",
    ]

    def __init__(self, db_path: Union[str, Path] = "nh3sofc.db"):
        """
        Initialize NH3SOFCDatabase.

        Parameters
        ----------
        db_path : str or Path
            Path to database file
        """
        self.db_path = Path(db_path)
        self.db = connect(str(self.db_path))

    def add_structure(
        self,
        atoms: Atoms,
        material: str,
        calc_type: str,
        miller: tuple = (0, 0, 1),
        termination: Optional[str] = None,
        adsorbate: Optional[str] = None,
        converged: bool = False,
        energy: Optional[float] = None,
        forces: Optional[Any] = None,
        vacancy_concentration: float = 0.0,
        nitrogen_fraction: float = 0.0,
        config_id: Optional[int] = None,
        data: Optional[Dict[str, Any]] = None,
        **kwargs,
    ) -> int:
        """
        Add a structure to the database.

        Parameters
        ----------
        atoms : Atoms
            Structure to add
        material : str
            Material name
        calc_type : str
            Calculation type (bulk, slab, relax, neb, freq)
        miller : tuple
            Miller indices
        termination : str, optional
            Surface termination
        adsorbate : str, optional
            Adsorbate species
        converged : bool
            Whether calculation converged
        energy : float, optional
            Total energy (eV)
        forces : array, optional
            Forces on atoms
        vacancy_concentration : float
            O vacancy concentration
        nitrogen_fraction : float
            N substitution fraction
        config_id : int, optional
            Configuration identifier
        data : dict, optional
            Additional data to store
        **kwargs : dict
            Additional key-value pairs

        Returns
        -------
        int
            Database row ID
        """
        # Generate calculation ID
        calc_id = generate_calc_id(material, miller, calc_type, adsorbate)

        # Prepare key-value pairs
        # Add "m" prefix to miller to avoid ASE DB interpreting as int
        miller_str = "m" + "".join(str(i) for i in miller)
        key_value_pairs = {
            "material": material,
            "miller": miller_str,
            "calc_type": calc_type,
            "converged": converged,
            "vacancy_concentration": vacancy_concentration,
            "nitrogen_fraction": nitrogen_fraction,
            "calc_id": calc_id,
        }

        if termination:
            key_value_pairs["termination"] = termination
        if adsorbate:
            key_value_pairs["adsorbate"] = adsorbate
        if config_id is not None:
            key_value_pairs["config_id"] = config_id

        # Add any additional key-value pairs (filter out reserved keys)
        reserved_keys = {"energy", "forces", "stress", "magmom", "magmoms", "charges"}
        for k, v in kwargs.items():
            if k not in reserved_keys:
                key_value_pairs[k] = v

        # Prepare data dictionary
        data_dict = data or {}
        if forces is not None:
            data_dict["forces"] = forces.tolist() if hasattr(forces, "tolist") else forces
        if energy is not None:
            data_dict["total_energy"] = energy

        # Create a copy of atoms and set single-point results if we have energy
        atoms_to_write = atoms.copy()
        if energy is not None:
            from ase.calculators.singlepoint import SinglePointCalculator
            calc = SinglePointCalculator(atoms_to_write, energy=energy)
            atoms_to_write.calc = calc

        # Add to database
        row_id = self.db.write(
            atoms_to_write,
            key_value_pairs=key_value_pairs,
            data=data_dict,
        )

        return row_id

    def add_relaxation_result(
        self,
        atoms: Atoms,
        material: str,
        miller: tuple = (0, 0, 1),
        termination: Optional[str] = None,
        adsorbate: Optional[str] = None,
        energy: Optional[float] = None,
        converged: bool = False,
        max_force: Optional[float] = None,
        n_steps: Optional[int] = None,
        **kwargs,
    ) -> int:
        """
        Add a relaxation result to the database.

        Parameters
        ----------
        atoms : Atoms
            Relaxed structure
        material : str
            Material name
        miller : tuple
            Miller indices
        termination : str, optional
            Surface termination
        adsorbate : str, optional
            Adsorbate species
        energy : float, optional
            Final energy
        converged : bool
            Whether optimization converged
        max_force : float, optional
            Maximum force on atoms
        n_steps : int, optional
            Number of optimization steps
        **kwargs : dict
            Additional parameters

        Returns
        -------
        int
            Database row ID
        """
        data = {}
        if max_force is not None:
            data["max_force"] = max_force
        if n_steps is not None:
            data["n_steps"] = n_steps

        return self.add_structure(
            atoms,
            material=material,
            calc_type="relax",
            miller=miller,
            termination=termination,
            adsorbate=adsorbate,
            converged=converged,
            energy=energy,
            data=data,
            **kwargs,
        )

    def add_neb_result(
        self,
        images: List[Atoms],
        material: str,
        miller: tuple = (0, 0, 1),
        initial_state: str = "",
        final_state: str = "",
        forward_barrier: Optional[float] = None,
        reverse_barrier: Optional[float] = None,
        reaction_energy: Optional[float] = None,
        converged: bool = False,
        **kwargs,
    ) -> List[int]:
        """
        Add NEB result to the database.

        Parameters
        ----------
        images : list
            List of NEB images
        material : str
            Material name
        miller : tuple
            Miller indices
        initial_state : str
            Initial state name
        final_state : str
            Final state name
        forward_barrier : float, optional
            Forward activation energy
        reverse_barrier : float, optional
            Reverse activation energy
        reaction_energy : float, optional
            Reaction energy
        converged : bool
            Whether NEB converged
        **kwargs : dict
            Additional parameters

        Returns
        -------
        list
            Database row IDs for each image
        """
        row_ids = []

        data = {
            "initial_state": initial_state,
            "final_state": final_state,
            "n_images": len(images),
        }

        if forward_barrier is not None:
            data["forward_barrier"] = forward_barrier
        if reverse_barrier is not None:
            data["reverse_barrier"] = reverse_barrier
        if reaction_energy is not None:
            data["reaction_energy"] = reaction_energy

        for i, image in enumerate(images):
            energy = image.get_potential_energy() if image.calc else None

            row_id = self.add_structure(
                image,
                material=material,
                calc_type="neb",
                miller=miller,
                converged=converged,
                energy=energy,
                config_id=i,
                data=data,
                neb_image=i,
                **kwargs,
            )
            row_ids.append(row_id)

        return row_ids

    def add_frequency_result(
        self,
        atoms: Atoms,
        material: str,
        frequencies: List[float],
        zpe: float,
        miller: tuple = (0, 0, 1),
        adsorbate: Optional[str] = None,
        electronic_energy: Optional[float] = None,
        **kwargs,
    ) -> int:
        """
        Add frequency calculation result.

        Parameters
        ----------
        atoms : Atoms
            Structure
        material : str
            Material name
        frequencies : list
            Vibrational frequencies (cm^-1)
        zpe : float
            Zero-point energy (eV)
        miller : tuple
            Miller indices
        adsorbate : str, optional
            Adsorbate species
        electronic_energy : float, optional
            Electronic energy
        **kwargs : dict
            Additional parameters

        Returns
        -------
        int
            Database row ID
        """
        data = {
            "frequencies": frequencies,
            "zpe": zpe,
            "n_frequencies": len(frequencies),
            "n_imaginary": sum(1 for f in frequencies if f < 0),
        }

        return self.add_structure(
            atoms,
            material=material,
            calc_type="freq",
            miller=miller,
            adsorbate=adsorbate,
            converged=True,
            energy=electronic_energy,
            data=data,
            **kwargs,
        )

    def query(
        self,
        material: Optional[str] = None,
        calc_type: Optional[str] = None,
        adsorbate: Optional[str] = None,
        converged: Optional[bool] = None,
        miller: Optional[str] = None,
        **kwargs,
    ) -> List[Any]:
        """
        Query the database.

        Parameters
        ----------
        material : str, optional
            Filter by material
        calc_type : str, optional
            Filter by calculation type
        adsorbate : str, optional
            Filter by adsorbate
        converged : bool, optional
            Filter by convergence status
        miller : str, optional
            Filter by Miller indices
        **kwargs : dict
            Additional filters

        Returns
        -------
        list
            Matching database rows
        """
        # Build selection string
        selections = []

        if material:
            selections.append(f"material={material}")
        if calc_type:
            selections.append(f"calc_type={calc_type}")
        if adsorbate:
            selections.append(f"adsorbate={adsorbate}")
        if converged is not None:
            selections.append(f"converged={converged}")
        if miller:
            selections.append(f"miller={miller}")

        for key, value in kwargs.items():
            if isinstance(value, str):
                selections.append(f"{key}={value}")
            else:
                selections.append(f"{key}={value}")

        selection = ",".join(selections) if selections else ""

        return list(self.db.select(selection))

    def get_structure(self, row_id: int) -> Atoms:
        """
        Get structure by row ID.

        Parameters
        ----------
        row_id : int
            Database row ID

        Returns
        -------
        Atoms
            Structure
        """
        return self.db.get_atoms(id=row_id)

    def get_lowest_energy(
        self,
        material: str,
        calc_type: str = "relax",
        adsorbate: Optional[str] = None,
        **kwargs,
    ) -> Optional[Any]:
        """
        Get lowest energy structure.

        Parameters
        ----------
        material : str
            Material name
        calc_type : str
            Calculation type
        adsorbate : str, optional
            Adsorbate filter
        **kwargs : dict
            Additional filters

        Returns
        -------
        row or None
            Database row with lowest energy
        """
        rows = self.query(
            material=material,
            calc_type=calc_type,
            adsorbate=adsorbate,
            converged=True,
            **kwargs,
        )

        if not rows:
            return None

        # Filter rows with energy
        rows_with_energy = [r for r in rows if r.energy is not None]

        if not rows_with_energy:
            return None

        return min(rows_with_energy, key=lambda r: r.energy)

    def get_decomposition_energies(
        self,
        material: str,
        miller: str = "001",
    ) -> Dict[str, float]:
        """
        Get decomposition pathway energies.

        Parameters
        ----------
        material : str
            Material name
        miller : str
            Miller indices

        Returns
        -------
        dict
            Step energies
        """
        steps = ["NH3", "NH2_H", "NH_2H", "N_3H"]
        energies = {}

        for step in steps:
            row = self.get_lowest_energy(
                material=material,
                calc_type="relax",
                adsorbate=step,
                miller=miller,
            )
            if row:
                energies[step] = row.energy

        return energies

    def export_to_json(
        self,
        filename: str,
        selection: Optional[str] = None,
    ) -> None:
        """
        Export database to JSON.

        Parameters
        ----------
        filename : str
            Output file path
        selection : str, optional
            Selection string
        """
        rows = list(self.db.select(selection or ""))

        data = []
        for row in rows:
            entry = {
                "id": row.id,
                "formula": row.formula,
                "energy": row.energy,
                "key_value_pairs": dict(row.key_value_pairs),
            }
            if row.data:
                entry["data"] = dict(row.data)
            data.append(entry)

        with open(filename, "w") as f:
            json.dump(data, f, indent=2)

        print(f"Exported {len(data)} entries to {filename}")

    def get_statistics(self) -> Dict[str, Any]:
        """
        Get database statistics.

        Returns
        -------
        dict
            Statistics summary
        """
        all_rows = list(self.db.select())

        stats = {
            "total_entries": len(all_rows),
            "materials": {},
            "calc_types": {},
            "adsorbates": {},
        }

        for row in all_rows:
            kvp = row.key_value_pairs

            mat = kvp.get("material", "unknown")
            stats["materials"][mat] = stats["materials"].get(mat, 0) + 1

            calc = kvp.get("calc_type", "unknown")
            stats["calc_types"][calc] = stats["calc_types"].get(calc, 0) + 1

            ads = kvp.get("adsorbate")
            if ads:
                stats["adsorbates"][ads] = stats["adsorbates"].get(ads, 0) + 1

        return stats

    def print_summary(self) -> None:
        """Print database summary."""
        stats = self.get_statistics()

        print("\n" + "=" * 50)
        print("NH3-SOFC Database Summary")
        print("=" * 50)

        print(f"\nTotal entries: {stats['total_entries']}")

        print("\nBy material:")
        for mat, count in sorted(stats["materials"].items()):
            print(f"  {mat}: {count}")

        print("\nBy calculation type:")
        for calc, count in sorted(stats["calc_types"].items()):
            print(f"  {calc}: {count}")

        if stats["adsorbates"]:
            print("\nBy adsorbate:")
            for ads, count in sorted(stats["adsorbates"].items()):
                print(f"  {ads}: {count}")

        print("=" * 50)


def create_database(db_path: str = "nh3sofc.db") -> NH3SOFCDatabase:
    """
    Create a new NH3-SOFC database.

    Parameters
    ----------
    db_path : str
        Database file path

    Returns
    -------
    NH3SOFCDatabase
        New database instance
    """
    return NH3SOFCDatabase(db_path)
