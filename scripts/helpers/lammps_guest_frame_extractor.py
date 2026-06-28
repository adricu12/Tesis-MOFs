"""
lammps_frame_extractor.py

Extract a selected frame from a LAMMPS trajectory containing:
    host atoms first, then guest atoms
where the host atom count is read automatically from a LAMMPS data.host file.

Main notebook function:
    data = extract_frame_data(...)

Then you can use:
    gen_inp_file(template, symbols=data["symbols"], coords=data["coords"],
                 lattice=data["lattice"], name="...", output="...")

Assumptions:
    - The trajectory dump has atom columns: id type x y z, or xu yu zu.
    - Atom IDs are sorted internally before slicing.
    - Host atoms are the atoms defined in data_host.
    - Any atoms beyond the host count in the trajectory are guest atoms.
    - The host data file is a replicated supercell, e.g. 2x2x2.
    - Guest elements can be inferred automatically from a LAMMPS data file
      containing the guest atom types in its Masses section, or from an explicit mapping.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from math import acos, degrees
from pathlib import Path
import re
from typing import Iterable

import numpy as np


# Used only for reconstructing/unwrapping small guest molecules.
COVALENT_RADII = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
    "Zr": 1.45,
}

MASS_TO_ELEMENT = [
    (1.008, "H"),
    (12.011, "C"),
    (14.007, "N"),
    (15.999, "O"),
    (18.998, "F"),
    (30.974, "P"),
    (32.06, "S"),
    (35.45, "Cl"),
    (79.904, "Br"),
    (91.224, "Zr"),
    (126.904, "I"),
]


@dataclass
class HostData:
    path: Path
    n_atoms: int
    type_to_element: dict[int, str]
    atom_ids: np.ndarray
    atom_types: np.ndarray
    atom_symbols: list[str]
    coords: np.ndarray | None
    origin: np.ndarray | None
    H_super: np.ndarray | None


@dataclass
class Frame:
    index: int
    timestep: int
    origin: np.ndarray
    H_super: np.ndarray
    ids: np.ndarray
    types: np.ndarray
    coords: np.ndarray
    columns: list[str]


def _clean_element_from_comment(comment: str) -> str | None:
    """Extract element from comments such as '# Zr8f4', '# C_R', '# O_3_f'."""
    if not comment:
        return None
    token = comment.strip().split()[0]
    match = re.match(r"([A-Z][a-z]?)", token)
    if match:
        return match.group(1)
    return None


def _infer_element_from_mass(mass: float, tolerance: float = 0.35) -> str | None:
    best = min(MASS_TO_ELEMENT, key=lambda x: abs(x[0] - mass))
    if abs(best[0] - mass) <= tolerance:
        return best[1]
    return None


def _parse_header_atom_count(lines: list[str]) -> int:
    for line in lines[:100]:
        m = re.match(r"\s*(\d+)\s+atoms\b", line)
        if m:
            return int(m.group(1))
    raise ValueError("Could not find '<N> atoms' in the LAMMPS data file header.")


def _parse_lammps_data_box(lines: list[str]) -> tuple[np.ndarray, np.ndarray] | tuple[None, None]:
    xlo = xhi = ylo = yhi = zlo = zhi = None
    xy = xz = yz = 0.0

    for line in lines[:120]:
        parts = line.split()
        if len(parts) >= 4 and parts[-2:] == ["xlo", "xhi"]:
            xlo, xhi = float(parts[0]), float(parts[1])
        elif len(parts) >= 4 and parts[-2:] == ["ylo", "yhi"]:
            ylo, yhi = float(parts[0]), float(parts[1])
        elif len(parts) >= 4 and parts[-2:] == ["zlo", "zhi"]:
            zlo, zhi = float(parts[0]), float(parts[1])
        elif len(parts) >= 6 and parts[-3:] == ["xy", "xz", "yz"]:
            xy, xz, yz = float(parts[0]), float(parts[1]), float(parts[2])

    if None in (xlo, xhi, ylo, yhi, zlo, zhi):
        return None, None

    origin = np.array([xlo, ylo, zlo], dtype=float)
    H = np.array(
        [
            [xhi - xlo, 0.0, 0.0],
            [xy, yhi - ylo, 0.0],
            [xz, yz, zhi - zlo],
        ],
        dtype=float,
    )
    return origin, H


def _find_section(lines: list[str], name: str) -> int | None:
    for i, line in enumerate(lines):
        if line.strip().split("#", 1)[0].strip() == name:
            return i
    return None


def _parse_masses(lines: list[str]) -> dict[int, str]:
    start = _find_section(lines, "Masses")
    if start is None:
        return {}

    type_to_element: dict[int, str] = {}
    i = start + 1

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        i += 1

        if not stripped:
            continue
        if re.match(r"^[A-Za-z]", stripped):
            break

        body, _, comment = line.partition("#")
        parts = body.split()
        if len(parts) < 2:
            continue
        try:
            atom_type = int(parts[0])
            mass = float(parts[1])
        except ValueError:
            continue

        element = _clean_element_from_comment(comment) or _infer_element_from_mass(mass)
        if element is not None:
            type_to_element[atom_type] = element

    return type_to_element


def _parse_atoms_section(lines: list[str], n_atoms: int, type_to_element: dict[int, str]):
    start = _find_section(lines, "Atoms")
    if start is None:
        return None, None, [], None

    ids: list[int] = []
    types: list[int] = []
    coords: list[list[float]] = []

    i = start + 1
    while i < len(lines) and len(ids) < n_atoms:
        line = lines[i].strip()
        i += 1
        if not line or line.startswith("#"):
            continue
        if re.match(r"^[A-Za-z]", line):
            break

        parts = line.split("#", 1)[0].split()
        if len(parts) < 6:
            continue

        # Handles common atom_style full:
        # id mol type q x y z [ix iy iz]
        # and atom_style charge:
        # id type q x y z
        try:
            atom_id = int(parts[0])
            if len(parts) >= 7:
                atom_type = int(parts[2])
                xyz = [float(parts[4]), float(parts[5]), float(parts[6])]
            else:
                atom_type = int(parts[1])
                xyz = [float(parts[3]), float(parts[4]), float(parts[5])]
        except ValueError:
            continue

        ids.append(atom_id)
        types.append(atom_type)
        coords.append(xyz)

    if not ids:
        return None, None, [], None

    ids_arr = np.array(ids, dtype=int)
    types_arr = np.array(types, dtype=int)
    coords_arr = np.array(coords, dtype=float)
    order = np.argsort(ids_arr)

    ids_arr = ids_arr[order]
    types_arr = types_arr[order]
    coords_arr = coords_arr[order]
    symbols = [type_to_element.get(int(t), "X") for t in types_arr]

    return ids_arr, types_arr, symbols, coords_arr


def read_host_data(data_file: str | Path) -> HostData:
    path = Path(data_file)
    lines = path.read_text().splitlines()

    n_atoms = _parse_header_atom_count(lines)
    type_to_element = _parse_masses(lines)
    origin, H_super = _parse_lammps_data_box(lines)
    atom_ids, atom_types, atom_symbols, coords = _parse_atoms_section(lines, n_atoms, type_to_element)

    if atom_ids is None:
        atom_ids = np.arange(1, n_atoms + 1, dtype=int)
        atom_types = np.array([], dtype=int)
        atom_symbols = []
        coords = None

    return HostData(
        path=path,
        n_atoms=n_atoms,
        type_to_element=type_to_element,
        atom_ids=atom_ids,
        atom_types=atom_types,
        atom_symbols=atom_symbols,
        coords=coords,
        origin=origin,
        H_super=H_super,
    )


def read_type_to_element_from_data(data_file: str | Path) -> dict[int, str]:
    """
    Read atom type -> element from the Masses section of a LAMMPS data file.

    Best case: Masses lines contain comments such as '# C_R', '# O_3', '# H'.
    Fallback: infer element from atomic mass, e.g. 12.011 -> C.

    This can be used with:
        - data_host: maps only host types
        - data_guest: maps only guest types
        - data_system / merged data file: maps host + guest types
    """
    lines = Path(data_file).read_text().splitlines()
    mapping = _parse_masses(lines)
    if not mapping:
        raise ValueError(f"Could not parse a Masses section from {data_file!s}.")
    return mapping


def merge_type_mappings(*mappings: dict[int, str] | None) -> dict[int, str]:
    """Merge atom-type mappings, preserving later values if duplicated."""
    out: dict[int, str] = {}
    for mapping in mappings:
        if mapping:
            out.update(mapping)
    return out


def lattice_from_lammps_bounds(bounds: list[list[float]]) -> tuple[np.ndarray, np.ndarray]:
    xlo_b, xhi_b, xy = bounds[0]
    ylo_b, yhi_b, xz = bounds[1]
    zlo_b, zhi_b, yz = bounds[2]

    xlo = xlo_b - min(0.0, xy, xz, xy + xz)
    xhi = xhi_b - max(0.0, xy, xz, xy + xz)
    ylo = ylo_b - min(0.0, yz)
    yhi = yhi_b - max(0.0, yz)
    zlo = zlo_b
    zhi = zhi_b

    origin = np.array([xlo, ylo, zlo], dtype=float)
    H = np.array(
        [
            [xhi - xlo, 0.0, 0.0],
            [xy, yhi - ylo, 0.0],
            [xz, yz, zhi - zlo],
        ],
        dtype=float,
    )
    return origin, H


def iter_lammps_frames(path: str | Path) -> Iterable[Frame]:
    path = Path(path)
    frame_index = 0

    with path.open("r") as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                raise ValueError(f"Expected ITEM: TIMESTEP, got: {line!r}")

            timestep = int(f.readline().strip())

            header = f.readline().strip()
            if not header.startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError(f"Expected ITEM: NUMBER OF ATOMS, got: {header!r}")
            n_atoms = int(f.readline().strip())

            box_header = f.readline().strip()
            if not box_header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Expected ITEM: BOX BOUNDS, got: {box_header!r}")
            bounds = [list(map(float, f.readline().split())) for _ in range(3)]
            origin, H_super = lattice_from_lammps_bounds(bounds)

            atom_header = f.readline().strip()
            if not atom_header.startswith("ITEM: ATOMS"):
                raise ValueError(f"Expected ITEM: ATOMS, got: {atom_header!r}")

            columns = atom_header.split()[2:]
            col = {name: i for i, name in enumerate(columns)}
            if "id" not in col or "type" not in col:
                raise ValueError("Dump must contain at least 'id' and 'type' columns.")

            # Prefer unwrapped coordinates if available; otherwise use wrapped x/y/z.
            if {"xu", "yu", "zu"}.issubset(col):
                xyz_cols = ["xu", "yu", "zu"]
            elif {"x", "y", "z"}.issubset(col):
                xyz_cols = ["x", "y", "z"]
            else:
                raise ValueError("Dump must contain either x y z or xu yu zu columns.")

            ids = np.empty(n_atoms, dtype=int)
            types = np.empty(n_atoms, dtype=int)
            coords = np.empty((n_atoms, 3), dtype=float)

            for i in range(n_atoms):
                parts = f.readline().split()
                ids[i] = int(parts[col["id"]])
                types[i] = int(parts[col["type"]])
                coords[i] = [float(parts[col[c]]) for c in xyz_cols]

            order = np.argsort(ids)
            yield Frame(
                index=frame_index,
                timestep=timestep,
                origin=origin,
                H_super=H_super,
                ids=ids[order],
                types=types[order],
                coords=coords[order],
                columns=columns,
            )
            frame_index += 1


def get_frame(path: str | Path, frame_index: int | None = None, timestep: int | None = None) -> Frame:
    if (frame_index is None) == (timestep is None):
        raise ValueError("Give exactly one of frame_index or timestep.")

    for frame in iter_lammps_frames(path):
        if frame_index is not None and frame.index == frame_index:
            return frame
        if timestep is not None and frame.timestep == timestep:
            return frame

    target = f"frame_index={frame_index}" if frame_index is not None else f"timestep={timestep}"
    raise ValueError(f"Could not find {target}.")


def cart_to_frac(cart: np.ndarray, origin: np.ndarray, H: np.ndarray) -> np.ndarray:
    return (cart - origin) @ np.linalg.inv(H)


def frac_to_cart(frac: np.ndarray, origin: np.ndarray, H: np.ndarray) -> np.ndarray:
    return frac @ H + origin


def minimum_image_delta_frac(delta_frac: np.ndarray) -> np.ndarray:
    return delta_frac - np.round(delta_frac)


def infer_bonds_from_distances(
    frac_wrapped: np.ndarray,
    elements: list[str],
    H: np.ndarray,
    scale: float = 1.25,
    tolerance: float = 0.25,
) -> list[tuple[int, int]]:
    bonds: list[tuple[int, int]] = []
    n = len(elements)
    for i in range(n):
        for j in range(i + 1, n):
            ei, ej = elements[i], elements[j]
            if ei == "H" and ej == "H":
                continue
            ri = COVALENT_RADII.get(ei, 0.75)
            rj = COVALENT_RADII.get(ej, 0.75)
            delta = minimum_image_delta_frac(frac_wrapped[j] - frac_wrapped[i])
            dist = np.linalg.norm(delta @ H)
            cutoff = scale * (ri + rj) + tolerance
            if dist <= cutoff:
                bonds.append((i, j))
    return bonds


def unwrap_guest_molecule(
    guest_cart: np.ndarray,
    guest_elements: list[str],
    origin: np.ndarray,
    H_super: np.ndarray,
) -> np.ndarray:
    frac = cart_to_frac(guest_cart, origin, H_super)
    frac_wrapped = frac % 1.0
    bonds = infer_bonds_from_distances(frac_wrapped, guest_elements, H_super)

    adjacency: dict[int, list[int]] = defaultdict(list)
    for i, j in bonds:
        adjacency[i].append(j)
        adjacency[j].append(i)

    n = len(guest_elements)
    unwrapped = np.full((n, 3), np.nan, dtype=float)
    visited = np.zeros(n, dtype=bool)
    unwrapped[0] = frac_wrapped[0]
    visited[0] = True
    queue = deque([0])

    while queue:
        i = queue.popleft()
        for j in adjacency[i]:
            if visited[j]:
                continue
            delta = minimum_image_delta_frac(frac_wrapped[j] - frac_wrapped[i])
            unwrapped[j] = unwrapped[i] + delta
            visited[j] = True
            queue.append(j)

    # Fallback for disconnected atoms/fragments.
    while not np.all(visited):
        missing = np.where(~visited)[0]
        present = np.where(visited)[0]
        best = None
        for j in missing:
            for i in present:
                delta = minimum_image_delta_frac(frac_wrapped[j] - frac_wrapped[i])
                dist = np.linalg.norm(delta @ H_super)
                if best is None or dist < best[0]:
                    best = (dist, i, j, delta)
        _, i, j, delta = best
        unwrapped[j] = unwrapped[i] + delta
        visited[j] = True

    return unwrapped


def lengths_angles_from_lattice(H: np.ndarray) -> tuple[float, float, float, float, float, float]:
    a_vec, b_vec, c_vec = H
    a = np.linalg.norm(a_vec)
    b = np.linalg.norm(b_vec)
    c = np.linalg.norm(c_vec)
    alpha = degrees(acos(np.clip(np.dot(b_vec, c_vec) / (b * c), -1.0, 1.0)))
    beta = degrees(acos(np.clip(np.dot(a_vec, c_vec) / (a * c), -1.0, 1.0)))
    gamma = degrees(acos(np.clip(np.dot(a_vec, b_vec) / (a * b), -1.0, 1.0)))
    return a, b, c, alpha, beta, gamma


def _detect_guest_tile(guest_frac_unwrapped: np.ndarray, nreplica: np.ndarray) -> np.ndarray:
    center = guest_frac_unwrapped.mean(axis=0)
    tile = np.floor(center * nreplica).astype(int)
    return np.mod(tile, nreplica)


def _select_host_tile(host_frac: np.ndarray, tile: np.ndarray, nreplica: np.ndarray, expected: int) -> np.ndarray:
    host_tiles = np.floor((host_frac % 1.0) * nreplica + 1e-9).astype(int)
    host_tiles = np.clip(host_tiles, 0, nreplica - 1)
    mask = np.all(host_tiles == tile, axis=1)
    if int(mask.sum()) != expected:
        raise ValueError(
            f"Selected {int(mask.sum())} host atoms for tile {tile.tolist()}, "
            f"but expected {expected}. Check nreplica={tuple(nreplica)}."
        )
    return mask


def _symbols_from_types(types: np.ndarray, type_to_element: dict[int, str], label: str) -> list[str]:
    symbols = [type_to_element.get(int(t), "X") for t in types]
    missing = sorted({int(t) for t, s in zip(types, symbols) if s == "X"})
    if missing:
        raise ValueError(
            f"Missing element mapping for {label} atom types: {missing}. "
            f"Pass a mapping such as guest_type_to_element={{10:'C', 11:'O', ...}}."
        )
    return symbols


def extract_frame_data(
    trajectory: str | Path,
    data_host: str | Path,
    *,
    frame_index: int | None = None,
    timestep: int | None = None,
    nreplica: tuple[int, int, int] = (2, 2, 2),
    guest_type_to_element: dict[int, str] | None = None,
    data_guest: str | Path | None = None,
    data_system: str | Path | None = None,
    include_host: bool = True,
    validate_system_atom_count: bool = True,
) -> dict:
    """
    Return one selected frame as a dictionary ready for gen_inp_file.

    The returned dictionary contains:
        symbols: list[str]
        coords: list[list[float]]      Cartesian coordinates in the selected unit-cell basis
        lattice: list[list[float]]     Unit-cell lattice vectors as rows
        frac_coords: np.ndarray        Fractional coordinates in the selected unit-cell basis
        host / guest nested dictionaries
        metadata: timestep, frame index, tile, atom counts
    """
    host_data = read_host_data(data_host)
    system_atom_count = None
    system_type_to_element = None
    if data_system is not None:
        system_data = read_host_data(data_system)
        system_atom_count = system_data.n_atoms
        system_type_to_element = system_data.type_to_element

    frame = get_frame(trajectory, frame_index=frame_index, timestep=timestep)

    nrep = np.array(nreplica, dtype=int)
    n_tiles = int(np.prod(nrep))
    host_atoms = host_data.n_atoms

    if frame.coords.shape[0] <= host_atoms:
        raise ValueError(
            f"Trajectory frame has {frame.coords.shape[0]} atoms, but host data has {host_atoms}. "
            "No extra guest atoms were found."
        )

    if validate_system_atom_count and system_atom_count is not None and frame.coords.shape[0] != system_atom_count:
        raise ValueError(
            f"Trajectory frame has {frame.coords.shape[0]} atoms, but data_system has {system_atom_count} atoms. "
            "Check that data.initial_before_min belongs to the same run as this trajectory."
        )

    if system_atom_count is not None and system_atom_count <= host_atoms:
        raise ValueError(
            f"data_system has {system_atom_count} atoms, but data_host has {host_atoms} atoms. "
            "The system file should contain host + guest atoms."
        )
    if host_atoms % n_tiles != 0:
        raise ValueError(
            f"Host atom count {host_atoms} is not divisible by prod(nreplica)={n_tiles}."
        )

    atoms_per_unit_cell = host_atoms // n_tiles
    guest_atoms = frame.coords.shape[0] - host_atoms

    host_types = frame.types[:host_atoms]
    host_cart = frame.coords[:host_atoms]
    guest_types = frame.types[host_atoms:]
    guest_cart = frame.coords[host_atoms:]

    host_symbols = _symbols_from_types(host_types, host_data.type_to_element, "host")

    # Build atom type -> element mapping for the guest.
    # Priority:
    #   1. explicit guest_type_to_element dictionary
    #   2. Masses section from a guest-only data file
    #   3. Masses section from a merged host+guest data file
    #   4. host mapping, only useful if type IDs are shared
    auto_guest_maps = []
    if data_guest is not None:
        auto_guest_maps.append(read_type_to_element_from_data(data_guest))
    if system_type_to_element is not None:
        auto_guest_maps.append(system_type_to_element)

    guest_type_to_element = merge_type_mappings(
        host_data.type_to_element,
        *auto_guest_maps,
        guest_type_to_element,
    )
    guest_symbols = _symbols_from_types(guest_types, guest_type_to_element, "guest")

    host_frac_super = cart_to_frac(host_cart, frame.origin, frame.H_super)
    guest_frac_super_unwrapped = unwrap_guest_molecule(
        guest_cart,
        guest_symbols,
        frame.origin,
        frame.H_super,
    )

    tile = _detect_guest_tile(guest_frac_super_unwrapped, nrep)
    host_mask = _select_host_tile(host_frac_super, tile, nrep, atoms_per_unit_cell)

    H_unit = frame.H_super / nrep.reshape(3, 1)

    host_frac_unit = host_frac_super[host_mask] * nrep - tile
    guest_frac_unit = guest_frac_super_unwrapped * nrep - tile

    selected_host_symbols = [s for s, keep in zip(host_symbols, host_mask) if keep]

    if include_host:
        symbols = selected_host_symbols + guest_symbols
        frac_coords = np.vstack([host_frac_unit % 1.0, guest_frac_unit])
    else:
        symbols = guest_symbols
        frac_coords = guest_frac_unit

    coords = frac_to_cart(frac_coords, np.zeros(3), H_unit)

    return {
        "symbols": symbols,
        "coords": coords.tolist(),
        "lattice": H_unit.tolist(),
        "frac_coords": frac_coords,
        "host": {
            "symbols": selected_host_symbols,
            "frac_coords": host_frac_unit % 1.0,
            "coords": frac_to_cart(host_frac_unit % 1.0, np.zeros(3), H_unit).tolist(),
            "atom_count": len(selected_host_symbols),
        },
        "guest": {
            "symbols": guest_symbols,
            "frac_coords": guest_frac_unit,
            "coords": frac_to_cart(guest_frac_unit, np.zeros(3), H_unit).tolist(),
            "atom_count": len(guest_symbols),
            "types": guest_types.tolist(),
            "type_to_element": {int(t): guest_type_to_element[int(t)] for t in sorted(set(map(int, guest_types)))},
        },
        "metadata": {
            "frame_index": frame.index,
            "timestep": frame.timestep,
            "tile": tile.tolist(),
            "host_atoms_in_data_file": host_atoms,
            "guest_atoms_in_trajectory": guest_atoms,
            "atoms_per_unit_cell": atoms_per_unit_cell,
            "nreplica": tuple(nreplica),
            "trajectory_atoms_total": int(frame.coords.shape[0]),
            "data_system_atoms_total": system_atom_count,
            "data_host": str(host_data.path),
            "data_system": str(data_system) if data_system is not None else None,
            "trajectory": str(trajectory),
        },
    }


def write_cif_from_data(data: dict, out_path: str | Path, name: str = "extracted_frame") -> None:
    out_path = Path(out_path)
    H = np.array(data["lattice"], dtype=float)
    frac = np.array(data["frac_coords"], dtype=float)
    symbols = data["symbols"]

    a, b, c, alpha, beta, gamma = lengths_angles_from_lattice(H)
    counters: dict[str, int] = defaultdict(int)

    with out_path.open("w") as f:
        f.write(f"data_{name}\n")
        f.write("_symmetry_space_group_name_H-M   'P1'\n")
        f.write("_symmetry_Int_Tables_number      1\n")
        f.write(f"_cell_length_a    {a:.8f}\n")
        f.write(f"_cell_length_b    {b:.8f}\n")
        f.write(f"_cell_length_c    {c:.8f}\n")
        f.write(f"_cell_angle_alpha {alpha:.8f}\n")
        f.write(f"_cell_angle_beta  {beta:.8f}\n")
        f.write(f"_cell_angle_gamma {gamma:.8f}\n")
        f.write("loop_\n")
        f.write("_symmetry_equiv_pos_as_xyz\n")
        f.write("'x, y, z'\n")
        f.write("loop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")

        host_n = data.get("host", {}).get("atom_count", 0)
        for i, (sym, xyz) in enumerate(zip(symbols, frac)):
            counters[sym] += 1
            label = f"{sym}{counters[sym]}"
            # Wrap host atoms; do not wrap guest atoms, to avoid splitting molecule.
            write_xyz = xyz % 1.0 if i < host_n else xyz
            f.write(f"{label:8s} {sym:3s} {write_xyz[0]: .8f} {write_xyz[1]: .8f} {write_xyz[2]: .8f}\n")


def list_frames(trajectory: str | Path, limit: int | None = 20) -> None:
    for frame in iter_lammps_frames(trajectory):
        print(f"frame_index={frame.index:6d}  timestep={frame.timestep}")
        if limit is not None and frame.index + 1 >= limit:
            break
