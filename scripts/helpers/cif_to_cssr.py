from pathlib import Path
import re


def _clean_cif_number(value: str) -> float:
    """
    Convert CIF numbers like '12.345(6)' into float 12.345.
    """
    value = value.strip().strip("'\"")
    value = re.sub(r"\(.+\)$", "", value)
    return float(value)


def _guess_element(label: str) -> str:
    """
    Guess element symbol from atom label.

    Examples
    --------
    C12 -> C
    Zr1 -> Zr
    O_carboxylate -> O
    """
    label = label.strip().strip("'\"")
    match = re.match(r"([A-Z][a-z]?)", label)
    if not match:
        raise ValueError(f"Could not guess element from atom label: {label}")
    return match.group(1)


def read_cif_cell_and_atoms(cif_file_path):
    """
    Read cell parameters and fractional atom coordinates from a CIF file.

    Returns
    -------
    cell_params : dict
        Keys: a, b, c, alpha, beta, gamma

    atoms : list[tuple]
        Each atom is (element, fract_x, fract_y, fract_z)
    """
    cif_file_path = Path(cif_file_path)

    cell_params = {}
    atoms = []

    lines = cif_file_path.read_text().splitlines()

    # Read cell parameters
    for raw_line in lines:
        line = raw_line.strip()

        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if len(parts) < 2:
            continue

        key = parts[0]
        value = parts[1]

        if key == "_cell_length_a":
            cell_params["a"] = _clean_cif_number(value)
        elif key == "_cell_length_b":
            cell_params["b"] = _clean_cif_number(value)
        elif key == "_cell_length_c":
            cell_params["c"] = _clean_cif_number(value)
        elif key == "_cell_angle_alpha":
            cell_params["alpha"] = _clean_cif_number(value)
        elif key == "_cell_angle_beta":
            cell_params["beta"] = _clean_cif_number(value)
        elif key == "_cell_angle_gamma":
            cell_params["gamma"] = _clean_cif_number(value)

    required_cell_keys = {"a", "b", "c", "alpha", "beta", "gamma"}
    missing = required_cell_keys - set(cell_params)
    if missing:
        raise ValueError(f"Missing cell parameters in {cif_file_path}: {sorted(missing)}")

    # Read atom-site loop
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.lower() != "loop_":
            i += 1
            continue

        headers = []
        i += 1

        while i < len(lines) and lines[i].strip().startswith("_"):
            headers.append(lines[i].strip())
            i += 1

        # Check whether this loop is the atom-site loop
        if not any(h.startswith("_atom_site_") for h in headers):
            continue

        try:
            x_idx = headers.index("_atom_site_fract_x")
            y_idx = headers.index("_atom_site_fract_y")
            z_idx = headers.index("_atom_site_fract_z")
        except ValueError:
            continue

        # Element can be stored as _atom_site_type_symbol or inferred from label
        element_idx = None
        label_idx = None

        if "_atom_site_type_symbol" in headers:
            element_idx = headers.index("_atom_site_type_symbol")

        if "_atom_site_label" in headers:
            label_idx = headers.index("_atom_site_label")

        if element_idx is None and label_idx is None:
            raise ValueError(
                f"Could not find _atom_site_type_symbol or _atom_site_label in {cif_file_path}"
            )

        # Data rows continue until another loop_, data_, or header starts
        while i < len(lines):
            atom_line = lines[i].strip()

            if (
                not atom_line
                or atom_line.startswith("#")
                or atom_line.startswith("_")
                or atom_line.lower().startswith("loop_")
                or atom_line.lower().startswith("data_")
            ):
                break

            parts = atom_line.split()

            if len(parts) < len(headers):
                i += 1
                continue

            if element_idx is not None:
                element = parts[element_idx].strip("'\"")
            else:
                element = _guess_element(parts[label_idx])

            x = _clean_cif_number(parts[x_idx])
            y = _clean_cif_number(parts[y_idx])
            z = _clean_cif_number(parts[z_idx])

            atoms.append((element, x, y, z))
            i += 1

        break

    if not atoms:
        raise ValueError(f"No fractional atom coordinates found in {cif_file_path}")

    return cell_params, atoms


def cif_to_cssr_direct(cif_file_path, cssr_file_path=None, verbose=True):
    """
    Convert a CIF file to CSSR format, preserving fractional coordinates.

    Parameters
    ----------
    cif_file_path : str or pathlib.Path
        Input CIF file.

    cssr_file_path : str or pathlib.Path, optional
        Output CSSR file. If None, uses the same name as the CIF file.

    verbose : bool
        Print summary information.

    Returns
    -------
    pathlib.Path
        Path to the created CSSR file.
    """
    cif_file_path = Path(cif_file_path)
    name = cif_file_path.stem

    if cssr_file_path is None:
        cssr_file_path = cif_file_path.with_suffix(".cssr")
    else:
        cssr_file_path = Path(cssr_file_path)

    cssr_file_path.parent.mkdir(parents=True, exist_ok=True)

    cell_params, atoms = read_cif_cell_and_atoms(cif_file_path)

    with cssr_file_path.open("w") as f:
        f.write(
            f"{cell_params['a']:10.4f}"
            f"{cell_params['b']:10.4f}"
            f"{cell_params['c']:10.4f}\n"
        )

        f.write(
            f" {cell_params['alpha']:4.0f}"
            f" {cell_params['beta']:4.0f}"
            f" {cell_params['gamma']:4.0f}"
            "  SPGR =  1 P 1         OPT = 1\n"
        )

        f.write(f"{len(atoms):4d}   0\n")
        f.write(f"0 {name:20s}: {name}\n")

        for atom_index, (element, x, y, z) in enumerate(atoms, start=1):
            elem = element.strip().capitalize().ljust(2)[:2]
            f.write(
                f"{atom_index:4d} {elem} "
                f"{x:10.6f} {y:10.6f} {z:10.6f} "
                "0 0 0 0 0 0 0 0 0.00\n"
            )

    if verbose:
        print(f"Converted: {cif_file_path} -> {cssr_file_path}")
        print(f"Structure name: {name}")
        print(f"Atoms: {len(atoms)}")
        print(
            "Cell: "
            f"a={cell_params['a']:.4f}, "
            f"b={cell_params['b']:.4f}, "
            f"c={cell_params['c']:.4f}"
        )
        print(
            "Angles: "
            f"alpha={cell_params['alpha']:.1f}, "
            f"beta={cell_params['beta']:.1f}, "
            f"gamma={cell_params['gamma']:.1f}"
        )

    return cssr_file_path


def batch_convert_cif_to_cssr(cif_directory=".", output_directory=None, pattern="*.cif", verbose=True):
    """
    Convert all CIF files in a directory to CSSR files.

    Parameters
    ----------
    cif_directory : str or pathlib.Path
        Folder containing CIF files.

    output_directory : str or pathlib.Path, optional
        Folder for CSSR outputs. If None, writes next to the CIF files.

    pattern : str
        File pattern, usually '*.cif'.

    verbose : bool
        Print conversion progress.

    Returns
    -------
    list[pathlib.Path]
        Paths to successfully created CSSR files.
    """
    cif_directory = Path(cif_directory)

    if output_directory is None:
        output_directory = cif_directory
    else:
        output_directory = Path(output_directory)

    output_directory.mkdir(parents=True, exist_ok=True)

    cif_files = sorted(cif_directory.glob(pattern))
    created_files = []

    if verbose:
        print(f"Found {len(cif_files)} CIF files in {cif_directory}")

    for cif_file in cif_files:
        cssr_file = output_directory / f"{cif_file.stem}.cssr"

        try:
            created = cif_to_cssr_direct(
                cif_file_path=cif_file,
                cssr_file_path=cssr_file,
                verbose=verbose,
            )
            created_files.append(created)

        except Exception as error:
            print(f"Failed: {cif_file.name}: {error}")

    return created_files


def verify_conversion(cif_file_path, cssr_file_path):
    """
    Compare the first atom in the CIF and CSSR files.
    """
    cif_file_path = Path(cif_file_path)
    cssr_file_path = Path(cssr_file_path)

    _, atoms = read_cif_cell_and_atoms(cif_file_path)

    first_cif_atom = atoms[0]

    cssr_lines = cssr_file_path.read_text().splitlines()
    if len(cssr_lines) < 5:
        raise ValueError(f"CSSR file has too few lines: {cssr_file_path}")

    first_cssr_parts = cssr_lines[4].split()

    print(
        "CIF first atom: "
        f"{first_cif_atom[0]} "
        f"({first_cif_atom[1]:.6f}, {first_cif_atom[2]:.6f}, {first_cif_atom[3]:.6f})"
    )

    print(
        "CSSR first atom: "
        f"{first_cssr_parts[1]} "
        f"({first_cssr_parts[2]}, {first_cssr_parts[3]}, {first_cssr_parts[4]})"
    )