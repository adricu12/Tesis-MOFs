#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
#####################################################################################
# Edit: 13 June 2026, by Adriana Ugarte                                            #
# HDF5 alternative to the per-complex JSON output produced by                     #
# host_guest.io.filetyper.append_json_atom.                                       #
#                                                                                   #
# A JSON file repeats the (often large) host coordinates once per docked          #
# complex, which is unmanageable at the scale of ~20k host-guest complexes.       #
# Here, for every base_name (complex_id) the host atoms, lattice and atom         #
# labels are stored once, and only the guest coordinates - which differ for       #
# every docking attempt and every loaded guest molecule - vary.                   #
#####################################################################################
import numpy as np
import h5py
from ase import Atoms


def append_hdf5_atom(base_name, complex_molecules, host_system, monomer, number_of_monomers, filename):
    '''
    Append the result of a docker.Dock run to an HDF5 file.

    The data for `base_name` is stored as a group with the following layout:

        /<base_name>/lattice_vectors          (3, 3)      float64
        /<base_name>/pbc                      (3,)        bool
        /<base_name>/complex_ids              (n_complexes,) string
        /<base_name>/host/positions           (n_host_atoms, 3)            float64
        /<base_name>/host/labels              (n_host_atoms,)              string
        /<base_name>/host/atomic_numbers      (n_host_atoms,)              int64
        /<base_name>/guest/positions          (n_complexes, n_loads, n_guest_atoms, 3) float64
        /<base_name>/guest/labels             (n_guest_atoms,)             string
        /<base_name>/guest/atomic_numbers     (n_guest_atoms,)             int64

    Host atoms, lattice and labels are identical for every docked complex of a
    given base_name (the host system is placed without rotation/translation by
    Dock), so they are written once. Only the guest positions vary across the
    `n_complexes` docking attempts and the `n_loads` (number_of_monomers) guest
    copies placed in each attempt.

    parameters
    ----------
    base_name: str
        Identifier for this host-guest pair, e.g. 'HOST_GUEST'. Stored as the
        top-level HDF5 group (complex_id).
    complex_molecules: dict
        Mapping of complex name (e.g. 'complex_0') to ase.Atoms, as returned
        by docker.Dock.
    host_system: ase.Atoms
        The host structure passed to docker.Dock.
    monomer: ase.Atoms
        The guest molecule passed to docker.Dock.
    number_of_monomers: int
        Number of guest copies placed per complex (n_loads).
    filename: str
        Path to the HDF5 file. Created if it does not exist; existing data
        for `base_name` is replaced.

    Note
    ----
    This currently assumes a single host (number_of_host == 1), which places
    the host atoms contiguously at the start of every complex. Complexes
    whose atom count does not match n_host_atoms + number_of_monomers *
    n_guest_atoms (e.g. when Dock failed to place all guest copies) are
    skipped with a warning.
    '''
    n_host_atoms = len(host_system)
    n_guest_atoms = len(monomer)
    expected_n_atoms = n_host_atoms + number_of_monomers * n_guest_atoms

    complex_ids = []
    guest_positions = []
    for complex_name, atoms in complex_molecules.items():
        if len(atoms) != expected_n_atoms:
            print(f"Warning: skipping {base_name}/{complex_name}: expected "
                  f"{expected_n_atoms} atoms (host + {number_of_monomers} guest "
                  f"loads), got {len(atoms)}")
            continue
        positions = atoms.get_positions()
        guest_positions.append(
            positions[n_host_atoms:].reshape(number_of_monomers, n_guest_atoms, 3))
        complex_ids.append(complex_name)

    if guest_positions:
        guest_positions = np.stack(guest_positions, axis=0)
    else:
        guest_positions = np.zeros((0, number_of_monomers, n_guest_atoms, 3))

    string_dtype = h5py.string_dtype(encoding='utf-8')

    with h5py.File(filename, 'a') as f_obj:
        if base_name in f_obj:
            del f_obj[base_name]
        grp = f_obj.create_group(base_name)

        grp.attrs['n_host_atoms'] = n_host_atoms
        grp.attrs['n_guest_atoms'] = n_guest_atoms
        grp.attrs['number_of_monomers'] = number_of_monomers
        grp.attrs['n_complexes'] = len(complex_ids)

        grp.create_dataset('lattice_vectors', data=np.array(host_system.get_cell()))
        grp.create_dataset('pbc', data=np.array(host_system.get_pbc()))
        grp.create_dataset('complex_ids', data=np.array(complex_ids, dtype=object),
                            dtype=string_dtype)

        host_grp = grp.create_group('host')
        host_grp.create_dataset('positions', data=host_system.get_positions())
        host_grp.create_dataset('labels', data=np.array(host_system.get_chemical_symbols(), dtype=object),
                                 dtype=string_dtype)
        host_grp.create_dataset('atomic_numbers', data=host_system.get_atomic_numbers())

        guest_grp = grp.create_group('guest')
        guest_grp.create_dataset('positions', data=guest_positions, compression='gzip')
        guest_grp.create_dataset('labels', data=np.array(monomer.get_chemical_symbols(), dtype=object),
                                  dtype=string_dtype)
        guest_grp.create_dataset('atomic_numbers', data=monomer.get_atomic_numbers())

    return


def read_hdf5_atom(filename, base_name):
    '''
    Read back the data written by `append_hdf5_atom` for a given base_name.

    Returns a dict with keys 'lattice_vectors', 'pbc', 'complex_ids', 'host'
    and 'guest', mirroring the on-disk layout (host/guest each contain
    'positions', 'labels' and 'atomic_numbers').
    '''
    with h5py.File(filename, 'r') as f_obj:
        grp = f_obj[base_name]
        data = {
            'lattice_vectors': grp['lattice_vectors'][()],
            'pbc': grp['pbc'][()],
            'complex_ids': [c for c in grp['complex_ids'].asstr()[()]],
            'host': {
                'positions': grp['host']['positions'][()],
                'labels': list(grp['host']['labels'].asstr()[()]),
                'atomic_numbers': grp['host']['atomic_numbers'][()],
            },
            'guest': {
                'positions': grp['guest']['positions'][()],
                'labels': list(grp['guest']['labels'].asstr()[()]),
                'atomic_numbers': grp['guest']['atomic_numbers'][()],
            },
        }
    return data


def hdf5_to_ase_complex(filename, base_name, complex_index):
    '''
    Reconstruct the full ase.Atoms object (host + all guest loads) for one
    docked complex, equivalent to complex_molecules[complex_id] in the JSON
    output.

    parameters
    ----------
    filename: str
        Path to the HDF5 file.
    base_name: str
        complex_id group to read from.
    complex_index: int
        Index into the n_complexes dimension of guest/positions (and
        complex_ids).
    '''
    data = read_hdf5_atom(filename, base_name)

    host_positions = data['host']['positions']
    host_labels = data['host']['labels']
    guest_positions = data['guest']['positions'][complex_index].reshape(-1, 3)
    guest_labels = data['guest']['labels'] * data['guest']['positions'].shape[1]

    positions = np.vstack([host_positions, guest_positions])
    labels = host_labels + guest_labels

    pbc = data['pbc']
    if pbc.any():
        atoms = Atoms(symbols=labels, positions=positions,
                       cell=data['lattice_vectors'], pbc=pbc)
    else:
        atoms = Atoms(symbols=labels, positions=positions)
    return atoms


def list_base_names(filename):
    '''
    Return the list of base_name (complex_id) groups stored in an HDF5 file.
    '''
    with h5py.File(filename, 'r') as f_obj:
        return list(f_obj.keys())
