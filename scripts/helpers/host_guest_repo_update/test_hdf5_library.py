import os
import time

import numpy as np
import pytest
from ase import Atoms

from host_guest.io import coords_library, filetyper, hdf5_library
from host_guest.energy import docker
from host_guest.setter import generate_complexes


TEST_DATA = os.path.join(os.path.dirname(__file__), 'test_data')


@pytest.fixture(scope='module')
def host_system():
    return coords_library.load_data_as_ase(os.path.join(TEST_DATA, 'EDUSIF.cif'))


@pytest.fixture(scope='module')
def monomer():
    return coords_library.load_data_as_ase(os.path.join(TEST_DATA, 'biphenyl.xyz'))


def test_append_and_read_single_load(tmp_path, host_system, monomer):
    _, complex_molecules = docker.Dock(
        host_system, monomer, number_of_host=1, number_of_monomers=1,
        number_of_complexes=3, energy=False, seed=42)

    filename = str(tmp_path / 'complexes.h5')
    hdf5_library.append_hdf5_atom(
        'EDUSIF_biphenyl', complex_molecules, host_system, monomer, 1, filename)

    data = hdf5_library.read_hdf5_atom(filename, 'EDUSIF_biphenyl')

    n_host = len(host_system)
    n_guest = len(monomer)

    assert data['host']['positions'].shape == (n_host, 3)
    np.testing.assert_allclose(data['host']['positions'], host_system.get_positions())
    assert data['host']['labels'] == host_system.get_chemical_symbols()
    np.testing.assert_array_equal(data['host']['atomic_numbers'], host_system.get_atomic_numbers())

    np.testing.assert_allclose(data['lattice_vectors'], np.array(host_system.get_cell()))
    np.testing.assert_array_equal(data['pbc'], host_system.get_pbc())

    assert data['guest']['labels'] == monomer.get_chemical_symbols()
    np.testing.assert_array_equal(data['guest']['atomic_numbers'], monomer.get_atomic_numbers())
    assert data['guest']['positions'].shape == (len(complex_molecules), 1, n_guest, 3)
    assert sorted(data['complex_ids']) == sorted(complex_molecules.keys())

    for i, complex_id in enumerate(data['complex_ids']):
        original = complex_molecules[complex_id]
        reconstructed = hdf5_library.hdf5_to_ase_complex(filename, 'EDUSIF_biphenyl', i)
        np.testing.assert_allclose(reconstructed.get_positions(), original.get_positions())
        assert reconstructed.get_chemical_symbols() == original.get_chemical_symbols()


def test_multiple_guest_loads(tmp_path, host_system, monomer):
    n_loads = 3
    _, complex_molecules = docker.Dock(
        host_system, monomer, number_of_host=1, number_of_monomers=n_loads,
        number_of_complexes=2, energy=False, seed=7)

    filename = str(tmp_path / 'complexes.h5')
    hdf5_library.append_hdf5_atom(
        'EDUSIF_biphenyl_x3', complex_molecules, host_system, monomer, n_loads, filename)

    data = hdf5_library.read_hdf5_atom(filename, 'EDUSIF_biphenyl_x3')
    n_guest = len(monomer)
    assert data['guest']['positions'].shape == (len(complex_molecules), n_loads, n_guest, 3)

    for i, complex_id in enumerate(data['complex_ids']):
        original = complex_molecules[complex_id]
        reconstructed = hdf5_library.hdf5_to_ase_complex(filename, 'EDUSIF_biphenyl_x3', i)
        np.testing.assert_allclose(reconstructed.get_positions(), original.get_positions())
        assert reconstructed.get_chemical_symbols() == original.get_chemical_symbols()
        # host block of the reconstructed complex must match host_system exactly
        np.testing.assert_allclose(
            reconstructed.get_positions()[:len(host_system)], host_system.get_positions())


def test_multiple_base_names_and_overwrite(tmp_path, host_system, monomer):
    filename = str(tmp_path / 'complexes.h5')

    _, complexes_a = docker.Dock(host_system, monomer, 1, 1, 2, energy=False, seed=1)
    _, complexes_b = docker.Dock(host_system, monomer, 1, 2, 2, energy=False, seed=2)

    hdf5_library.append_hdf5_atom('A_B', complexes_a, host_system, monomer, 1, filename)
    hdf5_library.append_hdf5_atom('C_D', complexes_b, host_system, monomer, 2, filename)

    assert sorted(hdf5_library.list_base_names(filename)) == ['A_B', 'C_D']

    # re-running for an existing base_name overwrites without leaving stale data
    _, complexes_a2 = docker.Dock(host_system, monomer, 1, 1, 5, energy=False, seed=3)
    hdf5_library.append_hdf5_atom('A_B', complexes_a2, host_system, monomer, 1, filename)

    assert sorted(hdf5_library.list_base_names(filename)) == ['A_B', 'C_D']
    data = hdf5_library.read_hdf5_atom(filename, 'A_B')
    assert data['guest']['positions'].shape[0] == len(complexes_a2)


def test_skips_malformed_complex(tmp_path, host_system, monomer):
    filename = str(tmp_path / 'complexes.h5')
    _, complex_molecules = docker.Dock(host_system, monomer, 1, 1, 3, energy=False, seed=4)

    # simulate a failed docking attempt that left an incomplete complex
    bad_key = next(iter(complex_molecules))
    complex_molecules[bad_key] = complex_molecules[bad_key][:-1]

    hdf5_library.append_hdf5_atom('BAD', complex_molecules, host_system, monomer, 1, filename)

    data = hdf5_library.read_hdf5_atom(filename, 'BAD')
    assert len(data['complex_ids']) == len(complex_molecules) - 1
    assert bad_key not in data['complex_ids']


def test_matches_json_encoding(tmp_path, host_system, monomer):
    _, complex_molecules = docker.Dock(host_system, monomer, 1, 1, 1, energy=False, seed=99)

    filename = str(tmp_path / 'complexes.h5')
    hdf5_library.append_hdf5_atom('X_Y', complex_molecules, host_system, monomer, 1, filename)

    complex_id = next(iter(complex_molecules))
    original = complex_molecules[complex_id]
    encoded = filetyper.AtomsEncoder().default(original)

    data = hdf5_library.read_hdf5_atom(filename, 'X_Y')
    idx = data['complex_ids'].index(complex_id)
    reconstructed = hdf5_library.hdf5_to_ase_complex(filename, 'X_Y', idx)

    np.testing.assert_allclose(reconstructed.get_positions(), np.array(encoded['positions']))
    assert reconstructed.get_chemical_symbols() == encoded['labels']
    np.testing.assert_allclose(np.array(reconstructed.get_cell()), np.array(encoded['lattice_vectors']))


def test_complexes_from_file_writes_hdf5_not_json(tmp_path):
    results_folder = str(tmp_path / 'results')
    host_file = os.path.join(TEST_DATA, 'EDUSIF.cif')
    monomer_file = os.path.join(TEST_DATA, 'biphenyl.xyz')

    generate_complexes.complexes_from_file(
        host_file, monomer_file, number_of_host=1, number_of_monomers=1,
        number_of_complexes=2, energy=False, results_folder=results_folder,
        seed=123, hdf5=True)

    h5_path = os.path.join(results_folder, 'complexes.h5')
    json_path = os.path.join(results_folder, 'complexes.json')
    assert os.path.exists(h5_path)
    assert not os.path.exists(json_path)

    host_base_name = os.path.basename(host_file).split('.')[0]
    monomer_base_name = os.path.basename(monomer_file).split('.')[0]
    base_name = host_base_name + '_' + monomer_base_name
    assert base_name in hdf5_library.list_base_names(h5_path)


def test_large_scale_synthetic(tmp_path):
    '''
    Stress-test append/read with many complexes and multiple guest loads,
    representative of a 20k-complex production run, without paying for the
    cost of the random docking search itself.
    '''
    rng = np.random.default_rng(0)
    n_host = 400
    n_guest = 22
    n_loads = 4
    n_complexes = 500

    host_positions = rng.random((n_host, 3)) * 25.832
    host_system = Atoms('C' * n_host, positions=host_positions,
                         cell=np.eye(3) * 25.832, pbc=True)
    monomer = Atoms('C' * n_guest, positions=rng.random((n_guest, 3)) * 5)

    complex_molecules = {}
    expected_guest_positions = np.zeros((n_complexes, n_loads, n_guest, 3))
    for i in range(n_complexes):
        guest_block = rng.random((n_loads * n_guest, 3)) * 25.832
        expected_guest_positions[i] = guest_block.reshape(n_loads, n_guest, 3)
        positions = np.vstack([host_positions, guest_block])
        symbols = host_system.get_chemical_symbols() + monomer.get_chemical_symbols() * n_loads
        complex_molecules[f'complex_{i}'] = Atoms(
            symbols=symbols, positions=positions, cell=host_system.get_cell(), pbc=True)

    filename = str(tmp_path / 'large.h5')

    t0 = time.time()
    hdf5_library.append_hdf5_atom(
        'LARGE_HOST_LARGE_GUEST', complex_molecules, host_system, monomer, n_loads, filename)
    write_time = time.time() - t0
    assert write_time < 30, f'writing {n_complexes} complexes took too long: {write_time:.1f}s'

    t0 = time.time()
    data = hdf5_library.read_hdf5_atom(filename, 'LARGE_HOST_LARGE_GUEST')
    read_time = time.time() - t0
    assert read_time < 30, f'reading {n_complexes} complexes took too long: {read_time:.1f}s'

    assert data['guest']['positions'].shape == (n_complexes, n_loads, n_guest, 3)
    np.testing.assert_allclose(data['guest']['positions'], expected_guest_positions)
    np.testing.assert_allclose(data['host']['positions'], host_positions)
    assert data['host']['positions'].nbytes < 1e6

    # host data is stored once, not once per complex
    file_size = os.path.getsize(filename)
    naive_json_like_size = n_complexes * (n_host + n_loads * n_guest) * 3 * 8
    assert file_size < naive_json_like_size
