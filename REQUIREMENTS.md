# Software requirements

Reproducing the computational pipeline requires several independent software packages spanning quantum chemistry, molecular dynamics, charge assignment, and data analysis. These packages are not managed by a single package manager and must be installed separately. The sections below list each dependency with its version, license, source, and any configuration notes required for compatibility with this pipeline.

> **HPC note.** All calculations in this thesis were performed on the ZIH HPC cluster at TU Dresden. LAMMPS was used as a pre-installed module on ZIH rather than a local build. Researchers reproducing this work on other HPC systems should verify that their LAMMPS installation includes the modules listed in the LAMMPS section below.

---

## Python environments

Two separate Python installations are required. They must not be mixed into a single environment.

### Python 3.12 (main analysis environment)

Used by all Jupyter notebooks, `pipeline_utils.py`, and the `mof-guest-toolkit`.

Recommended installation via [Miniforge](https://github.com/conda-forge/miniforge) (conda-forge channel).

```bash
conda create -n thesis-mof python=3.12
conda activate thesis-mof
conda install -c conda-forge rdkit h5py numpy pandas scipy jupyter ipykernel
pip install mof-guest-toolkit        # or pip install -e . from the cloned repo
```

| Package | Version used | Purpose |
|---------|-------------|---------|
| Python | 3.12 | Main interpreter |
| RDKit | 2024.03 | Guest molecular descriptors |
| h5py | 3.11 | HDF5 file I/O for SRD configurations and SPE results |
| NumPy | 1.26 | Numerical operations |
| pandas | 2.2 | Data handling and CSV I/O |
| SciPy | 1.13 | Pearson and Spearman correlation (`scipy.stats.pearsonr`, `scipy.stats.spearmanr`) |
| Jupyter | 7.x | Notebook execution |
| mof-guest-toolkit | see [repo](https://github.com/adricu12/mof-guest-toolkit) | CIF/XYZ conversion, Zeo++ interface, PubChem fetching |

### Python 2.7 (amber2lammps only)

Used exclusively to run `amber2lammps.py`, which converts Amber topology files to LAMMPS data format and is not compatible with Python 3.

> **Warning.** Python 2.7 reached end-of-life in January 2020 and is no longer maintained. It is not available via conda-forge. Options for obtaining a Python 2.7 installation: (1) system package manager on older Linux distributions (`apt install python2.7`); (2) [pyenv](https://github.com/pyenv/pyenv) (`pyenv install 2.7.18`); (3) port `amber2lammps.py` to Python 3 manually -- the script is short and the changes required are minimal (print statements, integer division, `urllib` imports).

---

## Quantum chemistry and semiempirical calculations

### AMS -- Amsterdam Modeling Suite

Version: 2024 (or the version available at your institution)\
License: Commercial. Academic licenses available from [SCM](https://www.scm.com/support/downloads/).\
Download: [https://www.scm.com/amsterdam-modeling-suite/](https://www.scm.com/amsterdam-modeling-suite/)

Used for: GFN1-xTB geometry optimization and single-point energy calculations on host frameworks, guest molecules, and host-guest complexes (via the AMS DFTB engine with D3(BJ) dispersion correction).

No special build configuration required. The AMS DFTB engine with GFN1-xTB parameters is included in all standard AMS distributions. Ensure `$AMSBIN` is set in your environment.

### ORCA

Version: 6.1.0\
License: Free for academic use. Registration required.\
Download: [https://orcaforum.kofo.mpg.de](https://orcaforum.kofo.mpg.de)

Used for: HF/6-31G(d) single-point energy calculations on guest molecules to obtain electron densities for RESP charge fitting.

After installation, ensure `orca` and `orca_2mkl` are on your `$PATH`. The `orca_2mkl` utility is required to generate `.molden.input` files for Multiwfn.

---

## Charge assignment

### EQeq (Extended Charge Equilibration)

Version: v1.00 (C++ implementation)\
License: MIT\
Source: [https://github.com/danieleongari/EQeq](https://github.com/danieleongari/EQeq)

Used for: Assigning partial charges to host framework atoms for MD simulations.

Compile with GCC:

```bash
git clone https://github.com/danieleongari/EQeq
cd EQeq
g++ EQeq_v1_00.cpp -o eqeq
```

The compiled binary `eqeq` must be accessible from the working directory or on your `$PATH`. Input CIF files require specific formatting before being passed to this code; this is handled by `mof_toolkit.cif_tools` (`simplify_cif`, `inject_eqeq_charges`), called from `scripts/modeling-and-simulation.ipynb`.

### Multiwfn

Version: 3.8\
License: Free for academic use.\
Download: [http://sobereva.com/multiwfn](http://sobereva.com/multiwfn)\
**Platform: Windows only** for this pipeline.

Used for: Two-stage RESP charge fitting from HF/6-31G(d) electron densities.

> **Platform note.** Multiwfn is available for Linux, macOS, and Windows. However, the Linux build is prone to segmentation faults if the system stack and memory settings are not configured strictly according to Section 2.1.2 of the Multiwfn manual; the official documentation explicitly recommends using the Windows version when available (Tian Lu, *Multiwfn quick start*, 2026). The RESP charge fitting in this pipeline was therefore performed using the Multiwfn 3.8 Windows build driven by a `.cmd` batch script, provided in [`data/resp_multiwfn_commands.txt`](data/resp_multiwfn_commands.txt). Researchers wishing to use the Linux build must follow Section 2.1.2 of the Multiwfn manual precisely before running the batch script. Alternatively, the RESP fitting step can be replicated using `resp` from AmberTools or `pyRESP`, both of which run natively on Linux.

---

## Molecular dynamics

### LAMMPS

Version: 23Jun2022\
License: GPL-2.0\
Download: [https://www.lammps.org/download.html](https://www.lammps.org/download.html)

Used for: Classical MD simulations of host-guest complexes for configurational sampling.

Required LAMMPS packages (must be enabled at compile time if building from source):

| Package | Required for |
|---------|-------------|
| `KSPACE` | Long-range electrostatics via PPPM |
| `MOLECULE` | Full atom style, bonds, angles, dihedrals |
| `EXTRA-MOLECULE` | `fourier` angle and improper styles |
| `RIGID` | Not used directly but recommended |

If using a pre-built module on an HPC system, verify that these packages are included:

```bash
lmp -h 2>&1 | grep -E "KSPACE|MOLECULE|EXTRA"
```

The LAMMPS input files in this pipeline use `atom_style full`, `pair_style lj/cut/coul/long`, `kspace_style pppm`, `angle_style hybrid fourier cosine/periodic harmonic`, `dihedral_style harmonic`, and `improper_style fourier`. Any LAMMPS build that does not include the packages above will fail to parse these inputs.

### AmberTools 25

License: GPL-2.0 (free)\
Download: [https://ambermd.org/AmberTools.php](https://ambermd.org/AmberTools.php)

Used for: Antechamber (GAFF2 atom-typing), parmchk2 (missing parameter generation), and tLEaP (topology assembly).

Installation via conda is the recommended route:

```bash
conda install -c conda-forge ambertools=25
```

The following AmberTools executables are called directly by the pipeline:

| Executable | Purpose |
|------------|---------|
| `antechamber` | GAFF2 atom-typing and RESP charge assignment |
| `parmchk2` | Identify and supplement missing GAFF2 parameters |
| `tleap` | Assemble Amber topology (`.prmtop`) and coordinate (`.inpcrd`) files |

### amber2lammps

Version: as retrieved from source (no version tag)\
License: see source repository\
Source: [https://collaborating.tuhh.de/m-29/software/maxentrdf/-/blob/6e76eadc7d94ee6b4432bc9aeed463d498982031/tools/amber2lmp/amber2lammps.py](https://collaborating.tuhh.de/m-29/software/maxentrdf/-/blob/6e76eadc7d94ee6b4432bc9aeed463d498982031/tools/amber2lmp/amber2lammps.py)

Used for: Converting Amber `.top` and `.crd` files to LAMMPS `data.[guest]` format.

**Requires Python 2.7.** See the Python 2.7 note above. The script is invoked as:

```bash
python2 amber2lammps.py [guest].top [guest].crd data.[guest]
```

A copy of the script at the exact commit used in this work is provided in `scripts/helpers/amber2lammps.py` to ensure the specific version is reproducible independently of the upstream repository.

---

## Pore geometry analysis

### Zeo++

Version: 0.3\
License: BSD\
Download: [https://www.zeoplusplus.org/download.html](https://www.zeoplusplus.org/download.html)

Used for: Computing host framework pore descriptors (pore limiting diameter, largest cavity diameter, accessible surface area, accessible volume, probe-occupiable volume, pore size distribution) with a 1.55 Å probe radius.

No special build configuration required. The compiled `network` binary must be accessible from the working directory or on your `$PATH`.

---

## Host-guest configuration generation

### host_guest (modified)

Base version: as available at [https://github.com/bafgreat/host_guest](https://github.com/bafgreat/host_guest)\
License: see upstream repository\
Modified version: [https://github.com/adricu12/mof-guest-toolkit](https://github.com/adricu12/mof-guest-toolkit)

Used for: Stochastic rigid-body docking (SRD) to generate host-guest configurations.

Two modifications were made to the original source for this work: user-defined random seed support and HDF5 output format replacing the original JSON output. The modified version is included in the `mof-guest-toolkit` package. Install via:

```bash
pip install mof-guest-toolkit
# or for development:
git clone https://github.com/adricu12/mof-guest-toolkit
pip install -e .
```

### lammps_interface

Source: [https://github.com/peteboyd/lammps_interface](https://github.com/peteboyd/lammps_interface)\
License: see upstream repository

Used for: Generating LAMMPS data files for host frameworks from EQeq-charged CIF files (UFF4MOF force field assignment).

Install via:

```bash
pip install lammps_interface
# or from source:
git clone https://github.com/peteboyd/lammps_interface
pip install -e .
```

---

## Summary table

| Software | Version | Platform | License | Install route |
|----------|---------|----------|---------|---------------|
| Python | 3.12 | Linux/macOS/Windows | PSF | conda-forge |
| Python | 2.7 | Linux | PSF | pyenv or system |
| RDKit | 2024.03 | Linux/macOS/Windows | BSD | conda-forge |
| h5py | 3.11 | Linux/macOS/Windows | BSD | conda-forge |
| NumPy | 1.26 | Linux/macOS/Windows | BSD | conda-forge |
| pandas | 2.2 | Linux/macOS/Windows | BSD | conda-forge |
| SciPy | 1.13 | Linux/macOS/Windows | BSD | conda-forge |
| AMS | 2024 | Linux/macOS/Windows | Commercial | scm.com |
| ORCA | 6.1.0 | Linux/macOS/Windows | Free academic | orcaforum |
| EQeq | 1.00 | Linux/macOS/Windows | MIT | compile from source |
| Multiwfn | 3.8 | Windows recommended (Linux possible with configuration) | Free academic | sobereva.com |
| LAMMPS | 23Jun2022 | Linux | GPL-2.0 | HPC module or source |
| AmberTools | 25 | Linux/macOS | GPL-2.0 | conda-forge |
| amber2lammps | -- | Linux | see source | copy from source |
| Zeo++ | 0.3 | Linux/macOS/Windows | BSD | zeoplusplus.org |
| host_guest (modified) | -- | Linux/macOS/Windows | see source | pip (mof-guest-toolkit) |
| lammps_interface | -- | Linux/macOS/Windows | see source | pip |
| mof-guest-toolkit | see repo | Linux/macOS/Windows | see repo | pip |