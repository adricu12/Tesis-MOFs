"""
pipeline_utils.py
=================
Thesis-specific orchestration helpers for the MOF framework-topology
selectivity study.

This module holds the glue that is *specific to this thesis* — the binding
energy definition, the SPE-pool filtering and candidate-selection rules, and
the selectivity metrics (topology ranking, descriptor correlations, CBD-relative
binding energies). Generic, reusable utilities (CIF/CSSR/Zeo++/LAMMPS/AMS file
handling) live in the installable ``mof-guest-toolkit`` package and are imported
from there, not re-implemented here.

The functions map onto the methodology sections:

* §6  Minimum binding energy calculation  →  :func:`binding_energy`,
  :func:`steric_clash_filter`, :func:`top_n_per_method`, :func:`pool_statistics`,
  :func:`reported_minimum_be`
* §7  Selectivity analysis  →  :func:`topology_ranking`, :func:`rank1_counts`,
  :func:`descriptor_correlations`, :func:`cbd_relative`

See ``../energetic-data-analysis/energetic-data-analysis.md`` for the equations
and rationale behind each step.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mof_toolkit import get_energy

# CODATA Hartree -> kcal/mol. AMS reports total energies in hartree; binding
# energies in this thesis are reported in kcal/mol.
HARTREE_TO_KCAL = 627.5094740631

# Topology / framework labels used consistently across the analysis.
TOPOLOGIES = ["shp", "csq", "scu", "sqc"]
FRAMEWORKS = {"shp": "PCN-223", "csq": "MOF-545", "scu": "NU-902", "sqc": "PCN-225"}

# Configurations whose raw SPE-derived BE exceeds this (kcal/mol) are treated as
# near-contact artifacts of the semiempirical Hamiltonian, not physical poses
# (see energetic-data-analysis §6.1).
STERIC_CLASH_THRESHOLD = 10.0


# ---------------------------------------------------------------------------
# Loading the committed result tables
# ---------------------------------------------------------------------------

def load_results(correlation_dir, descriptors_dir):
    """
    Load and clean the committed binding-energy and guest-descriptor tables.

    ``minimum_binding_energies.csv`` stores guests inconsistently
    (pharmaceuticals by full name, cannabinoids by abbreviation, with a couple of
    latin-1 mojibake artifacts where ``β``/``α`` were mangled). This reconciles
    those against ``guests_descriptors.csv`` and returns everything keyed by the
    canonical guest **abbreviation**.

    Parameters
    ----------
    correlation_dir : str or pathlib.Path
        Folder holding ``minimum_binding_energies.csv`` (i.e. ``070_CORRELATION``).
    descriptors_dir : str or pathlib.Path
        Folder holding ``guests_descriptors.csv`` (i.e. ``040_DESCRIPTORS``).

    Returns
    -------
    be : pandas.DataFrame
        Cleaned long table with columns ``host, guest, minBE, abbr``.
    be_table : pandas.DataFrame
        Minimum BE pivoted to guests (abbreviation) × hosts.
    guest_type : pandas.Series
        Abbreviation → ``"Pharmaceutical"`` / ``"Cannabinoid"``.
    gdesc : pandas.DataFrame
        Guest descriptor table (includes ``Abbreviation``, ``Guest_Type``).
    """
    correlation_dir = Path(correlation_dir)
    descriptors_dir = Path(descriptors_dir)

    be = pd.read_csv(correlation_dir / "minimum_binding_energies.csv", encoding="latin-1")
    gdesc = pd.read_csv(descriptors_dir / "guests_descriptors.csv")

    # Repair known latin-1 mojibake in the two estradiol names.
    be["guest"] = (
        be["guest"]
        .str.replace("17-?-Estradiol", "17-beta-Estradiol", regex=False)
        .str.replace("17-\xe0-Ethinylestradiol", "17-alpha-Ethinylestradiol", regex=False)
    )

    def _norm(s):  # case/separator-insensitive key for matching
        return str(s).lower().replace("-", "").replace(" ", "").replace("_", "")

    key_to_abbr = {}
    for row in gdesc.itertuples():
        key_to_abbr[_norm(row.Name)] = row.Abbreviation
        key_to_abbr[_norm(row.Abbreviation)] = row.Abbreviation

    be["abbr"] = be["guest"].map(lambda g: key_to_abbr.get(_norm(g)))
    unmatched = sorted(be.loc[be["abbr"].isna(), "guest"].unique())
    if unmatched:
        print("WARNING: unmatched guests in minimum_binding_energies.csv:", unmatched)

    be_table = (
        be.dropna(subset=["abbr"])
        .pivot_table(index="abbr", columns="host", values="minBE", aggfunc="min")
    )
    be_table.index.name = "guest"
    guest_type = gdesc.set_index("Abbreviation")["Guest_Type"]
    return be, be_table, guest_type, gdesc


# ---------------------------------------------------------------------------
# §6 — Binding energy
# ---------------------------------------------------------------------------

def total_energy_hartree(ams_out_path) -> float:
    """
    Return the GFN1-xTB total energy (hartree) from an AMS ``.out`` file.

    Thin wrapper over :func:`mof_toolkit.get_energy` so the notebooks never
    have to scan the raw output text. Works for both ``SinglePoint`` and
    ``GeometryOptimization`` tasks.
    """
    return get_energy(str(ams_out_path))["energy_hartree"]


def binding_energy(e_complex, e_host, e_guest, to_kcal: bool = True):
    """
    Binding energy ``BE = E_complex - (E_host + E_guest)``.

    Parameters
    ----------
    e_complex, e_host, e_guest : float or array-like
        Total energies (hartree) of the host-guest complex and the isolated
        host and guest, all at the same level of theory.
    to_kcal : bool, optional
        If True (default), convert the result from hartree to kcal/mol.

    Returns
    -------
    float or numpy.ndarray
        Binding energy. More negative = stronger binding.
    """
    e_complex = np.asarray(e_complex, dtype=float)
    be = e_complex - (np.asarray(e_host, dtype=float) + np.asarray(e_guest, dtype=float))
    if to_kcal:
        be = be * HARTREE_TO_KCAL
    return be if be.ndim else float(be)


def noise_level(be_values, threshold: float = STERIC_CLASH_THRESHOLD) -> float:
    """
    Percentage of configurations excluded by the steric-clash threshold (%NL).

    ``%NL = (configurations with BE > threshold) / N * 100`` — a diagnostic of
    how many sampled poses were near-contact artifacts (energetic-data-analysis
    §6.4).
    """
    be = np.asarray(be_values, dtype=float)
    if be.size == 0:
        return float("nan")
    return float((be > threshold).sum() / be.size * 100.0)


def steric_clash_filter(df: pd.DataFrame, be_col: str = "BE_kcal", threshold: float = STERIC_CLASH_THRESHOLD) -> pd.DataFrame:
    """
    Drop configurations whose raw BE exceeds the steric-clash threshold.

    Returns a copy of *df* keeping only rows with ``df[be_col] <= threshold``
    (energetic-data-analysis §6.1).
    """
    return df[df[be_col] <= threshold].copy()


def top_n_per_method(
    df: pd.DataFrame,
    n: int = 20,
    be_col: str = "BE_kcal",
    method_col: str = "method",
) -> pd.DataFrame:
    """
    Select the *n* lowest-BE configurations from each sampling method.

    SRD and MD candidates are selected independently (energetic-data-analysis
    §6.1) because their raw SPE values are not directly comparable before
    relaxation. Returns the concatenated top-*n* rows per method, sorted by BE.
    """
    return (
        df.sort_values(be_col)
        .groupby(method_col, group_keys=False)
        .head(n)
        .reset_index(drop=True)
    )


def pool_statistics(be_values, threshold: float = STERIC_CLASH_THRESHOLD) -> dict:
    """
    Descriptive statistics over a configuration pool's binding energies.

    Returns ``{n, mean, std, min, max, range, noise_level_pct}`` (population
    standard deviation, matching energetic-data-analysis §6.4). Statistics are
    computed over the steric-clash-filtered pool; ``noise_level_pct`` is
    reported over the *unfiltered* pool.
    """
    raw = np.asarray(be_values, dtype=float)
    kept = raw[raw <= threshold]
    if kept.size == 0:
        return {"n": 0, "mean": float("nan"), "std": float("nan"),
                "min": float("nan"), "max": float("nan"), "range": float("nan"),
                "noise_level_pct": noise_level(raw, threshold)}
    return {
        "n": int(kept.size),
        "mean": float(kept.mean()),
        "std": float(kept.std(ddof=0)),
        "min": float(kept.min()),
        "max": float(kept.max()),
        "range": float(kept.max() - kept.min()),
        "noise_level_pct": noise_level(raw, threshold),
    }


def reported_minimum_be(be_values) -> float:
    """Reported binding energy for a host-guest pair: the minimum over all relaxed candidates (§6.3)."""
    be = np.asarray(be_values, dtype=float)
    return float(be.min()) if be.size else float("nan")


# ---------------------------------------------------------------------------
# §7 — Selectivity analysis
# ---------------------------------------------------------------------------

def topology_ranking(be_table: pd.DataFrame, ascending: bool = True) -> pd.DataFrame:
    """
    Rank the topologies per guest by binding energy (rank 1 = most favourable).

    Parameters
    ----------
    be_table : pandas.DataFrame
        Reported minimum BEs with guests as rows and topologies as columns.
    ascending : bool, optional
        If True (default), the most negative BE (strongest binding) gets rank 1.

    Returns
    -------
    pandas.DataFrame
        Integer ranks, same shape as *be_table* (energetic-data-analysis §7.2).
    """
    return be_table.rank(axis=1, ascending=ascending, method="min").astype("Int64")


def rank1_counts(be_table: pd.DataFrame) -> pd.Series:
    """
    Count, per topology, how many guests rank it first (strongest binding).

    The aggregate figure of merit for multi-contaminant sorbent selection
    (energetic-data-analysis §7.2): how broad a range of guests each framework
    binds most strongly.
    """
    ranks = topology_ranking(be_table)
    return (ranks == 1).sum(axis=0).astype(int)


def descriptor_correlations(
    be_series: pd.Series,
    descriptor_df: pd.DataFrame,
    methods=("pearson", "spearman"),
) -> pd.DataFrame:
    """
    Correlate binding energy against a set of structural descriptors.

    Computes Pearson (linear) and/or Spearman (rank, monotonic) correlation
    coefficients and their two-sided p-values between *be_series* and each
    column of *descriptor_df*, aligned on a shared index (energetic-data-analysis
    §7.1). Spearman is included because the BE-descriptor relationship is not
    assumed to be linear.

    Returns
    -------
    pandas.DataFrame
        Indexed by descriptor name, with ``<method>_r`` and ``<method>_p``
        columns for each requested method.
    """
    from scipy.stats import pearsonr, spearmanr

    funcs = {"pearson": pearsonr, "spearman": spearmanr}
    rows = {}
    for descriptor in descriptor_df.columns:
        joined = pd.concat([be_series, descriptor_df[descriptor]], axis=1).dropna()
        if len(joined) < 3:
            continue
        x = joined.iloc[:, 0].to_numpy()
        y = joined.iloc[:, 1].to_numpy()
        record = {}
        for method in methods:
            coeff, pval = funcs[method](x, y)
            record[f"{method}_r"] = coeff
            record[f"{method}_p"] = pval
        rows[descriptor] = record
    return pd.DataFrame.from_dict(rows, orient="index")


def within_host_correlations(
    be_table: pd.DataFrame,
    descriptor_df: pd.DataFrame,
    guests,
    topologies,
    min_n: int = 4,
) -> pd.DataFrame:
    """
    Per-host Spearman correlation of every descriptor with BE, within a guest subset.

    This is the engine behind the **class-separated** correlation analysis (§7):
    restricting to one chemical class and correlating *within each host* removes
    the two-class size confound that dominates a pooled correlation. A descriptor
    is a robust within-group driver if its correlation is **sign-consistent**
    across the four hosts and significant in several of them.

    Parameters
    ----------
    be_table : pandas.DataFrame
        Binding energies, guests (index) × topologies (columns).
    descriptor_df : pandas.DataFrame
        Guest descriptors indexed by the same guest labels as *be_table*.
    guests : iterable
        Guest labels to include (e.g. one chemical class).
    topologies : iterable
        Column labels of *be_table* to iterate over.
    min_n : int, optional
        Minimum number of guests with both BE and a varying descriptor required to
        correlate within a host (default 4).

    Returns
    -------
    pandas.DataFrame
        Indexed by descriptor, sorted by ``mean_abs_rho`` descending, with one
        ``rho_<topology>`` column per host plus ``mean_abs_rho``, ``n_sig``
        (hosts with p < 0.05) and ``sign_consistent``.
    """
    from scipy.stats import spearmanr

    guests = [g for g in guests if g in be_table.index]
    rows = {}
    for descriptor in descriptor_df.columns:
        rhos, pvals, ok = {}, {}, True
        for t in topologies:
            col = be_table.loc[guests, t].dropna()
            x = descriptor_df[descriptor].reindex(col.index)
            if x.notna().sum() < min_n or x.var() == 0:
                ok = False
                break
            r, p = spearmanr(x.to_numpy(), col.to_numpy())
            rhos[t], pvals[t] = r, p
        if not ok:
            continue
        vals = np.array(list(rhos.values()))
        rows[descriptor] = {
            **{f"rho_{t}": rhos[t] for t in topologies},
            "mean_abs_rho": float(np.abs(vals).mean()),
            "n_sig": int(sum(p < 0.05 for p in pvals.values())),
            "sign_consistent": len(set(np.sign(vals))) == 1,
        }
    out = pd.DataFrame.from_dict(rows, orient="index")
    return out.sort_values("mean_abs_rho", ascending=False) if len(out) else out


def cbd_relative(be_table: pd.DataFrame, reference: str = "CBD") -> pd.DataFrame:
    """
    CBD-referenced relative binding energies ``ΔBE_rel = BE_guest - BE_CBD``.

    For each topology (column), subtract the reference guest's BE from every
    other guest's BE (energetic-data-analysis §7.3). A negative value means the
    guest binds more strongly than CBD in that topology.

    Parameters
    ----------
    be_table : pandas.DataFrame
        Cannabinoid BEs with guests as rows and topologies as columns.
    reference : str, optional
        Reference guest (row label), default ``"CBD"``.

    Returns
    -------
    pandas.DataFrame
        ΔBE_rel for every non-reference guest across all topologies.
    """
    if reference not in be_table.index:
        raise KeyError(f"Reference guest '{reference}' not found in the BE table index.")
    delta = be_table.subtract(be_table.loc[reference], axis=1)
    return delta.drop(index=reference)
