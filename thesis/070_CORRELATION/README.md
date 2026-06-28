# 070_CORRELATION — final results

The primary committed outputs of the thesis. Every file here is regenerable: the
**derived** tables are rebuilt from `minimum_binding_energies.csv` + the descriptors
by [`scripts/selectivity-analysis.ipynb`](../../scripts/selectivity-analysis.ipynb)
(methodology §7); the **diagnostic** tables are computed from the raw per-configuration
/ per-frame energy pools (too large to commit) by
[`scripts/energetic-data-analysis.ipynb`](../../scripts/energetic-data-analysis.ipynb)
§6.4. The figures are written to [`../../figures/`](../../figures/) and interpreted in
[`results/selectivity-assessments-results.md`](../../results/selectivity-assessments-results.md).

| File | Source | Contents |
|------|--------|----------|
| `minimum_binding_energies.csv` | pipeline (§6) | Reported minimum GFN1-xTB binding energy per host-guest pair |
| `ranking.csv` | derived (§7) | Per-guest topology ranking (1 = most negative BE) across the four frameworks |
| `diffBE_CBD_relative.csv` | derived (§7) | CBD-referenced relative binding energy $\Delta BE_{\text{rel}} = BE_{\text{guest}} - BE_{\text{CBD}}$ (cannabinoids) |
| `pharma_host_correlations.csv` | derived (§7) | Per-host Pearson/Spearman of each guest descriptor vs minBE (pharmaceuticals) |
| `spearman_pearson_coeff.csv` | derived (§7) | Within-host Spearman consistency summary per guest set (mean ρ, sign consistency, # significant hosts) |
| `csq_sampling_summary.csv` | diagnostic (§6.4) | Wall-contact statistics for the open-pore csq host (mean/median distance, % within 5/7 Å) |
| `md_convergence_summary.csv` | diagnostic (§6.4) | First-half vs second-half BE drift per (topology, guest, T) MD trajectory |
| `be_outliers.csv` | diagnostic (§6.4) | Raw minima beyond the physical bulk range, flagged as GFN1-xTB fragment misassignment (mechanism annotated) |
| `selectivity_explorer.html` | derived (§7) | Self-contained **interactive** dashboard (Plotly): BE heatmap, descriptor explorer, acid/neutral, host geometry |

### Schema — `minimum_binding_energies.csv`

Long format, one row per host-guest pair:

| Column | Description |
|--------|-------------|
| `host` | Framework name (`MOF-545`, `NU-902`, `PCN-223`, `PCN-225`) |
| `guest` | Guest name or abbreviation (pharmaceuticals by full name, cannabinoids by abbreviation) |
| `minBE` | Reported minimum binding energy (kcal/mol; more negative = stronger) |

> **Coverage.** 136 of the 148 nominal rows (4 hosts × 37 guests). THCV and THCVA
> are pending (no final BE yet) and the cannabinoid acids did not converge in MOF-545 (csq). The loader
> [`pipeline_utils.load_results`](../../scripts/pipeline_utils.py) reconciles the
> mixed guest naming against `040_DESCRIPTORS/guests_descriptors.csv` (keying
> everything by abbreviation) and never imputes missing values.
