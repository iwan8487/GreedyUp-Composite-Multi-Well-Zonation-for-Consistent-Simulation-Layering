# GreedyUp: Composite Multi-Well Zonation for Consistent Simulation Layering

GreedyUp builds **one laterally consistent set of simulation layers for several wells directly from core or log data**, before any 3-D geocellular model exists. It extends a greedy, threshold-based single-well zonation algorithm to a composite formulation. The local heterogeneity of log2-permeability in every well is mapped to a common relative stratigraphic coordinate and aggregated across wells. The aggregated curve is segmented once, and the resulting boundaries are transferred proportionally to each well's measured depth.

The repository contains:

* an interactive, single-file browser application (no installation, works offline);
* a Python reference implementation that reproduces the results of the accompanying paper;
* the result tables and figures of the paper.

> **Paper.** Alfadli, M. K., & Natasia, N. *Composite multi-well zonation for laterally consistent simulation layering: dynamic validation and the decisive role of the permeability-averaging rule.* [Journal / preprint, year, DOI — to be completed]
>
> **Previous work (single-well algorithm).** Natasia, N., & Alfadli, M. K. (2025). *Data-Driven Zonation and Upscaling of Permeability Profiles Utilizing Log-Space Homogenization for Reservoir Modelling* [Preprint]. SSRN. https://doi.org/10.2139/ssrn.5295005

![Independent versus composite zonation]({{artifact:art_a17d745e-b210-4b2b-ba53-503694452d20}})

*Independent zonation of four wells (grey) places boundaries inconsistently (mean pairwise agreement 34%). The composite method (blue) yields one shared 33-zone framework.*

---

## Contents

1. [Method in brief](#method-in-brief)
2. [Main results](#main-results)
3. [Repository structure](#repository-structure)
4. [Interactive application](#interactive-application)
5. [Python implementation](#python-implementation)
6. [Input data format](#input-data-format)
7. [Reproducing the paper](#reproducing-the-paper)
8. [Known issues and differences between implementations](#known-issues-and-differences-between-implementations)
9. [Citation, license and contact](#citation-license-and-contact)

---

## Method in brief

| Step | Operation | Default |
|---|---|---|
| 1 | Map each well to the relative coordinate t = (MD − MD_top) / (MD_base − MD_top) and resample log2K onto N_T points (clamped linear interpolation) | N_T = 1,000 |
| 2 | Compute the local heterogeneity H_w(t), i.e. the moving-window standard deviation of log2K | window f_W = 2% of the interval |
| 3 | Aggregate across wells: H(t) = max / P90 / mean of H_w(t) | max |
| 4 | Greedy segmentation of H(t): grow each zone while the RMSE from the zone mean is below E_z = c · sd(H), within L_min ≤ width ≤ L_max | c = 0.5, L_min = 5, L_max = 40 samples |
| 5 | Transfer the boundaries to each well: MD_z,w = MD_top,w + t_z (MD_base,w − MD_top,w) | proportional |
| 6 | Zone properties from the native samples: porosity is the arithmetic mean; permeability is geometric, or arithmetic k_h with harmonic k_v | see [Known issues](#known-issues-and-differences-between-implementations) |
| 7 | Parameter search over c × f_W (8 × 8) and Pareto front of shared zone count versus worst-well Pearson r | 64 combinations |
| 8 | Dynamic validation with a 2-D single-phase finite-volume model against a 382-layer proportional benchmark, for all producer-well combinations | 400 columns, 2,000 days |

Every well receives the same number of zones, and the relative thickness of each zone is the same in every well. This is a property of the construction, valid under the assumption that layering is proportional between the two picked surfaces.

## Main results

Four cored wells (HW-7, HW-29, HW-32, HW-28; 0.1-ft core data, H2 interval):

| Result | Value |
|---|---|
| Independent per-well zonation | 33–38 zones per well; boundary agreement 14–75% (mean 34%) |
| Composite zonation (default) | 33 shared zones (−91.4% vs. 382-layer benchmark); worst-well r = 0.835 |
| Zones controlled by L_max rather than by heterogeneity | 18 of 33 (55%); 100% at 25 zones, i.e. identical to uniform layering |
| Parameter search | 64 combinations, 25–120 zones, 19 Pareto-optimal |
| Producer-pressure error, arithmetic k_h / harmonic k_v | ≤ 1.16% of drawdown for every grid; ≤ 0.29% from 61 zones |
| Composite vs. uniform layering, geometric averaging | 20–72% lower error from 29 zones (Wilcoxon p ≤ 0.002) |
| Composite vs. uniform layering, arithmetic k_h / harmonic k_v | no difference up to 43 zones (p = 0.30–0.72); 38–62% lower error from 61 zones (p = 6.1 × 10⁻⁵) |
| Mechanism under geometric averaging | the error ranks exactly with the loss of column flow capacity (Spearman ρ = 1.00, 16 grids) |
| Generalisation to all 17 wells | worst-well r at default settings 0.73; zones needed for r ≥ 0.85: 49 (vs. 35 for 4 wells) |

![Effect of the permeability-averaging rule]({{artifact:art_a9e6a8d7-adef-4f08-9ad9-de8f1043fe74}})

## Repository structure

```
.
├── README.md
├── CITATION.cff
├── requirements.txt
├── app/
│   └── consistent_grid_explorer.html     # interactive application (single file, offline)
├── python/
│   ├── zonation_core.py                  # zonation, upscaling and flow-simulation library
│   └── reproduce_results.py              # command-line pipeline for the paper's results
├── data/                                 # LAS files and well tops (see "Input data format")
├── results/                              # CSV tables behind every figure and table of the paper
├── figures/                              # paper figures (PNG, 300 dpi)
└── docs/
    └── APP_PROMPT.md                     # specification prompt for maintaining or extending the app
```

## Interactive application

Open `app/consistent_grid_explorer.html` in a recent desktop browser (Chrome, Edge or Firefox). No server, installation or internet connection is required. The application starts with the four-well demonstration data set.

Typical workflow:

1. **Data source.** Keep the demo data, or add LAS files (multiple selection) and a well-top file, choose the top and base surfaces, and click *Load data & recompute*. Depth, permeability and porosity curves are detected from their mnemonics (`DEPT`/`DEPTH`/`MD`, `CORE_PERM`/`PERM`/`K`…, `CORE_POR`/`PHI`/`POR`…).
2. **Well map.** Select the wells of the cross-section on the map and click *Use this section*. Wells are ordered west to east by X, and the section distance is the cumulative path length between consecutive wells.
3. **Shared zonation.** Set E_z (× composite σ), the window (% of interval), L_min, L_max and the aggregator (Max, P90, Mean). The panels show the composite heterogeneity curve, the consistent cross-section, the zone-thickness pattern and the layer-count comparison with the fine grid.
4. **Parameter optimization.** *Auto-optimize (Grid Search)* evaluates the E_z × window grid and plots the Pareto front of zone count versus worst-well fidelity.
5. **Grid comparison set and flow simulation.** Add up to 8 parameter sets to the comparison set, choose the producer, rate, duration, lateral resolution and lateral extrapolation, and click *Run flow simulation*. Each grid is compared with the fine-grid benchmark.
6. **Export.** Figures can be exported as PNG (1×–6×, 150–1,200 dpi) or SVG. The upscaling-ratio table can be exported as CSV, and the current parameter set as JSON.

## Python implementation

```bash
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
cd python

# exact reproduction of the paper (uses the data set embedded in the application)
python reproduce_results.py --app-html ../app/consistent_grid_explorer.html --out ../results

# from LAS files and a well-top file (also runs the 17-well generalisation)
python reproduce_results.py --las-dir ../data/las --welltops ../data/welltop.txt --out ../results_las

# fast smoke test (3 grids, 80 columns, 12 time steps; about 15 s)
python reproduce_results.py --las-dir ../data/las --welltops ../data/welltop.txt --out ../results_quick --quick
```

Options: `--wells HW-7,HW-29,HW-32,HW-28` selects the section wells; `--skip-flow` runs only the static part. The full run (2 averaging rules × 17 grids × 15 producer sets on a 400-column section) took about 5 minutes on a 12-core laptop and ran within about 1 GB of free memory. The `--app-html` run reproduces the flow results of the paper to machine precision.

Outputs written to `--out`: `gridsearch.csv`, `independent_zonation_agreement.csv`, `zones_needed_<N>wells.csv`, `flow_batch.csv` (one row per grid × producer set × well) and `flow_summary.csv` (normalised RMSE, advantage over uniform layering and Wilcoxon p-value per zone count and averaging rule).

Library use:

```python
import zonation_core as zc
res = zc.run_composite(wells, order, wf=0.02, ezf=0.5, lmin=5, lmax=40, agg="max", avg="arith")
res["nZones"], res["tedges"], res["zt"]["HW-7"][0]      # zone count, boundaries in t, zone table
```

`wells` is a dict `{name: {"depth", "k", "por", "h2_top", "h2_base"}}` with NumPy arrays.

## Input data format

**LAS files** (LAS 2.0): one file per well, named after the well (for example `HW-7.las`), with a depth curve in feet, core permeability in mD and core porosity as a fraction. Null values (for example −999.25) are handled.

**Well-top file**: delimited text (tab, comma or whitespace) with one row per well and surface:

```
Well identifier	Surface	X	Y	Z	MD
HW-7	H2 Top	34526	80000	-8162.46	8529.46
HW-7	H2 Base	34526	80000	-8311	8678
```

## Reproducing the paper

The tables below are the result files archived with the paper (folder `results/`); the script regenerates equivalent files under the names listed above.

| Paper item | File in `results/` | Notes |
|---|---|---|
| Fig. 1, independent vs. composite zonation | `independent_zonation_agreement.csv` | ±5 samples of t (±0.5% of interval) |
| Fig. 4, Table 4, Pareto front | `gridsearch_replication.csv`, `heterogeneity_metrics.csv` | 64 combinations |
| Fig. 5, Fig. 6, Table 5, flow validation | `batch_taligned.csv`, `retest_averaging_summary.csv`, `kh_deficit_geometric.csv` | both averaging rules |
| Fig. 7, Table 6, 17 wells | `generalization_17wells.csv`, `gridsearch_4v17_aggregators.csv` | LAS preprocessing |
| Table 2, well statistics | `table2_corrected.csv` | analysed H2 interval only |

## Known issues and differences between implementations

These points matter when you compare numbers from the browser application with the paper.

1. **Fine-grid benchmark thickness (application).** `buildFineZoneTableForWell` sets the layer thickness to `depth[e] − depth[s]`, which omits one sample interval per layer. The benchmark is then 26–28% thinner than the real interval (HW-7: 110.3 ft instead of 148.5 ft), and its pressure errors do not converge with zone count. The Python implementation uses the exact thickness and t-aligned (proportional) benchmark layers.
2. **Permeability averaging (application).** The application uses the geometric mean for both horizontal and vertical transmissibility. The paper shows that this rule inflates both the errors and the apparent benefit of heterogeneity-aware layering; the Python implementation also offers arithmetic k_h with harmonic k_v (`avg="arith"`).
3. **P90 aggregator.** P90 selects the ⌈0.9 N⌉-th value and equals the maximum for N ≤ 9 wells.
4. **Configuration import.** The application exports the parameter set as JSON but cannot yet import it.
5. **Data paths.** The demo data embedded in the application (used for all four-well numbers in the paper) differ slightly from a fresh import of the LAS files. The sample counts of HW-7 and HW-32 differ by 1–2 samples and some interior values of HW-29 differ. With the LAS path, the default grid has worst-well r = 0.848 instead of 0.835, and the search gives 15 instead of 19 Pareto points. Use `--app-html` for exact reproduction.

`docs/APP_PROMPT.md` contains a specification prompt that lists these fixes for the next version of the application.

## Citation, license and contact

If you use this software, please cite the paper and the software.

* License: Redistribution of the core data in `data/` requires the permission of the data owner.
* Contact: M. Kurniawan Alfadli, m.kurniawan@unpad.ac.id , Universitas Padjadjaran.

---
