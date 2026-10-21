# Position estimation vs codebook size: wide range

Codebook: phase steps spread evenly over -180..180 deg/element (beams over all directions). Users 3200, anchors 480 (15%), noiseless covariance with the 1e-6 regularization, as the paper experiment (Main_Multi_Slot_Fusion); only the codebook changes. Generated 2026-10-21 18:30 by `experiments/Main_Position_CodebookSize.m`.

| K | phase steps per EMS [deg/element] | configurations | local fusion mean / p90 error [m] | global fusion mean [m] | best single mean [m] | typical single mean [m] | local fusion TW / CT |
|---|---|---|---|---|---|---|---|
| 1x1 | 0 | 1 | 13.5 / 29.4 | 13.5 | 13.5 | 13.5 | 0.878 / 0.875 |
| 2x2 | [-180 0] | 4 | 10.4 / 24.4 | 11.9 | 13.4 | 14.0 | 0.913 / 0.917 |
| 3x3 | [-120 0 120] | 9 | 9.2 / 22.0 | 10.7 | 11.1 | 13.4 | 0.930 / 0.934 |
| 4x4 | [-180 -90 0 90] | 16 | 9.1 / 20.2 | 10.8 | 12.3 | 13.7 | 0.931 / 0.934 |
| 5x5 | [-144 -72 0 72 144] | 25 | 9.0 / 20.8 | 10.9 | 11.6 | 13.9 | 0.934 / 0.937 |

## Definitions
- Errors are over the unlabeled users. `typical_single` = average over the K^2 single-configuration charts (what one random configuration gives); `best_single` = the configuration with the lowest mean error (oracle, as `evaluate_best_codebook`); `global_fusion` / `local_fusion` = the repo's `global_weighting` / `local_weighting` over the K^2 charts. Fusion gain = typical_single - fused error.
- TW / CT = trustworthiness / continuity of the chart vs the true positions (K = 5% of the users), mean over the unlabeled users (p10 also in the CSVs). Note: the repo's `compute_all_metrics` and `run_channel_charting` store these with swapped names (its `TW_*` is column 2 of `drscore_pointwise` = continuity); here they are named correctly. The fusion weights use the repo's quality (column 2) unchanged.
- Files: `trend.png` / `trend.csv` (error, TW, CT vs codebook size), `K<k>/summary.csv`, `K<k>/per_configuration.csv` (each single configuration), `K<k>/cdf.png`.
