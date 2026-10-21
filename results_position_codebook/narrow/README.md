# Position estimation vs codebook size: narrow range

Codebook: phase steps 0..40 deg/element (the range of the paper's 11x11 codebook). Users 3200, anchors 480 (15%), noiseless covariance with the 1e-6 regularization, as the paper experiment (Main_Multi_Slot_Fusion); only the codebook changes. Generated 2026-10-21 18:26 by `experiments/Main_Position_CodebookSize.m`.

| K | phase steps per EMS [deg/element] | configurations | local fusion mean / p90 error [m] | global fusion mean [m] | best single mean [m] | typical single mean [m] | local fusion TW / CT |
|---|---|---|---|---|---|---|---|
| 1x1 | 0 | 1 | 13.5 / 29.4 | 13.5 | 13.5 | 13.5 | 0.878 / 0.875 |
| 2x2 | [0 40] | 4 | 9.0 / 19.6 | 10.5 | 11.5 | 12.7 | 0.933 / 0.938 |
| 3x3 | [0 20 40] | 9 | 8.0 / 16.2 | 9.6 | 11.5 | 12.4 | 0.943 / 0.948 |
| 4x4 | [0 13.33 26.67 40] | 16 | 7.8 / 15.4 | 9.4 | 11.1 | 12.4 | 0.943 / 0.947 |
| 5x5 | [0 10 20 30 40] | 25 | 7.5 / 15.1 | 9.1 | 11.2 | 12.3 | 0.948 / 0.952 |

## Definitions
- Errors are over the unlabeled users. `typical_single` = average over the K^2 single-configuration charts (what one random configuration gives); `best_single` = the configuration with the lowest mean error (oracle, as `evaluate_best_codebook`); `global_fusion` / `local_fusion` = the repo's `global_weighting` / `local_weighting` over the K^2 charts. Fusion gain = typical_single - fused error.
- TW / CT = trustworthiness / continuity of the chart vs the true positions (K = 5% of the users), mean over the unlabeled users (p10 also in the CSVs). Note: the repo's `compute_all_metrics` and `run_channel_charting` store these with swapped names (its `TW_*` is column 2 of `drscore_pointwise` = continuity); here they are named correctly. The fusion weights use the repo's quality (column 2) unchanged.
- Files: `trend.png` / `trend.csv` (error, TW, CT vs codebook size), `K<k>/summary.csv`, `K<k>/per_configuration.csv` (each single configuration), `K<k>/cdf.png`.
