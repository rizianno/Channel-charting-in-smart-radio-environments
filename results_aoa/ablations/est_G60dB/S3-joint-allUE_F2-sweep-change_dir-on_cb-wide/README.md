# S3-joint-allUE_F2-sweep-change_dir-on_cb-wide

Noise model: `est_G60dB`. Codebook sizes: 2x2, 3x3, 4x4, 5x5. Generated 2026-10-21 19:50 by `experiments/Main_AoA_Ablations.m`.

## What this simulation does
- **Users and charts (S3).** Joint 2-D chart as S1, but t-SNE sees ALL users, including those behind one surface (they have no valid joint label, so they are unlabelled points only). Anchors: the same as S1. Also scored on all users against their geometric azimuth (beyond +-90 deg for users behind a surface, error wrapped to +-180 deg) as the worst case.
  - chart `S3`: 3200 users, 305 anchors, label = azimuth from EMS1 and EMS2
- **Input feature (F2).** Sweep changes: per user, the log-covariance of each configuration minus that user's mean over the configurations, stacked; one chart. Removes, in the log domain, everything that does not change with the surface phases (direct path, static scattering). Joint charts centre over all K^2 configurations; per-surface charts centre along their own surface's sweep with the other surface fixed. Needs K >= 2.
- **Direct path.** **on**: the direct BS-UE channel (LoS + one bounce, `H_LOS_1bounce`) is included, as in the reference.
- **Noise (est_G60dB).** Estimated channels as est_G0dB, with a 60 dB pilot processing gain (noise variance on each estimate divided by 10^6).
- **Codebook range (wide).** phase steps spread evenly over -180..180 deg/element (beams spread evenly over all directions).
- **Unchanged everywhere:** dataset, surface geometry, t-SNE (`tsne_sem`, perplexity 150), Log-Euclidean distance, 15% supervision, fusion rule (where there is fusion).

| K | phase steps per EMS [deg/element] | acquisitions |
|---|---|---|
| 2x2 | [-180 0] | 4 |
| 3x3 | [-120 0 120] | 9 |
| 4x4 | [-180 -90 0 90] | 16 |
| 5x5 | [-144 -72 0 72 144] | 25 |

## How it differs from the other scenarios
- **vs the reference** (S1 joint, F0 per-sweep covariance + fusion, direct on, narrow codebook): users/charts S3-joint-allUE instead of S1-joint; feature F2-sweep-change instead of F0-cov-per-sweep; wide codebook (-180..180 deg/element) instead of narrow (0..40 deg/element).
- **Scenarios that differ in exactly one factor** (compare with these to isolate it):
  - `S3-joint-allUE_F2-sweep-change_dir-on`: narrow codebook (0..40 deg/element) instead of wide (-180..180 deg/element)
  - `S1-joint_F2-sweep-change_dir-on_cb-wide`: users/charts S1-joint instead of S3-joint-allUE
  - `S2-perEMS_F2-sweep-change_dir-on_cb-wide`: users/charts S2-perEMS instead of S3-joint-allUE
  - `S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide`: feature F0-cov-per-sweep instead of F2-sweep-change
  - `S3-joint-allUE_F1-stacked-cov_dir-on_cb-wide`: feature F1-stacked-cov instead of F2-sweep-change
  - `S3-joint-allUE_F2-sweep-change_dir-off_cb-wide`: direct path off instead of on
  - `S3-joint-allUE_F3-off-subtract_dir-on_cb-wide`: feature F3-off-subtract instead of F2-sweep-change
- Older reference (no noise, old fixed 1e-6 floor): `results_aoa/coarse_3x3` (3x3) and `results_aoa/full_11x11`.

## Files
- `K<k>/summary.csv`: median / mean / p90 / RMSE of |error| [deg] per method, user subset and EMS (EMS1, EMS2, pooled = both sets of errors together, joint = sqrt(e1^2+e2^2) for joint charts).
- `K<k>/per_user_*.csv`: per-user truth, estimate and error of the main method, with flags (front of each surface, direct path present, reachable from each surface, anchor).
- `K<k>/cdf.png`: CDFs per EMS, pooled and (joint charts) joint. `K<k>/map_per_ems*.png`: error maps per EMS (absolute 0-30 deg and signed +-30 deg; S3 also all users, 0-90 deg).
- `trend.csv`, `trend.png`: error vs codebook size (number of acquisitions).

## Definitions
- **Error:** estimated minus true azimuth, wrapped to [-180, 180) deg; azimuth measured from each surface normal, positive toward +y.
- **Subsets:** `common` = in front of both surfaces and not an anchor in any scenario (identical users everywhere: use it to compare scenarios); `natural` = every user this scenario estimates in front of that surface; `reach` = natural with a non-zero channel to that surface; `all` (S3) = every user, also behind the surface (worst case).
- **Method:** `chart` = the single chart built from all acquisitions (no fusion); its improvement with K comes only from the richer combined feature.
