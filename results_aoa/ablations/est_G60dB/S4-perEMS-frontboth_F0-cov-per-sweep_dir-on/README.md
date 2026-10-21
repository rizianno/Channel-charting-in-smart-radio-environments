# S4-perEMS-frontboth_F0-cov-per-sweep_dir-on

Noise model: `est_G60dB`. Codebook sizes: 3x3. Generated 2026-10-21 17:27 by `experiments/Main_AoA_Ablations.m`.

## What this simulation does
- **Users and charts (S4).** Two 1-D per-surface charts as S2, but on exactly the S1 users and anchors: isolates "joint vs separate charts" with identical users.
  - chart `S4e1`: 2033 users, 305 anchors, label = azimuth from EMS1
  - chart `S4e2`: 2033 users, 305 anchors, label = azimuth from EMS2
- **Input feature (F0).** Full covariance of each acquisition: one chart per EMS configuration, then fusion across configurations (global and local quality-weighted), as in the reference.
- **Direct path.** **on**: the direct BS-UE channel (LoS + one bounce, `H_LOS_1bounce`) is included, as in the reference.
- **Noise (est_G60dB).** Estimated channels as est_G0dB, with a 60 dB pilot processing gain (noise variance on each estimate divided by 10^6).
- **Unchanged everywhere:** dataset, surface geometry, t-SNE (`tsne_sem`, perplexity 150), Log-Euclidean distance, 15% supervision, fusion rule (where there is fusion).

| K | phase gradients per EMS [deg] | acquisitions |
|---|---|---|
| 3x3 | [0 20 40] | 9 |

## How it differs from the other scenarios
- **vs the reference** (S1 joint, F0 per-sweep covariance + fusion, direct on): users/charts S4-perEMS-frontboth instead of S1-joint.
- **Scenarios that differ in exactly one factor** (compare with these to isolate it):
  - `S1-joint_F0-cov-per-sweep_dir-on`: users/charts S1-joint instead of S4-perEMS-frontboth
  - `S2-perEMS_F0-cov-per-sweep_dir-on`: users/charts S2-perEMS instead of S4-perEMS-frontboth
  - `S3-joint-allUE_F0-cov-per-sweep_dir-on`: users/charts S3-joint-allUE instead of S4-perEMS-frontboth
- Older reference (no noise, old fixed 1e-6 floor): `results_aoa/coarse_3x3` (3x3) and `results_aoa/full_11x11`.

## Files
- `K<k>/summary.csv`: median / mean / p90 / RMSE of |error| [deg] per method, user subset and EMS (EMS1, EMS2, pooled = both sets of errors together, joint = sqrt(e1^2+e2^2) for joint charts).
- `K<k>/per_user_*.csv`: per-user truth, estimate and error of the main method, with flags (front of each surface, direct path present, reachable from each surface, anchor).
- `K<k>/cdf.png`: CDFs per EMS, pooled and (joint charts) joint. `K<k>/map_per_ems*.png`: error maps per EMS (absolute 0-30 deg and signed +-30 deg; S3 also all users, 0-90 deg).
- `trend.csv`, `trend.png`: error vs codebook size (number of acquisitions).

## Definitions
- **Error:** estimated minus true azimuth, wrapped to [-180, 180) deg; azimuth measured from each surface normal, positive toward +y.
- **Subsets:** `common` = in front of both surfaces and not an anchor in any scenario (identical users everywhere: use it to compare scenarios); `natural` = every user this scenario estimates in front of that surface; `reach` = natural with a non-zero channel to that surface; `all` (S3) = every user, also behind the surface (worst case).
- **Methods:** `typical_single` = average over the configurations of a single-configuration chart; `best_single` = best configuration (oracle, chosen on the scored users); `global_*` / `local_*` = fusion of all configurations. `*_orig` weights each chart by its quality measured against the **true** labels of all users, exactly as the original code (so it uses the unlabeled users' ground truth); `*_free` measures the same quality against the input dissimilarities instead (no labels needed). Fusion gain = typical_single minus the fused error.
- Note: the original code calls this quality "TW" but it is column 2 of `drscore_pointwise`, i.e. the continuity; kept as is for comparability.
