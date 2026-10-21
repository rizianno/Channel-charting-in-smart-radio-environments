# AoA ablations: noise calibration

Reference scenario only (S1 joint chart, F0 covariance per sweep + fusion, direct path on, 3x3 codebook) under each noise model, to choose the noise model of the full ablation matrix.

Users: 3200 (in front of both surfaces: 2033; EMS1: 2262; EMS2: 2971). Runtime 13 min.

## Factors
- **S1** (`S1-joint`): One joint 2-D t-SNE chart whose coordinates are the azimuths from EMS1 and EMS2. Users: those in front of both surfaces; anchors: 15% of them. Same structure as the reference runs (results_aoa/coarse_3x3, full_11x11).
- **S2** (`S2-perEMS`): Two separate 1-D t-SNE charts, one per surface, each labelled only with that surface's azimuth. Chart e uses every user in front of surface e (also those behind the other surface); anchors: 15% of them. No coupling between the two angles (no triangulation).
- **S3** (`S3-joint-allUE`): Joint 2-D chart as S1, but t-SNE sees ALL users, including those behind one surface (they have no valid joint label, so they are unlabelled points only). Anchors: the same as S1. Also scored on all users against their geometric azimuth (beyond +-90 deg for users behind a surface, error wrapped to +-180 deg) as the worst case.
- **S4** (`S4-perEMS-frontboth`): Two 1-D per-surface charts as S2, but on exactly the S1 users and anchors: isolates "joint vs separate charts" with identical users.
- **F0** (`F0-cov-per-sweep`): Full covariance of each acquisition: one chart per EMS configuration, then fusion across configurations (global and local quality-weighted), as in the reference.
- **F1** (`F1-stacked-cov`): All acquisitions combined into one feature: the log-covariances of all configurations are stacked and give one chart (no fusion). Control for F2/F3: isolates "one chart from all sweeps" vs "per-sweep charts + fusion". (For K = 1 it equals F0.)
- **F2** (`F2-sweep-change`): Sweep changes: per user, the log-covariance of each configuration minus that user's mean over the configurations, stacked; one chart. Removes, in the log domain, everything that does not change with the surface phases (direct path, static scattering). Joint charts centre over all K^2 configurations; per-surface charts centre along their own surface's sweep with the other surface fixed. Needs K >= 2.
- **F3** (`F3-off-subtract`): OFF subtraction: per configuration, the log-covariance with the surface(s) on minus the log-covariance of an extra acquisition with the surface(s) ideally off (absorbing), stacked; one chart. Joint charts switch both surfaces off (1 extra acquisition); per-surface charts switch off only their own surface (one extra acquisition per phase of the other surface). This is the ON/OFF protocol of RIS channel estimation, used here at covariance level.
- **Direct path** `dir-on` / `dir-off` (ideal removal).
- **Noise `noiseless`**: Exact covariance (no estimation noise), as all earlier runs. A floor of 2*sigma_n^2/sigma_s^2 (sigma_s^2 = 23 dBm, sigma_n^2 = -92 dBm) is added only to keep it invertible; this replaces the old fixed 1e-6, which was 55 dB higher.
- **Noise `est_G0dB`**: Estimated channels: each of the N_real = 500 snapshots is H*g + E, E ~ CN(0, sigma_n^2/sigma_s^2), i.e. one pilot symbol per snapshot (no pilot gain).
- **Noise `est_G30dB`**: Estimated channels as est_G0dB, with a 30 dB pilot processing gain (noise variance on each estimate divided by 1000).
- **Noise `est_G60dB`**: Estimated channels as est_G0dB, with a 60 dB pilot processing gain (noise variance on each estimate divided by 10^6).

## Results at a glance (pooled EMS1+EMS2 median |error| [deg], subset `common`)

| noise | scenario | method | K | acquisitions | pooled median | pooled mean | EMS1 median | EMS2 median |
|---|---|---|---|---|---|---|---|---|
| noiseless | [S1-joint_F0-cov-per-sweep_dir-on](noiseless/S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 3 | 9 | 4.1 | 9.1 | 4.0 | 4.1 |
| noiseless | [S1-joint_F0-cov-per-sweep_dir-on](noiseless/S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 3 | 9 | 4.1 | 9.9 | 4.1 | 4.2 |
| est_G0dB | [S1-joint_F0-cov-per-sweep_dir-on](est_G0dB/S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 3 | 9 | 20.5 | 30.2 | 20.4 | 20.7 |
| est_G0dB | [S1-joint_F0-cov-per-sweep_dir-on](est_G0dB/S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 3 | 9 | 21.2 | 32.0 | 21.1 | 21.2 |
| est_G30dB | [S1-joint_F0-cov-per-sweep_dir-on](est_G30dB/S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 3 | 9 | 5.4 | 12.4 | 6.3 | 4.7 |
| est_G30dB | [S1-joint_F0-cov-per-sweep_dir-on](est_G30dB/S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 3 | 9 | 5.4 | 13.4 | 6.5 | 4.6 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](est_G60dB/S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 3 | 9 | 5.1 | 11.2 | 5.6 | 4.8 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](est_G60dB/S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 3 | 9 | 5.2 | 11.5 | 5.7 | 4.8 |

Each scenario folder has its own README (what changes and which scenarios isolate each factor), `trend.png` (error vs codebook size) and one `K<k>/` folder per codebook size. `comparison_common.png` / `comparison_natural.png` compare all scenarios.
