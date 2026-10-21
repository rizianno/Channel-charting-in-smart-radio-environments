# AoA ablations, noise model `est_G60dB`

Noise: Estimated channels as est_G0dB, with a 60 dB pilot processing gain (noise variance on each estimate divided by 10^6).

Users: 3200 (in front of both surfaces: 2033; EMS1: 2262; EMS2: 2971). Runtime 81 min.

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
- **Noise `est_G60dB`**: Estimated channels as est_G0dB, with a 60 dB pilot processing gain (noise variance on each estimate divided by 10^6).

## Results at a glance (pooled EMS1+EMS2 median |error| [deg], subset `common`)

| noise | scenario | method | K | acquisitions | pooled median | pooled mean | EMS1 median | EMS2 median |
|---|---|---|---|---|---|---|---|---|
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 1 | 1 | 5.7 | 12.9 | 6.4 | 5.1 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 1 | 1 | 5.7 | 12.9 | 6.4 | 5.1 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 2 | 4 | 5.5 | 11.6 | 6.2 | 4.9 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 2 | 4 | 5.5 | 11.9 | 6.3 | 4.9 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 3 | 9 | 5.1 | 11.2 | 5.6 | 4.8 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 3 | 9 | 5.2 | 11.5 | 5.7 | 4.8 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 4 | 16 | 5.2 | 11.2 | 5.7 | 4.6 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 4 | 16 | 5.3 | 11.6 | 5.8 | 4.7 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 5 | 25 | 5.1 | 11.0 | 5.4 | 4.7 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 5 | 25 | 5.1 | 11.5 | 5.5 | 4.7 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-off](S1-joint_F0-cov-per-sweep_dir-off/) | local_orig | 1 | 1 | 40.0 | 50.1 | 37.7 | 41.6 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-off](S1-joint_F0-cov-per-sweep_dir-off/) | local_free | 1 | 1 | 40.0 | 50.1 | 37.7 | 41.6 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-off](S1-joint_F0-cov-per-sweep_dir-off/) | local_orig | 2 | 4 | 32.5 | 37.5 | 31.3 | 33.8 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-off](S1-joint_F0-cov-per-sweep_dir-off/) | local_free | 2 | 4 | 35.2 | 40.0 | 34.0 | 36.4 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-off](S1-joint_F0-cov-per-sweep_dir-off/) | local_orig | 3 | 9 | 31.3 | 34.5 | 30.2 | 32.3 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-off](S1-joint_F0-cov-per-sweep_dir-off/) | local_free | 3 | 9 | 34.6 | 37.4 | 33.6 | 36.1 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-off](S1-joint_F0-cov-per-sweep_dir-off/) | local_orig | 4 | 16 | 30.6 | 34.1 | 31.2 | 29.8 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-off](S1-joint_F0-cov-per-sweep_dir-off/) | local_free | 4 | 16 | 33.7 | 37.1 | 33.9 | 33.7 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-off](S1-joint_F0-cov-per-sweep_dir-off/) | local_orig | 5 | 25 | 30.5 | 33.6 | 30.2 | 30.8 |
| est_G60dB | [S1-joint_F0-cov-per-sweep_dir-off](S1-joint_F0-cov-per-sweep_dir-off/) | local_free | 5 | 25 | 33.8 | 36.7 | 33.3 | 34.4 |
| est_G60dB | [S1-joint_F1-stacked-cov_dir-on](S1-joint_F1-stacked-cov_dir-on/) | chart | 2 | 4 | 5.8 | 12.3 | 7.2 | 4.8 |
| est_G60dB | [S1-joint_F1-stacked-cov_dir-on](S1-joint_F1-stacked-cov_dir-on/) | chart | 3 | 9 | 5.2 | 11.7 | 5.9 | 4.7 |
| est_G60dB | [S1-joint_F1-stacked-cov_dir-on](S1-joint_F1-stacked-cov_dir-on/) | chart | 4 | 16 | 5.1 | 11.6 | 5.8 | 4.6 |
| est_G60dB | [S1-joint_F1-stacked-cov_dir-on](S1-joint_F1-stacked-cov_dir-on/) | chart | 5 | 25 | 5.0 | 11.5 | 5.4 | 4.7 |
| est_G60dB | [S1-joint_F1-stacked-cov_dir-off](S1-joint_F1-stacked-cov_dir-off/) | chart | 2 | 4 | 31.1 | 40.6 | 28.0 | 34.1 |
| est_G60dB | [S1-joint_F1-stacked-cov_dir-off](S1-joint_F1-stacked-cov_dir-off/) | chart | 3 | 9 | 26.8 | 36.5 | 27.5 | 26.5 |
| est_G60dB | [S1-joint_F1-stacked-cov_dir-off](S1-joint_F1-stacked-cov_dir-off/) | chart | 4 | 16 | 24.4 | 39.8 | 22.6 | 26.7 |
| est_G60dB | [S1-joint_F1-stacked-cov_dir-off](S1-joint_F1-stacked-cov_dir-off/) | chart | 5 | 25 | 21.1 | 36.0 | 18.1 | 24.2 |
| est_G60dB | [S1-joint_F2-sweep-change_dir-on](S1-joint_F2-sweep-change_dir-on/) | chart | 2 | 4 | 37.4 | 50.9 | 35.9 | 38.8 |
| est_G60dB | [S1-joint_F2-sweep-change_dir-on](S1-joint_F2-sweep-change_dir-on/) | chart | 3 | 9 | 38.2 | 50.9 | 36.1 | 39.6 |
| est_G60dB | [S1-joint_F2-sweep-change_dir-on](S1-joint_F2-sweep-change_dir-on/) | chart | 4 | 16 | 39.7 | 48.5 | 36.3 | 43.4 |
| est_G60dB | [S1-joint_F2-sweep-change_dir-on](S1-joint_F2-sweep-change_dir-on/) | chart | 5 | 25 | 40.9 | 50.8 | 38.7 | 42.8 |
| est_G60dB | [S1-joint_F2-sweep-change_dir-off](S1-joint_F2-sweep-change_dir-off/) | chart | 2 | 4 | 35.3 | 43.9 | 33.0 | 37.8 |
| est_G60dB | [S1-joint_F2-sweep-change_dir-off](S1-joint_F2-sweep-change_dir-off/) | chart | 3 | 9 | 29.4 | 38.3 | 28.4 | 30.4 |
| est_G60dB | [S1-joint_F2-sweep-change_dir-off](S1-joint_F2-sweep-change_dir-off/) | chart | 4 | 16 | 26.6 | 37.4 | 24.0 | 28.7 |
| est_G60dB | [S1-joint_F2-sweep-change_dir-off](S1-joint_F2-sweep-change_dir-off/) | chart | 5 | 25 | 23.2 | 37.6 | 21.0 | 25.4 |
| est_G60dB | [S1-joint_F3-off-subtract_dir-on](S1-joint_F3-off-subtract_dir-on/) | chart | 1 | 2 | 45.9 | 56.5 | 45.0 | 46.7 |
| est_G60dB | [S1-joint_F3-off-subtract_dir-on](S1-joint_F3-off-subtract_dir-on/) | chart | 2 | 5 | 39.2 | 50.7 | 40.2 | 38.5 |
| est_G60dB | [S1-joint_F3-off-subtract_dir-on](S1-joint_F3-off-subtract_dir-on/) | chart | 3 | 10 | 37.0 | 49.7 | 36.1 | 37.9 |
| est_G60dB | [S1-joint_F3-off-subtract_dir-on](S1-joint_F3-off-subtract_dir-on/) | chart | 4 | 17 | 36.3 | 49.8 | 34.5 | 37.9 |
| est_G60dB | [S1-joint_F3-off-subtract_dir-on](S1-joint_F3-off-subtract_dir-on/) | chart | 5 | 26 | 38.5 | 50.3 | 39.6 | 37.5 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-on](S2-perEMS_F0-cov-per-sweep_dir-on/) | local_orig | 1 | 1 | 5.9 | 14.4 | 6.5 | 5.4 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-on](S2-perEMS_F0-cov-per-sweep_dir-on/) | local_free | 1 | 1 | 5.9 | 14.4 | 6.5 | 5.4 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-on](S2-perEMS_F0-cov-per-sweep_dir-on/) | local_orig | 2 | 4 | 5.7 | 11.9 | 6.1 | 5.4 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-on](S2-perEMS_F0-cov-per-sweep_dir-on/) | local_free | 2 | 4 | 5.9 | 13.0 | 6.1 | 5.5 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-on](S2-perEMS_F0-cov-per-sweep_dir-on/) | local_orig | 3 | 9 | 5.8 | 11.1 | 6.2 | 5.4 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-on](S2-perEMS_F0-cov-per-sweep_dir-on/) | local_free | 3 | 9 | 6.0 | 12.4 | 6.4 | 5.7 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-on](S2-perEMS_F0-cov-per-sweep_dir-on/) | local_orig | 4 | 16 | 5.8 | 11.1 | 6.0 | 5.5 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-on](S2-perEMS_F0-cov-per-sweep_dir-on/) | local_free | 4 | 16 | 5.9 | 12.5 | 6.1 | 5.7 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-on](S2-perEMS_F0-cov-per-sweep_dir-on/) | local_orig | 5 | 25 | 5.8 | 10.8 | 6.0 | 5.4 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-on](S2-perEMS_F0-cov-per-sweep_dir-on/) | local_free | 5 | 25 | 6.0 | 12.3 | 6.3 | 5.7 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-off](S2-perEMS_F0-cov-per-sweep_dir-off/) | local_orig | 1 | 1 | 39.6 | 45.5 | 35.7 | 42.8 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-off](S2-perEMS_F0-cov-per-sweep_dir-off/) | local_free | 1 | 1 | 39.6 | 45.5 | 35.7 | 42.8 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-off](S2-perEMS_F0-cov-per-sweep_dir-off/) | local_orig | 2 | 4 | 33.7 | 38.7 | 32.6 | 35.4 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-off](S2-perEMS_F0-cov-per-sweep_dir-off/) | local_free | 2 | 4 | 37.0 | 41.0 | 35.7 | 38.5 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-off](S2-perEMS_F0-cov-per-sweep_dir-off/) | local_orig | 3 | 9 | 32.0 | 35.9 | 31.6 | 32.7 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-off](S2-perEMS_F0-cov-per-sweep_dir-off/) | local_free | 3 | 9 | 35.8 | 38.8 | 35.4 | 36.2 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-off](S2-perEMS_F0-cov-per-sweep_dir-off/) | local_orig | 4 | 16 | 33.3 | 36.1 | 34.0 | 32.1 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-off](S2-perEMS_F0-cov-per-sweep_dir-off/) | local_free | 4 | 16 | 36.4 | 38.9 | 36.8 | 35.5 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-off](S2-perEMS_F0-cov-per-sweep_dir-off/) | local_orig | 5 | 25 | 32.0 | 35.7 | 31.7 | 32.2 |
| est_G60dB | [S2-perEMS_F0-cov-per-sweep_dir-off](S2-perEMS_F0-cov-per-sweep_dir-off/) | local_free | 5 | 25 | 35.8 | 38.6 | 35.7 | 35.8 |
| est_G60dB | [S2-perEMS_F1-stacked-cov_dir-on](S2-perEMS_F1-stacked-cov_dir-on/) | chart | 2 | 4 | 4.9 | 11.9 | 4.7 | 5.0 |
| est_G60dB | [S2-perEMS_F1-stacked-cov_dir-on](S2-perEMS_F1-stacked-cov_dir-on/) | chart | 3 | 9 | 4.7 | 11.3 | 4.4 | 5.0 |
| est_G60dB | [S2-perEMS_F1-stacked-cov_dir-on](S2-perEMS_F1-stacked-cov_dir-on/) | chart | 4 | 16 | 4.6 | 11.1 | 4.3 | 5.0 |
| est_G60dB | [S2-perEMS_F1-stacked-cov_dir-on](S2-perEMS_F1-stacked-cov_dir-on/) | chart | 5 | 25 | 4.6 | 11.4 | 4.1 | 5.1 |
| est_G60dB | [S2-perEMS_F1-stacked-cov_dir-off](S2-perEMS_F1-stacked-cov_dir-off/) | chart | 2 | 4 | 28.9 | 39.3 | 25.1 | 33.9 |
| est_G60dB | [S2-perEMS_F1-stacked-cov_dir-off](S2-perEMS_F1-stacked-cov_dir-off/) | chart | 3 | 9 | 23.2 | 32.7 | 22.1 | 24.1 |
| est_G60dB | [S2-perEMS_F1-stacked-cov_dir-off](S2-perEMS_F1-stacked-cov_dir-off/) | chart | 4 | 16 | 23.7 | 38.6 | 21.9 | 26.1 |
| est_G60dB | [S2-perEMS_F1-stacked-cov_dir-off](S2-perEMS_F1-stacked-cov_dir-off/) | chart | 5 | 25 | 21.0 | 35.4 | 19.0 | 23.3 |
| est_G60dB | [S2-perEMS_F2-sweep-change_dir-on](S2-perEMS_F2-sweep-change_dir-on/) | chart | 2 | 4 | 45.0 | 49.9 | 35.2 | 49.9 |
| est_G60dB | [S2-perEMS_F2-sweep-change_dir-on](S2-perEMS_F2-sweep-change_dir-on/) | chart | 3 | 9 | 43.5 | 48.5 | 34.6 | 50.8 |
| est_G60dB | [S2-perEMS_F2-sweep-change_dir-on](S2-perEMS_F2-sweep-change_dir-on/) | chart | 4 | 16 | 43.8 | 48.4 | 31.9 | 51.0 |
| est_G60dB | [S2-perEMS_F2-sweep-change_dir-on](S2-perEMS_F2-sweep-change_dir-on/) | chart | 5 | 25 | 45.5 | 53.9 | 38.4 | 54.3 |
| est_G60dB | [S2-perEMS_F2-sweep-change_dir-off](S2-perEMS_F2-sweep-change_dir-off/) | chart | 2 | 4 | 43.2 | 47.7 | 37.4 | 49.0 |
| est_G60dB | [S2-perEMS_F2-sweep-change_dir-off](S2-perEMS_F2-sweep-change_dir-off/) | chart | 3 | 9 | 35.2 | 44.1 | 24.5 | 47.3 |
| est_G60dB | [S2-perEMS_F2-sweep-change_dir-off](S2-perEMS_F2-sweep-change_dir-off/) | chart | 4 | 16 | 33.7 | 43.1 | 23.4 | 45.9 |
| est_G60dB | [S2-perEMS_F2-sweep-change_dir-off](S2-perEMS_F2-sweep-change_dir-off/) | chart | 5 | 25 | 37.5 | 51.6 | 28.7 | 48.1 |
| est_G60dB | [S2-perEMS_F3-off-subtract_dir-on](S2-perEMS_F3-off-subtract_dir-on/) | chart | 1 | 3 | 50.1 | 53.6 | 47.3 | 51.9 |
| est_G60dB | [S2-perEMS_F3-off-subtract_dir-on](S2-perEMS_F3-off-subtract_dir-on/) | chart | 2 | 8 | 44.5 | 49.7 | 36.9 | 50.0 |
| est_G60dB | [S2-perEMS_F3-off-subtract_dir-on](S2-perEMS_F3-off-subtract_dir-on/) | chart | 3 | 15 | 44.1 | 49.3 | 33.5 | 51.2 |
| est_G60dB | [S2-perEMS_F3-off-subtract_dir-on](S2-perEMS_F3-off-subtract_dir-on/) | chart | 4 | 24 | 42.9 | 49.6 | 32.2 | 50.5 |
| est_G60dB | [S2-perEMS_F3-off-subtract_dir-on](S2-perEMS_F3-off-subtract_dir-on/) | chart | 5 | 35 | 42.8 | 50.3 | 32.6 | 51.2 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_orig | 1 | 1 | 5.7 | 13.2 | 6.5 | 5.1 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_free | 1 | 1 | 5.7 | 13.2 | 6.5 | 5.1 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_orig | 2 | 4 | 5.5 | 11.3 | 6.9 | 4.8 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_free | 2 | 4 | 5.5 | 11.9 | 7.0 | 4.9 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_orig | 3 | 9 | 5.5 | 11.3 | 6.7 | 4.7 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_free | 3 | 9 | 5.5 | 11.8 | 6.7 | 4.8 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_orig | 4 | 16 | 5.5 | 11.2 | 7.0 | 4.7 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_free | 4 | 16 | 5.6 | 11.8 | 7.2 | 4.8 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_orig | 5 | 25 | 4.9 | 10.1 | 6.3 | 4.0 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_free | 5 | 25 | 5.0 | 10.9 | 6.4 | 4.0 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-off](S3-joint-allUE_F0-cov-per-sweep_dir-off/) | local_orig | 1 | 1 | 40.3 | 50.3 | 39.6 | 40.9 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-off](S3-joint-allUE_F0-cov-per-sweep_dir-off/) | local_free | 1 | 1 | 40.3 | 50.3 | 39.6 | 40.9 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-off](S3-joint-allUE_F0-cov-per-sweep_dir-off/) | local_orig | 2 | 4 | 32.5 | 36.8 | 31.7 | 33.5 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-off](S3-joint-allUE_F0-cov-per-sweep_dir-off/) | local_free | 2 | 4 | 35.5 | 39.5 | 34.5 | 37.1 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-off](S3-joint-allUE_F0-cov-per-sweep_dir-off/) | local_orig | 3 | 9 | 29.7 | 34.1 | 28.0 | 31.5 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-off](S3-joint-allUE_F0-cov-per-sweep_dir-off/) | local_free | 3 | 9 | 33.1 | 37.1 | 31.6 | 34.8 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-off](S3-joint-allUE_F0-cov-per-sweep_dir-off/) | local_orig | 4 | 16 | 30.0 | 33.6 | 29.4 | 30.7 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-off](S3-joint-allUE_F0-cov-per-sweep_dir-off/) | local_free | 4 | 16 | 33.6 | 36.6 | 32.7 | 34.4 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-off](S3-joint-allUE_F0-cov-per-sweep_dir-off/) | local_orig | 5 | 25 | 29.8 | 33.1 | 29.2 | 30.8 |
| est_G60dB | [S3-joint-allUE_F0-cov-per-sweep_dir-off](S3-joint-allUE_F0-cov-per-sweep_dir-off/) | local_free | 5 | 25 | 33.6 | 36.2 | 32.6 | 34.2 |
| est_G60dB | [S3-joint-allUE_F1-stacked-cov_dir-on](S3-joint-allUE_F1-stacked-cov_dir-on/) | chart | 2 | 4 | 4.0 | 10.9 | 4.5 | 3.5 |
| est_G60dB | [S3-joint-allUE_F1-stacked-cov_dir-on](S3-joint-allUE_F1-stacked-cov_dir-on/) | chart | 3 | 9 | 3.9 | 10.8 | 4.3 | 3.6 |
| est_G60dB | [S3-joint-allUE_F1-stacked-cov_dir-on](S3-joint-allUE_F1-stacked-cov_dir-on/) | chart | 4 | 16 | 3.8 | 10.2 | 4.1 | 3.5 |
| est_G60dB | [S3-joint-allUE_F1-stacked-cov_dir-on](S3-joint-allUE_F1-stacked-cov_dir-on/) | chart | 5 | 25 | 3.8 | 10.7 | 4.1 | 3.6 |
| est_G60dB | [S3-joint-allUE_F1-stacked-cov_dir-off](S3-joint-allUE_F1-stacked-cov_dir-off/) | chart | 2 | 4 | 32.5 | 41.1 | 29.6 | 35.2 |
| est_G60dB | [S3-joint-allUE_F1-stacked-cov_dir-off](S3-joint-allUE_F1-stacked-cov_dir-off/) | chart | 3 | 9 | 28.1 | 37.3 | 27.5 | 28.6 |
| est_G60dB | [S3-joint-allUE_F1-stacked-cov_dir-off](S3-joint-allUE_F1-stacked-cov_dir-off/) | chart | 4 | 16 | 24.5 | 39.4 | 22.1 | 26.8 |
| est_G60dB | [S3-joint-allUE_F1-stacked-cov_dir-off](S3-joint-allUE_F1-stacked-cov_dir-off/) | chart | 5 | 25 | 21.1 | 35.7 | 17.8 | 23.4 |
| est_G60dB | [S3-joint-allUE_F2-sweep-change_dir-on](S3-joint-allUE_F2-sweep-change_dir-on/) | chart | 2 | 4 | 37.5 | 49.6 | 36.3 | 38.8 |
| est_G60dB | [S3-joint-allUE_F2-sweep-change_dir-on](S3-joint-allUE_F2-sweep-change_dir-on/) | chart | 3 | 9 | 35.3 | 48.8 | 32.5 | 38.2 |
| est_G60dB | [S3-joint-allUE_F2-sweep-change_dir-on](S3-joint-allUE_F2-sweep-change_dir-on/) | chart | 4 | 16 | 39.5 | 48.2 | 36.1 | 43.9 |
| est_G60dB | [S3-joint-allUE_F2-sweep-change_dir-on](S3-joint-allUE_F2-sweep-change_dir-on/) | chart | 5 | 25 | 40.7 | 50.7 | 37.9 | 43.8 |
| est_G60dB | [S3-joint-allUE_F2-sweep-change_dir-off](S3-joint-allUE_F2-sweep-change_dir-off/) | chart | 2 | 4 | 35.7 | 44.3 | 31.4 | 40.7 |
| est_G60dB | [S3-joint-allUE_F2-sweep-change_dir-off](S3-joint-allUE_F2-sweep-change_dir-off/) | chart | 3 | 9 | 29.0 | 38.4 | 27.5 | 30.3 |
| est_G60dB | [S3-joint-allUE_F2-sweep-change_dir-off](S3-joint-allUE_F2-sweep-change_dir-off/) | chart | 4 | 16 | 26.9 | 41.4 | 24.3 | 30.9 |
| est_G60dB | [S3-joint-allUE_F2-sweep-change_dir-off](S3-joint-allUE_F2-sweep-change_dir-off/) | chart | 5 | 25 | 22.5 | 37.1 | 20.8 | 24.8 |
| est_G60dB | [S3-joint-allUE_F3-off-subtract_dir-on](S3-joint-allUE_F3-off-subtract_dir-on/) | chart | 1 | 2 | 45.0 | 54.9 | 45.5 | 44.3 |
| est_G60dB | [S3-joint-allUE_F3-off-subtract_dir-on](S3-joint-allUE_F3-off-subtract_dir-on/) | chart | 2 | 5 | 39.2 | 49.8 | 37.5 | 40.5 |
| est_G60dB | [S3-joint-allUE_F3-off-subtract_dir-on](S3-joint-allUE_F3-off-subtract_dir-on/) | chart | 3 | 10 | 37.2 | 48.6 | 36.1 | 37.9 |
| est_G60dB | [S3-joint-allUE_F3-off-subtract_dir-on](S3-joint-allUE_F3-off-subtract_dir-on/) | chart | 4 | 17 | 35.8 | 49.0 | 35.8 | 35.9 |
| est_G60dB | [S3-joint-allUE_F3-off-subtract_dir-on](S3-joint-allUE_F3-off-subtract_dir-on/) | chart | 5 | 26 | 35.9 | 49.0 | 35.6 | 36.4 |
| est_G60dB | [S4-perEMS-frontboth_F0-cov-per-sweep_dir-on](S4-perEMS-frontboth_F0-cov-per-sweep_dir-on/) | local_orig | 3 | 9 | 5.4 | 10.3 | 5.9 | 4.6 |
| est_G60dB | [S4-perEMS-frontboth_F0-cov-per-sweep_dir-on](S4-perEMS-frontboth_F0-cov-per-sweep_dir-on/) | local_free | 3 | 9 | 5.4 | 11.1 | 6.1 | 4.7 |

Each scenario folder has its own README (what changes and which scenarios isolate each factor), `trend.png` (error vs codebook size) and one `K<k>/` folder per codebook size. `comparison_common.png` / `comparison_natural.png` compare all scenarios.
