# AoA ablations, noise model `noiseless_1e6`

Noise: Exact covariance (no estimation noise) with the OLD fixed 1e-6 regularization, exactly as the original code and the paper's position experiment (Main_Position_CodebookSize on main): for like-for-like comparisons with them.

Users: 3200 (in front of both surfaces: 2033; EMS1: 2262; EMS2: 2971). Runtime 11 min.

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
- **Codebook range `narrow`** (folder suffix ``): phase steps 0..40 deg/element (beam steered only within a ~16-28 deg sector away from the BS), as all earlier runs.
- **Codebook range `wide`** (folder suffix `_cb-wide`): phase steps spread evenly over -180..180 deg/element (beams spread evenly over all directions).
- **Noise `noiseless_1e6`**: Exact covariance (no estimation noise) with the OLD fixed 1e-6 regularization, exactly as the original code and the paper's position experiment (Main_Position_CodebookSize on main): for like-for-like comparisons with them.

## Results at a glance (pooled EMS1+EMS2 median |error| [deg], subset `common`)

| noise | scenario | method | K | acquisitions | pooled median | pooled mean | EMS1 median | EMS2 median |
|---|---|---|---|---|---|---|---|---|
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 1 | 1 | 5.0 | 12.5 | 5.1 | 4.9 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 1 | 1 | 5.0 | 12.5 | 5.1 | 4.9 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 2 | 4 | 4.3 | 9.5 | 4.4 | 4.3 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 2 | 4 | 4.4 | 10.4 | 4.4 | 4.4 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 3 | 9 | 4.1 | 9.0 | 4.0 | 4.1 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 3 | 9 | 4.1 | 9.9 | 4.1 | 4.2 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 4 | 16 | 3.9 | 9.0 | 3.7 | 4.0 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 4 | 16 | 4.0 | 9.9 | 3.9 | 4.1 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_orig | 5 | 25 | 3.8 | 8.9 | 3.8 | 3.9 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on](S1-joint_F0-cov-per-sweep_dir-on/) | local_free | 5 | 25 | 3.9 | 9.8 | 3.9 | 3.9 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_orig | 1 | 1 | 3.8 | 11.6 | 3.9 | 3.7 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_free | 1 | 1 | 3.8 | 11.6 | 3.9 | 3.7 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_orig | 2 | 4 | 3.6 | 9.0 | 3.6 | 3.6 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_free | 2 | 4 | 3.6 | 9.9 | 3.6 | 3.6 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_orig | 3 | 9 | 3.4 | 8.6 | 3.4 | 3.5 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_free | 3 | 9 | 3.5 | 9.5 | 3.4 | 3.6 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_orig | 4 | 16 | 3.3 | 8.6 | 3.2 | 3.4 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_free | 4 | 16 | 3.4 | 9.5 | 3.3 | 3.4 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_orig | 5 | 25 | 3.3 | 8.5 | 3.2 | 3.4 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on](S3-joint-allUE_F0-cov-per-sweep_dir-on/) | local_free | 5 | 25 | 3.3 | 9.5 | 3.3 | 3.4 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on_cb-wide](S1-joint_F0-cov-per-sweep_dir-on_cb-wide/) | local_orig | 1 | 1 | 5.0 | 12.5 | 5.1 | 4.9 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on_cb-wide](S1-joint_F0-cov-per-sweep_dir-on_cb-wide/) | local_free | 1 | 1 | 5.0 | 12.5 | 5.1 | 4.9 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on_cb-wide](S1-joint_F0-cov-per-sweep_dir-on_cb-wide/) | local_orig | 2 | 4 | 4.4 | 10.1 | 4.5 | 4.3 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on_cb-wide](S1-joint_F0-cov-per-sweep_dir-on_cb-wide/) | local_free | 2 | 4 | 4.4 | 10.9 | 4.6 | 4.3 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on_cb-wide](S1-joint_F0-cov-per-sweep_dir-on_cb-wide/) | local_orig | 3 | 9 | 4.1 | 9.6 | 4.2 | 4.0 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on_cb-wide](S1-joint_F0-cov-per-sweep_dir-on_cb-wide/) | local_free | 3 | 9 | 4.1 | 10.5 | 4.2 | 4.0 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on_cb-wide](S1-joint_F0-cov-per-sweep_dir-on_cb-wide/) | local_orig | 4 | 16 | 3.9 | 9.4 | 4.0 | 3.9 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on_cb-wide](S1-joint_F0-cov-per-sweep_dir-on_cb-wide/) | local_free | 4 | 16 | 4.0 | 10.4 | 4.1 | 3.9 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on_cb-wide](S1-joint_F0-cov-per-sweep_dir-on_cb-wide/) | local_orig | 5 | 25 | 4.0 | 9.1 | 4.1 | 3.9 |
| noiseless_1e6 | [S1-joint_F0-cov-per-sweep_dir-on_cb-wide](S1-joint_F0-cov-per-sweep_dir-on_cb-wide/) | local_free | 5 | 25 | 4.0 | 10.3 | 4.2 | 3.9 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide](S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide/) | local_orig | 1 | 1 | 3.8 | 11.6 | 3.9 | 3.7 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide](S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide/) | local_free | 1 | 1 | 3.8 | 11.6 | 3.9 | 3.7 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide](S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide/) | local_orig | 2 | 4 | 3.5 | 9.5 | 3.7 | 3.4 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide](S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide/) | local_free | 2 | 4 | 3.6 | 10.3 | 3.7 | 3.4 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide](S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide/) | local_orig | 3 | 9 | 3.4 | 9.1 | 3.4 | 3.4 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide](S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide/) | local_free | 3 | 9 | 3.4 | 10.1 | 3.5 | 3.4 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide](S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide/) | local_orig | 4 | 16 | 3.4 | 9.0 | 3.4 | 3.4 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide](S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide/) | local_free | 4 | 16 | 3.4 | 10.0 | 3.5 | 3.3 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide](S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide/) | local_orig | 5 | 25 | 3.4 | 8.8 | 3.6 | 3.3 |
| noiseless_1e6 | [S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide](S3-joint-allUE_F0-cov-per-sweep_dir-on_cb-wide/) | local_free | 5 | 25 | 3.5 | 10.1 | 3.7 | 3.3 |

Each scenario folder has its own README (what changes and which scenarios isolate each factor), `trend.png` (error vs codebook size) and one `K<k>/` folder per codebook size. `comparison_common.png` / `comparison_natural.png` compare all scenarios.
