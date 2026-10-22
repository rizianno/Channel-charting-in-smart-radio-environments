# AoA from estimated positions vs direct AoA estimation (test)

The position charts of `Main_Position_CodebookSize` (no new t-SNE) are fused per codebook size, the estimated positions are converted into the azimuth seen from EMS1 and EMS2 ("position proxy"), and compared with direct AoA charts (azimuths as t-SNE labels, per-sweep charts + local fusion, direct path kept), on 1726 users (in front of both surfaces, anchor in none of the runs). Generated 2026-10-22 08:32 by `experiments/Position_to_AoA.m`.

- `proxy_vs_direct.png`: **fair comparison**, all noiseless with the 1e-6 regularization (the position run's setting).
- `direct_noiseless_vs_noisy.png`: direct AoA, noiseless vs estimated channels (60 dB pilot gain).
- Rows of each figure: joint error sqrt(e1^2+e2^2) and per-surface error (EMS1 and EMS2 pooled); columns: median, mean, 90th percentile, RMSE; solid = narrow phase steps (0..40 deg/element), dashed = wide (-180..180).

## Joint error, local fusion [deg]: median / mean / p90 / RMSE

| source | users | range | 1x1 | 2x2 | 3x3 | 4x4 | 5x5 |
|---|---|---|---|---|---|---|---|
| position proxy (noiseless, 1e-6) | all users | narrow | 8.8 / 24.5 / 76.1 / 47.6 | 7.4 / 20.7 / 54.7 / 41.0 | 7.0 / 19.2 / 48.6 / 37.9 | 6.9 / 19.6 / 53.0 / 39.4 | 6.8 / 18.7 / 47.7 / 37.0 |
| position proxy (noiseless, 1e-6) | all users | wide | 8.8 / 24.5 / 76.1 / 47.6 | 7.8 / 20.3 / 58.4 / 38.7 | 7.6 / 19.5 / 54.2 / 37.3 | 7.3 / 19.1 / 53.9 / 36.1 | 7.4 / 18.5 / 49.6 / 34.8 |
| direct AoA (noiseless, 1e-6) | users in front of both | narrow | 8.8 / 19.7 / 42.2 / 39.3 | 7.4 / 15.1 / 29.9 / 29.2 | 6.9 / 14.4 / 28.3 / 28.1 | 6.6 / 14.4 / 29.6 / 28.5 | 6.5 / 14.2 / 28.7 / 28.0 |
| direct AoA (noiseless, 1e-6) | users in front of both | wide | 8.8 / 19.7 / 42.2 / 39.3 | 7.3 / 16.0 / 34.4 / 31.3 | 7.0 / 15.3 / 31.7 / 30.3 | 6.8 / 15.0 / 31.2 / 29.3 | 6.7 / 14.5 / 31.4 / 27.8 |
| direct AoA (noiseless, 1e-6) | all users | narrow | 6.9 / 18.3 / 39.9 / 38.6 | 6.2 / 14.3 / 29.3 / 28.9 | 5.9 / 13.7 / 28.6 / 27.8 | 5.7 / 13.7 / 28.4 / 28.1 | 5.7 / 13.6 / 28.6 / 27.6 |
| direct AoA (noiseless, 1e-6) | all users | wide | 6.9 / 18.3 / 39.9 / 38.6 | 6.4 / 15.1 / 33.0 / 30.8 | 6.1 / 14.5 / 31.4 / 29.8 | 5.9 / 14.3 / 31.6 / 28.9 | 5.9 / 14.0 / 30.7 / 27.8 |
| direct AoA (est. channels, 60 dB) | users in front of both | narrow | 10.4 / 20.6 / 43.4 / 37.8 | 9.9 / 18.7 / 39.3 / 33.2 | 9.6 / 18.1 / 39.1 / 32.6 | 9.4 / 18.1 / 39.6 / 32.7 | 9.3 / 17.8 / 39.6 / 31.5 |
| direct AoA (est. channels, 60 dB) | users in front of both | wide | 10.4 / 20.6 / 43.4 / 37.8 | 10.1 / 19.2 / 40.7 / 34.5 | 10.5 / 19.1 / 40.3 / 34.2 | 10.3 / 19.1 / 41.1 / 33.8 | 10.3 / 18.5 / 40.4 / 32.2 |
| direct AoA (est. channels, 60 dB) | all users | narrow | 10.7 / 21.4 / 46.9 / 37.8 | 10.7 / 18.1 / 39.2 / 30.1 | 10.5 / 18.2 / 38.8 / 31.2 | 10.6 / 17.9 / 38.3 / 30.1 | 9.6 / 16.2 / 33.7 / 28.3 |
| direct AoA (est. channels, 60 dB) | all users | wide | 10.7 / 21.4 / 46.9 / 37.8 | 10.0 / 17.7 / 35.8 / 30.9 | 10.5 / 19.1 / 39.9 / 33.6 | 9.8 / 17.8 / 35.4 / 31.1 | 8.7 / 15.3 / 30.5 / 27.2 |

All statistics (also global fusion and the typical single configuration of the position proxy, joint and pooled errors) are in `summary.csv`.
