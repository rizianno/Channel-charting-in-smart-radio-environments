# AoA from estimated positions vs direct AoA estimation (test)

The position charts of `Main_Position_CodebookSize` (no new t-SNE) are fused per codebook size, the estimated positions are converted into the azimuth seen from EMS1 and EMS2, and scored on 1726 users (in front of both surfaces, anchor in none of the compared runs). Generated 2026-10-22 08:12 by `experiments/Position_to_AoA.m`.

**Fair comparison** (same noiseless setting, 1e-6 regularization, same 3x3 narrow codebook): `position -> AoA, narrow` K=3 vs `direct AoA, narrow (old run, same settings)`. The ablation rows use estimated channels (60 dB pilot gain), so they are a secondary comparison.

| method | K | fusion | EMS1 RMSE | EMS2 RMSE | joint RMSE | EMS1 median | EMS2 median | joint median |
|---|---|---|---|---|---|---|---|---|
| position -> AoA, narrow | 1x1 | local_fusion | 37.2 | 29.7 | 47.6 | 4.7 | 4.8 | 8.8 |
| position -> AoA, narrow | 1x1 | global_fusion | 37.2 | 29.7 | 47.6 | 4.7 | 4.8 | 8.8 |
| position -> AoA, narrow | 1x1 | typical_single | 37.2 | 29.7 | 47.6 | 4.7 | 4.8 | 8.8 |
| position -> AoA, narrow | 2x2 | local_fusion | 32.6 | 24.8 | 41.0 | 4.1 | 4.3 | 7.4 |
| position -> AoA, narrow | 2x2 | global_fusion | 35.3 | 26.6 | 44.1 | 4.2 | 4.3 | 7.5 |
| position -> AoA, narrow | 2x2 | typical_single | 35.8 | 29.0 | 46.1 | 4.7 | 4.9 | 8.9 |
| position -> AoA, narrow | 3x3 | local_fusion | 29.8 | 23.5 | 37.9 | 3.8 | 4.1 | 7.0 |
| position -> AoA, narrow | 3x3 | global_fusion | 32.8 | 25.5 | 41.6 | 3.9 | 4.1 | 7.3 |
| position -> AoA, narrow | 3x3 | typical_single | 34.9 | 28.9 | 45.3 | 4.6 | 5.0 | 8.8 |
| position -> AoA, narrow | 4x4 | local_fusion | 31.5 | 23.6 | 39.4 | 3.7 | 4.0 | 6.9 |
| position -> AoA, narrow | 4x4 | global_fusion | 34.4 | 25.5 | 42.8 | 3.9 | 4.2 | 7.3 |
| position -> AoA, narrow | 4x4 | typical_single | 36.1 | 29.2 | 46.5 | 4.6 | 4.9 | 8.6 |
| position -> AoA, narrow | 5x5 | local_fusion | 29.1 | 22.8 | 37.0 | 3.7 | 3.9 | 6.8 |
| position -> AoA, narrow | 5x5 | global_fusion | 32.8 | 25.1 | 41.3 | 3.7 | 4.0 | 7.0 |
| position -> AoA, narrow | 5x5 | typical_single | 34.7 | 29.0 | 45.2 | 4.6 | 4.9 | 8.7 |
| position -> AoA, wide | 1x1 | local_fusion | 37.2 | 29.7 | 47.6 | 4.7 | 4.8 | 8.8 |
| position -> AoA, wide | 1x1 | global_fusion | 37.2 | 29.7 | 47.6 | 4.7 | 4.8 | 8.8 |
| position -> AoA, wide | 1x1 | typical_single | 37.2 | 29.7 | 47.6 | 4.7 | 4.8 | 8.8 |
| position -> AoA, wide | 2x2 | local_fusion | 30.4 | 24.0 | 38.7 | 4.3 | 4.2 | 7.8 |
| position -> AoA, wide | 2x2 | global_fusion | 34.1 | 25.9 | 42.8 | 4.4 | 4.2 | 7.9 |
| position -> AoA, wide | 2x2 | typical_single | 35.4 | 29.4 | 46.0 | 4.8 | 4.9 | 9.0 |
| position -> AoA, wide | 3x3 | local_fusion | 28.5 | 24.1 | 37.3 | 4.1 | 4.2 | 7.6 |
| position -> AoA, wide | 3x3 | global_fusion | 33.2 | 25.7 | 42.0 | 4.1 | 4.2 | 7.7 |
| position -> AoA, wide | 3x3 | typical_single | 34.8 | 28.8 | 45.2 | 4.9 | 4.9 | 9.0 |
| position -> AoA, wide | 4x4 | local_fusion | 28.1 | 22.6 | 36.1 | 4.1 | 4.1 | 7.3 |
| position -> AoA, wide | 4x4 | global_fusion | 33.4 | 24.9 | 41.7 | 4.1 | 4.1 | 7.4 |
| position -> AoA, wide | 4x4 | typical_single | 35.4 | 29.0 | 45.8 | 4.8 | 4.9 | 9.0 |
| position -> AoA, wide | 5x5 | local_fusion | 27.1 | 21.8 | 34.8 | 3.9 | 4.1 | 7.4 |
| position -> AoA, wide | 5x5 | global_fusion | 32.6 | 24.9 | 41.1 | 4.0 | 4.0 | 7.5 |
| position -> AoA, wide | 5x5 | typical_single | 35.1 | 28.6 | 45.3 | 4.8 | 4.9 | 8.9 |
| direct AoA, narrow (old run, same settings) | 3x3 | local_fusion | 22.1 | 18.0 | 28.5 | 3.8 | 4.1 | 6.7 |
| direct AoA, narrow (old run, same settings) | 3x3 | global_fusion | 23.5 | 20.6 | 31.2 | 3.9 | 4.2 | 6.8 |
| direct AoA (est. channels 60 dB), users in front of both, narrow | 1x1 | local_fusion | 26.4 | 27.0 | 37.8 | 6.4 | 5.1 | 10.4 |
| direct AoA (est. channels 60 dB), users in front of both, narrow | 2x2 | local_fusion | 23.6 | 23.4 | 33.2 | 6.2 | 4.9 | 9.9 |
| direct AoA (est. channels 60 dB), users in front of both, narrow | 3x3 | local_fusion | 23.4 | 22.7 | 32.6 | 5.6 | 4.8 | 9.6 |
| direct AoA (est. channels 60 dB), users in front of both, narrow | 4x4 | local_fusion | 23.5 | 22.7 | 32.7 | 5.6 | 4.6 | 9.4 |
| direct AoA (est. channels 60 dB), users in front of both, narrow | 5x5 | local_fusion | 23.1 | 21.5 | 31.5 | 5.4 | 4.7 | 9.3 |
| direct AoA (est. channels 60 dB), users in front of both, wide | 1x1 | local_fusion | 26.4 | 27.0 | 37.8 | 6.4 | 5.1 | 10.4 |
| direct AoA (est. channels 60 dB), users in front of both, wide | 2x2 | local_fusion | 24.6 | 24.2 | 34.5 | 6.1 | 4.8 | 10.1 |
| direct AoA (est. channels 60 dB), users in front of both, wide | 3x3 | local_fusion | 24.4 | 24.0 | 34.2 | 6.1 | 4.9 | 10.5 |
| direct AoA (est. channels 60 dB), users in front of both, wide | 4x4 | local_fusion | 24.3 | 23.5 | 33.8 | 6.4 | 4.8 | 10.3 |
| direct AoA (est. channels 60 dB), users in front of both, wide | 5x5 | local_fusion | 23.8 | 21.6 | 32.2 | 6.1 | 5.0 | 10.3 |
| direct AoA (est. channels 60 dB), all users in t-SNE, narrow | 1x1 | local_fusion | 25.6 | 27.7 | 37.8 | 6.5 | 5.1 | 10.7 |
| direct AoA (est. channels 60 dB), all users in t-SNE, narrow | 2x2 | local_fusion | 22.5 | 20.0 | 30.1 | 6.9 | 4.8 | 10.7 |
| direct AoA (est. channels 60 dB), all users in t-SNE, narrow | 3x3 | local_fusion | 22.4 | 21.7 | 31.2 | 6.7 | 4.7 | 10.5 |
| direct AoA (est. channels 60 dB), all users in t-SNE, narrow | 4x4 | local_fusion | 22.8 | 19.8 | 30.1 | 7.0 | 4.7 | 10.6 |
| direct AoA (est. channels 60 dB), all users in t-SNE, narrow | 5x5 | local_fusion | 21.2 | 18.7 | 28.3 | 6.3 | 4.0 | 9.6 |
| direct AoA (est. channels 60 dB), all users in t-SNE, wide | 1x1 | local_fusion | 25.6 | 27.7 | 37.8 | 6.5 | 5.1 | 10.7 |
| direct AoA (est. channels 60 dB), all users in t-SNE, wide | 2x2 | local_fusion | 22.8 | 20.8 | 30.9 | 6.6 | 4.2 | 10.0 |
| direct AoA (est. channels 60 dB), all users in t-SNE, wide | 3x3 | local_fusion | 23.2 | 24.3 | 33.6 | 6.5 | 4.9 | 10.5 |
| direct AoA (est. channels 60 dB), all users in t-SNE, wide | 4x4 | local_fusion | 22.4 | 21.6 | 31.1 | 6.1 | 4.4 | 9.8 |
| direct AoA (est. channels 60 dB), all users in t-SNE, wide | 5x5 | local_fusion | 21.0 | 17.3 | 27.2 | 5.4 | 3.8 | 8.7 |

Angles in degrees. Joint = sqrt(e1^2 + e2^2) per user. Position run: all 3200 users, 480 anchors (first users in dataset order); direct AoA runs: users in front of both surfaces (or all users), 305 anchors.
