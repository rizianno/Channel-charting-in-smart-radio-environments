# Full-NLOS users: EMS beam selection vs number of sweeps, surfaces 60..240 (noise `est_G60dB`)

Generated 2026-10-23 10:54 by `experiments/Main_AoA_Serving_NLoS_Size.m` (branch aoa-serving), runtime 6 min.

## Setup
- **Users:** full NLOS = no BS-UE channel at all (885 users); the EMSs are their only link.
- **t-SNE:** one 1-D chart per EMS with the users of this population it reaches (EMS1: 419 users, 63 anchors, perplexity 84; EMS2: 737, 111, 147). Scored: reachable by at least one EMS and anchor in neither chart.
- **Surfaces:** M x M elements, lambda/2, steering along one axis only (vertical: same phase on every column). `M60_raytraced` = the dataset; `M60/120/180/240` = synthesized with spherical waves from the UE, EMS and BS positions, per-user strengths fitted on the 60 x 60 data. Codebook from the beamwidth: N = 68/136/204/272.
- **Schemes vs number of sweeps:** sequential, spread, hierarchical (4 spread entries, then halving), spread + R neighbours, t-SNE AoA + LUT (K = [1 2 3 4 5 6 7 8 9 10] spread training entries; sweep k sets both EMSs to training entry k, K training sweeps), t-SNE + LUT entry + R neighbours (R = [2 4 6]), exhaustive (N), true AoA + LUT.
- Measurement-based selection uses each EMS path's own gain (noise-free, in favour of the sweeping baselines); the served SNR uses the true channel; 23 dBm, -92 dBm.

## Validation of the synthesized channels (60 x 60 vs ray-traced)

- EMS1: element grid orientation (m, n) = (+1, +1), phase sign +1; |correlation| of the EMS->BS channel with the model 0.998; UE->EMS channels: median |correlation| 0.999, above 0.9 for 100% of the 419 reachable full-NLOS users.
- EMS2: element grid orientation (m, n) = (+1, -1), phase sign +1; |correlation| of the EMS->BS channel with the model 0.999; UE->EMS channels: median |correlation| 0.999, above 0.9 for 100% of the 737 reachable full-NLOS users.
- validation EMS1 (419 users): best entry same as ray-traced for 100%, within 1 entry for 100%; ray-traced gain at the synthetic best entry: median 0.0 dB below the ray-traced best; best gain synthetic vs ray-traced: median -0.1 dB
- validation EMS2 (737 users): best entry same as ray-traced for 99%, within 1 entry for 100%; ray-traced gain at the synthetic best entry: median 0.0 dB below the ray-traced best; best gain synthetic vs ray-traced: median -0.0 dB

## Median SNR gain vs zero phase [dB] / 10th percentile [dB], full-NLOS users

**M60_raytraced** (codebook 68, 736 users)
| scheme | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -6.0 / -16.3 | -4.3 / -13.5 | -3.2 / -12.3 | -2.6 / -11.4 | -2.0 / -10.5 | -2.0 / -10.0 | -2.0 / -9.9 | -1.9 / -9.7 | -1.5 / -9.3 | -0.6 / -8.4 |
| spread | 0.0 / 0.0 | 2.2 / -6.0 | 0.3 / 0.0 | 15.7 / 3.4 | 12.9 / 4.3 | 14.2 / -2.7 | 16.4 / 4.6 | 19.9 / 9.7 | 22.3 / 8.8 | 24.6 / 13.0 |
| hierarchical | 0.0 / 0.0 | 2.2 / -6.0 | 0.3 / 0.0 | 15.7 / 3.4 | 15.9 / 4.3 | 18.5 / 9.8 | 23.8 / 12.8 | 28.0 / 17.1 | 30.7 / 23.0 | 30.7 / 23.0 |
| spread+refine (R=6) | - | - | - | - | - | 9.3 / -3.3 | 28.7 / 13.8 | 23.9 / 9.8 | 29.7 / 20.6 | 30.2 / 19.8 |
| tsne | - | - | - | - | - | - | - | - | - | - |
| tsne+refine (R=2) | - | - | - | - | - | - | - | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | - | - | - | - | - |
| true AoA + LUT (0 sweeps) | 32.1 / 24.0 | | | | | | | | | |
| exhaustive (68 sweeps) | 32.2 / 24.2 | | | | | | | | | |
Zero-phase median SNR -69.6 dB, exhaustive median SNR -36.2 dB.
t-SNE angle error / LUT: K=1 avg: AoA median |err| 5.3/47.9 deg, RMSE 18.6/49.9; K=1 sharp: AoA median |err| 5.3/47.9 deg, RMSE 18.6/49.9; K=1 select: AoA median |err| 5.3/47.9 deg, RMSE 18.6/49.9; K=2 avg: AoA median |err| 5.9/47.1 deg, RMSE 19.5/49.8; K=2 sharp: AoA median |err| 5.5/47.0 deg, RMSE 18.7/50.0; K=2 select: AoA median |err| 5.5/47.0 deg, RMSE 18.7/50.0; K=3 avg: AoA median |err| 6.1/42.3 deg, RMSE 18.6/45.9; K=3 sharp: AoA median |err| 5.8/38.3 deg, RMSE 18.2/43.8; K=3 select: AoA median |err| 5.3/40.0 deg, RMSE 18.1/43.7; K=4 avg: AoA median |err| 5.9/44.6 deg, RMSE 18.0/45.8; K=4 sharp: AoA median |err| 4.7/37.4 deg, RMSE 17.4/43.8; K=4 select: AoA median |err| 4.1/37.1 deg, RMSE 16.9/43.0; K=5 avg: AoA median |err| 6.2/40.0 deg, RMSE 17.8/42.4; K=5 sharp: AoA median |err| 4.8/33.3 deg, RMSE 17.4/39.4; K=5 select: AoA median |err| 4.2/5.8 deg, RMSE 16.9/36.7; K=6 avg: AoA median |err| 5.7/37.9 deg, RMSE 17.5/40.4; K=6 sharp: AoA median |err| 5.3/28.3 deg, RMSE 17.1/35.8; K=6 select: AoA median |err| 4.8/3.7 deg, RMSE 16.9/29.8; K=7 avg: AoA median |err| 5.9/34.2 deg, RMSE 17.3/37.9; K=7 sharp: AoA median |err| 5.0/25.4 deg, RMSE 16.9/33.3; K=7 select: AoA median |err| 3.9/3.8 deg, RMSE 16.7/30.3; K=8 avg: AoA median |err| 5.9/36.8 deg, RMSE 17.3/39.0; K=8 sharp: AoA median |err| 5.2/28.2 deg, RMSE 16.9/34.6; K=8 select: AoA median |err| 3.9/3.8 deg, RMSE 16.7/30.0; K=9 avg: AoA median |err| 5.9/38.4 deg, RMSE 17.3/39.3; K=9 sharp: AoA median |err| 5.2/30.0 deg, RMSE 17.0/34.9; K=9 select: AoA median |err| 3.9/3.7 deg, RMSE 16.7/29.4; K=10 avg: AoA median |err| 5.8/36.9 deg, RMSE 17.3/38.0; K=10 sharp: AoA median |err| 5.2/28.5 deg, RMSE 17.0/33.5; K=10 select: AoA median |err| 3.9/2.9 deg, RMSE 16.7/28.2

**M60** (codebook 68, 736 users)
| scheme | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -6.0 / -16.1 | -4.3 / -13.5 | -3.2 / -12.3 | -2.6 / -11.4 | -2.0 / -10.5 | -1.9 / -9.9 | -1.9 / -9.9 | -1.8 / -9.7 | -1.4 / -9.3 | -0.5 / -8.3 |
| spread | 0.0 / 0.0 | 2.2 / -6.1 | 0.4 / 0.0 | 16.0 / 4.6 | 13.0 / 4.3 | 14.4 / -1.9 | 16.4 / 5.0 | 20.1 / 10.1 | 22.5 / 8.9 | 24.6 / 13.0 |
| hierarchical | 0.0 / 0.0 | 2.2 / -6.1 | 0.4 / 0.0 | 16.0 / 4.6 | 16.1 / 4.7 | 18.7 / 10.2 | 23.9 / 12.9 | 28.0 / 17.5 | 30.7 / 23.2 | 30.7 / 23.2 |
| spread+refine (R=6) | - | - | - | - | - | 9.3 / -3.3 | 28.7 / 14.3 | 23.9 / 10.1 | 29.8 / 21.0 | 30.2 / 20.2 |
| tsne | - | - | - | - | - | - | - | - | - | - |
| tsne+refine (R=2) | - | - | - | - | - | - | - | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | - | - | - | - | - |
| true AoA + LUT (0 sweeps) | 32.1 / 24.0 | | | | | | | | | |
| exhaustive (68 sweeps) | 32.2 / 24.2 | | | | | | | | | |
Zero-phase median SNR -69.6 dB, exhaustive median SNR -36.2 dB.
t-SNE angle error / LUT: K=1 avg: AoA median |err| 5.2/46.3 deg, RMSE 17.8/48.0; K=1 sharp: AoA median |err| 5.2/46.3 deg, RMSE 17.8/48.0; K=1 select: AoA median |err| 5.2/46.3 deg, RMSE 17.8/48.0; K=2 avg: AoA median |err| 5.7/44.6 deg, RMSE 17.5/47.7; K=2 sharp: AoA median |err| 5.4/47.3 deg, RMSE 17.7/47.2; K=2 select: AoA median |err| 5.4/47.3 deg, RMSE 17.7/47.2; K=3 avg: AoA median |err| 6.0/39.8 deg, RMSE 17.5/43.7; K=3 sharp: AoA median |err| 5.9/25.6 deg, RMSE 17.6/41.5; K=3 select: AoA median |err| 5.9/7.9 deg, RMSE 17.8/41.4; K=4 avg: AoA median |err| 5.6/44.0 deg, RMSE 17.3/44.6; K=4 sharp: AoA median |err| 5.1/32.4 deg, RMSE 17.1/42.1; K=4 select: AoA median |err| 5.0/7.3 deg, RMSE 17.0/41.0; K=5 avg: AoA median |err| 5.8/37.1 deg, RMSE 17.3/41.5; K=5 sharp: AoA median |err| 5.2/24.1 deg, RMSE 17.2/38.1; K=5 select: AoA median |err| 4.8/4.9 deg, RMSE 17.2/35.1; K=6 avg: AoA median |err| 5.8/37.8 deg, RMSE 17.2/39.7; K=6 sharp: AoA median |err| 5.1/25.9 deg, RMSE 17.0/34.9; K=6 select: AoA median |err| 4.8/3.2 deg, RMSE 17.0/28.4; K=7 avg: AoA median |err| 6.1/35.1 deg, RMSE 17.2/38.8; K=7 sharp: AoA median |err| 5.6/24.9 deg, RMSE 17.0/34.3; K=7 select: AoA median |err| 5.2/3.0 deg, RMSE 17.0/27.5; K=8 avg: AoA median |err| 5.8/36.5 deg, RMSE 17.3/40.0; K=8 sharp: AoA median |err| 5.4/27.8 deg, RMSE 17.1/35.7; K=8 select: AoA median |err| 5.2/3.0 deg, RMSE 17.0/27.5; K=9 avg: AoA median |err| 6.0/37.8 deg, RMSE 17.3/40.8; K=9 sharp: AoA median |err| 5.7/29.9 deg, RMSE 17.2/36.7; K=9 select: AoA median |err| 5.2/3.1 deg, RMSE 17.0/27.5; K=10 avg: AoA median |err| 5.8/38.7 deg, RMSE 17.5/39.3; K=10 sharp: AoA median |err| 5.7/31.2 deg, RMSE 17.2/34.9; K=10 select: AoA median |err| 5.2/2.3 deg, RMSE 17.0/22.0

**M120** (codebook 136, 736 users)
| scheme | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -6.2 / -15.5 | -4.7 / -13.5 | -3.5 / -12.0 | -2.7 / -11.0 | -2.0 / -10.0 | -1.9 / -9.6 | -1.9 / -9.5 | -1.9 / -9.5 | -1.7 / -9.5 | -1.6 / -9.5 |
| spread | 0.0 / 0.0 | 2.4 / -5.7 | 0.3 / 0.0 | 15.2 / 6.6 | 11.9 / 4.0 | 14.0 / 6.2 | 17.4 / 7.0 | 19.1 / 10.6 | 22.1 / 7.0 | 23.2 / 14.3 |
| hierarchical | 0.0 / 0.0 | 2.4 / -5.7 | 0.3 / 0.0 | 15.2 / 6.6 | 15.6 / 7.3 | 20.3 / 12.1 | 24.2 / 15.3 | 29.2 / 19.9 | 36.6 / 30.1 | 36.6 / 30.1 |
| spread+refine (R=6) | - | - | - | - | - | 8.3 / -3.1 | 29.8 / 11.1 | 19.7 / 8.7 | 28.1 / 10.3 | 34.4 / 19.5 |
| tsne | - | - | - | - | - | - | - | - | - | - |
| tsne+refine (R=2) | - | - | - | - | - | - | - | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | - | - | - | - | - |
| true AoA + LUT (0 sweeps) | 36.9 / 30.2 | | | | | | | | | |
| exhaustive (136 sweeps) | 37.6 / 31.0 | | | | | | | | | |
Zero-phase median SNR -66.3 dB, exhaustive median SNR -26.6 dB.
t-SNE angle error / LUT: K=1 avg: AoA median |err| 5.2/46.8 deg, RMSE 17.8/49.3; K=1 sharp: AoA median |err| 5.2/46.8 deg, RMSE 17.8/49.3; K=1 select: AoA median |err| 5.2/46.8 deg, RMSE 17.8/49.3; K=2 avg: AoA median |err| 5.4/46.7 deg, RMSE 17.8/48.8; K=2 sharp: AoA median |err| 5.3/46.8 deg, RMSE 17.2/48.8; K=2 select: AoA median |err| 5.3/46.8 deg, RMSE 17.2/48.8; K=3 avg: AoA median |err| 5.8/41.1 deg, RMSE 17.1/44.3; K=3 sharp: AoA median |err| 5.4/37.4 deg, RMSE 17.2/42.5; K=3 select: AoA median |err| 5.3/36.8 deg, RMSE 17.3/41.9; K=4 avg: AoA median |err| 5.5/41.7 deg, RMSE 17.1/43.7; K=4 sharp: AoA median |err| 4.5/31.0 deg, RMSE 17.1/41.2; K=4 select: AoA median |err| 4.1/9.2 deg, RMSE 17.0/40.5; K=5 avg: AoA median |err| 5.7/36.8 deg, RMSE 17.8/40.5; K=5 sharp: AoA median |err| 4.8/25.3 deg, RMSE 17.1/36.9; K=5 select: AoA median |err| 4.1/4.8 deg, RMSE 17.0/34.8; K=6 avg: AoA median |err| 5.3/38.3 deg, RMSE 17.1/40.0; K=6 sharp: AoA median |err| 4.9/27.1 deg, RMSE 16.2/35.8; K=6 select: AoA median |err| 4.5/3.6 deg, RMSE 17.1/32.4; K=7 avg: AoA median |err| 5.0/32.3 deg, RMSE 16.7/36.4; K=7 sharp: AoA median |err| 4.4/21.1 deg, RMSE 15.9/32.1; K=7 select: AoA median |err| 3.9/3.1 deg, RMSE 16.0/32.5; K=8 avg: AoA median |err| 5.3/34.8 deg, RMSE 16.8/37.6; K=8 sharp: AoA median |err| 4.5/24.5 deg, RMSE 16.1/33.3; K=8 select: AoA median |err| 3.9/3.1 deg, RMSE 16.0/32.2; K=9 avg: AoA median |err| 5.4/33.9 deg, RMSE 16.8/37.6; K=9 sharp: AoA median |err| 4.6/24.0 deg, RMSE 16.1/33.5; K=9 select: AoA median |err| 3.9/3.1 deg, RMSE 16.0/31.6; K=10 avg: AoA median |err| 5.4/32.5 deg, RMSE 16.9/36.6; K=10 sharp: AoA median |err| 4.8/23.4 deg, RMSE 16.3/32.6; K=10 select: AoA median |err| 3.9/2.8 deg, RMSE 16.0/30.4

**M180** (codebook 204, 736 users)
| scheme | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -5.9 / -16.0 | -4.4 / -13.8 | -3.4 / -12.2 | -2.7 / -11.1 | -2.1 / -9.9 | -2.0 / -9.7 | -2.0 / -9.7 | -2.0 / -9.7 | -2.0 / -9.7 | -2.0 / -9.7 |
| spread | 0.0 / 0.0 | 2.4 / -6.0 | 0.0 / 0.0 | 15.8 / 5.9 | 11.8 / 6.9 | 15.8 / -0.7 | 16.6 / 5.9 | 19.5 / 9.3 | 22.0 / 12.2 | 24.9 / 14.7 |
| hierarchical | 0.0 / 0.0 | 2.4 / -6.0 | 0.0 / 0.0 | 15.8 / 5.9 | 16.2 / 6.5 | 21.7 / 11.2 | 27.0 / 15.5 | 32.2 / 21.1 | 39.3 / 31.8 | 39.3 / 31.8 |
| spread+refine (R=6) | - | - | - | - | - | 7.7 / -3.2 | 25.7 / 10.8 | 18.4 / 9.2 | 26.2 / 6.6 | 34.0 / 18.5 |
| tsne | - | - | - | - | - | - | - | - | - | - |
| tsne+refine (R=2) | - | - | - | - | - | - | - | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | - | - | - | - | - |
| true AoA + LUT (0 sweeps) | 40.0 / 32.8 | | | | | | | | | |
| exhaustive (204 sweeps) | 40.8 / 34.0 | | | | | | | | | |
Zero-phase median SNR -64.5 dB, exhaustive median SNR -21.9 dB.
t-SNE angle error / LUT: K=1 avg: AoA median |err| 5.1/45.4 deg, RMSE 17.2/49.1; K=1 sharp: AoA median |err| 5.1/45.4 deg, RMSE 17.2/49.1; K=1 select: AoA median |err| 5.1/45.4 deg, RMSE 17.2/49.1; K=2 avg: AoA median |err| 5.1/41.7 deg, RMSE 17.6/46.9; K=2 sharp: AoA median |err| 5.0/42.6 deg, RMSE 17.1/45.9; K=2 select: AoA median |err| 5.0/42.6 deg, RMSE 17.1/45.9; K=3 avg: AoA median |err| 5.2/38.1 deg, RMSE 19.4/43.1; K=3 sharp: AoA median |err| 5.0/31.8 deg, RMSE 19.2/40.7; K=3 select: AoA median |err| 5.1/16.2 deg, RMSE 16.7/39.6; K=4 avg: AoA median |err| 5.1/35.7 deg, RMSE 18.5/41.3; K=4 sharp: AoA median |err| 4.3/24.6 deg, RMSE 17.7/37.6; K=4 select: AoA median |err| 3.7/6.1 deg, RMSE 17.3/36.4; K=5 avg: AoA median |err| 5.4/33.2 deg, RMSE 18.0/38.2; K=5 sharp: AoA median |err| 4.5/23.3 deg, RMSE 17.4/33.8; K=5 select: AoA median |err| 3.7/4.0 deg, RMSE 17.3/30.9; K=6 avg: AoA median |err| 5.0/34.2 deg, RMSE 16.0/36.8; K=6 sharp: AoA median |err| 4.5/23.6 deg, RMSE 16.0/31.6; K=6 select: AoA median |err| 4.4/3.1 deg, RMSE 17.3/26.0; K=7 avg: AoA median |err| 5.2/29.6 deg, RMSE 15.9/34.3; K=7 sharp: AoA median |err| 4.6/19.9 deg, RMSE 16.0/29.6; K=7 select: AoA median |err| 3.3/3.3 deg, RMSE 16.6/27.9; K=8 avg: AoA median |err| 5.4/31.8 deg, RMSE 16.1/35.6; K=8 sharp: AoA median |err| 4.8/22.5 deg, RMSE 16.2/30.9; K=8 select: AoA median |err| 3.3/3.3 deg, RMSE 16.6/27.9; K=9 avg: AoA median |err| 5.4/32.5 deg, RMSE 16.3/35.3; K=9 sharp: AoA median |err| 4.8/23.4 deg, RMSE 16.3/30.6; K=9 select: AoA median |err| 3.3/3.3 deg, RMSE 16.6/27.0; K=10 avg: AoA median |err| 5.3/31.4 deg, RMSE 16.4/34.4; K=10 sharp: AoA median |err| 4.8/22.4 deg, RMSE 16.4/30.0; K=10 select: AoA median |err| 3.3/2.6 deg, RMSE 16.6/26.1

**M240** (codebook 272, 736 users)
| scheme | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -6.2 / -14.9 | -4.7 / -13.2 | -3.5 / -12.1 | -2.8 / -11.2 | -1.7 / -9.9 | -1.7 / -9.5 | -1.7 / -9.5 | -1.7 / -9.5 | -1.7 / -9.5 | -1.6 / -9.5 |
| spread | 0.0 / 0.0 | 2.2 / -6.0 | 0.6 / 0.0 | 16.2 / 6.7 | 13.1 / 5.2 | 15.0 / 7.0 | 18.1 / 5.8 | 20.3 / 10.8 | 22.5 / 13.2 | 25.5 / 15.3 |
| hierarchical | 0.0 / 0.0 | 2.2 / -6.0 | 0.6 / 0.0 | 16.2 / 6.7 | 16.6 / 8.2 | 22.4 / 13.7 | 29.1 / 18.9 | 34.4 / 24.8 | 39.5 / 31.1 | 42.2 / 35.4 |
| spread+refine (R=6) | - | - | - | - | - | 7.8 / -3.4 | 25.2 / 9.9 | 17.6 / 9.0 | 26.2 / 9.8 | 33.3 / 15.6 |
| tsne | - | - | - | - | - | - | - | - | - | - |
| tsne+refine (R=2) | - | - | - | - | - | - | - | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | - | - | - | - | - |
| true AoA + LUT (0 sweeps) | 41.8 / 34.7 | | | | | | | | | |
| exhaustive (272 sweeps) | 42.6 / 35.6 | | | | | | | | | |
Zero-phase median SNR -63.2 dB, exhaustive median SNR -19.4 dB.
t-SNE angle error / LUT: K=1 avg: AoA median |err| 5.2/6.9 deg, RMSE 18.2/54.4; K=1 sharp: AoA median |err| 5.2/6.9 deg, RMSE 18.2/54.4; K=1 select: AoA median |err| 5.2/6.9 deg, RMSE 18.2/54.4; K=2 avg: AoA median |err| 5.6/25.9 deg, RMSE 17.8/42.1; K=2 sharp: AoA median |err| 5.1/6.1 deg, RMSE 17.1/49.7; K=2 select: AoA median |err| 5.1/6.1 deg, RMSE 17.1/49.7; K=3 avg: AoA median |err| 5.5/22.5 deg, RMSE 17.3/38.7; K=3 sharp: AoA median |err| 5.2/12.7 deg, RMSE 17.3/40.5; K=3 select: AoA median |err| 5.1/5.3 deg, RMSE 17.4/45.7; K=4 avg: AoA median |err| 5.1/21.4 deg, RMSE 17.4/37.7; K=4 sharp: AoA median |err| 4.2/16.3 deg, RMSE 17.1/36.7; K=4 select: AoA median |err| 4.3/5.0 deg, RMSE 17.8/41.4; K=5 avg: AoA median |err| 5.7/21.4 deg, RMSE 17.6/36.6; K=5 sharp: AoA median |err| 4.5/15.1 deg, RMSE 17.2/34.9; K=5 select: AoA median |err| 4.5/4.3 deg, RMSE 17.9/39.6; K=6 avg: AoA median |err| 5.2/22.3 deg, RMSE 15.7/35.7; K=6 sharp: AoA median |err| 4.5/14.8 deg, RMSE 14.9/32.5; K=6 select: AoA median |err| 4.8/3.0 deg, RMSE 17.9/34.4; K=7 avg: AoA median |err| 5.0/19.5 deg, RMSE 15.7/32.9; K=7 sharp: AoA median |err| 4.3/12.1 deg, RMSE 15.1/29.8; K=7 select: AoA median |err| 4.0/3.2 deg, RMSE 16.6/32.5; K=8 avg: AoA median |err| 5.2/22.2 deg, RMSE 15.9/33.9; K=8 sharp: AoA median |err| 4.2/14.2 deg, RMSE 15.3/30.5; K=8 select: AoA median |err| 4.0/3.2 deg, RMSE 16.6/32.5; K=9 avg: AoA median |err| 5.2/21.6 deg, RMSE 16.1/34.1; K=9 sharp: AoA median |err| 4.5/14.6 deg, RMSE 15.5/30.8; K=9 select: AoA median |err| 4.0/3.2 deg, RMSE 16.6/32.5; K=10 avg: AoA median |err| 5.4/21.1 deg, RMSE 16.2/33.4; K=10 sharp: AoA median |err| 4.8/14.2 deg, RMSE 15.7/30.1; K=10 select: AoA median |err| 4.0/2.7 deg, RMSE 16.6/31.4

Figures: `sweeps_<surface>.png` per surface, `sizes_comparison.png` across surfaces. All numbers: `summary.csv`, `validation_and_aoa.csv`.
