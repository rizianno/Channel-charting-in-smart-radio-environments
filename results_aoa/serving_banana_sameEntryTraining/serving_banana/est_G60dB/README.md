# Banana experiment: charts from all users WITHOUT LoS, beam selection vs number of sweeps (noise `est_G60dB`)

Generated 2026-10-22 23:37 by `experiments/Main_AoA_Serving_NLoS_Size.m` (branch aoa-serving), runtime 31 min.

## Setup
- **Users (banana):** every user without LoS (2389): reflection only (1504, BS link only via a wall: direct-channel power more than 36 dB below free space) and full NLOS (885, no BS-UE channel). Their BS channel is part of the t-SNE input and of the served SNR. The full-NLOS results compare directly with results_aoa/serving_nlos_size (charts from full-NLOS users only) and results_aoa/serving_sweeps (charts from all users).
- **t-SNE:** one 1-D chart per EMS with the users of this population it reaches (EMS1: 1353 users, 203 anchors, perplexity 150; EMS2: 1951, 293, 150). Scored: reachable by at least one EMS and anchor in neither chart.
- **Surfaces:** M x M elements, lambda/2, steering along one axis only (vertical: same phase on every column). `M60_raytraced` = the dataset; `M60/120/180/240` = synthesized with spherical waves from the UE, EMS and BS positions, per-user strengths fitted on the 60 x 60 data. Codebook from the beamwidth: N = 68/136/204/272.
- **Schemes vs number of sweeps:** sequential, spread, hierarchical (4 spread entries, then halving), spread + R neighbours, t-SNE AoA + LUT (K = [1 2 3 4 5 6 7 8] spread training sweeps, both EMSs on the same entry), t-SNE + LUT entry + R neighbours (R = [2 4 6]), exhaustive (N), true AoA + LUT.
- Measurement-based selection uses each EMS path's own gain (noise-free, in favour of the sweeping baselines); the served SNR uses the true channel; 23 dBm, -92 dBm.

## Validation of the synthesized channels (60 x 60 vs ray-traced)

- EMS1: element grid orientation (m, n) = (+1, +1), phase sign +1; |correlation| of the EMS->BS channel with the model 0.998; UE->EMS channels: median |correlation| 0.998, above 0.9 for 100% of the 1353 reachable full-NLOS users.
- EMS2: element grid orientation (m, n) = (+1, -1), phase sign +1; |correlation| of the EMS->BS channel with the model 0.999; UE->EMS channels: median |correlation| 0.999, above 0.9 for 100% of the 1951 reachable full-NLOS users.
- validation EMS1 (1353 users): best entry same as ray-traced for 98%, within 1 entry for 100%; ray-traced gain at the synthetic best entry: median 0.0 dB below the ray-traced best; best gain synthetic vs ray-traced: median -0.0 dB
- validation EMS2 (1951 users): best entry same as ray-traced for 98%, within 1 entry for 100%; ray-traced gain at the synthetic best entry: median 0.0 dB below the ray-traced best; best gain synthetic vs ray-traced: median +0.1 dB

## Median SNR gain vs zero phase [dB] / 10th percentile [dB], full-NLOS users

**M60_raytraced** (codebook 68, 737 users)
| scheme | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -6.1 / -16.3 | -4.4 / -13.5 | -3.2 / -12.3 | -2.7 / -11.4 | -2.1 / -10.5 | -2.0 / -10.0 | -2.0 / -9.9 | -1.9 / -9.7 | -1.5 / -9.3 | -0.8 / -8.4 |
| spread | 0.0 / 0.0 | 2.1 / -6.1 | 0.3 / 0.0 | 15.7 / 3.4 | 12.9 / 4.2 | 14.1 / -2.7 | 16.4 / 4.6 | 19.9 / 9.7 | 22.4 / 8.5 | 24.6 / 13.0 |
| hierarchical | 0.0 / 0.0 | 2.1 / -6.1 | 0.3 / 0.0 | 15.7 / 3.4 | 15.8 / 4.3 | 18.5 / 9.8 | 23.8 / 12.8 | 28.0 / 17.1 | 30.6 / 23.0 | 30.6 / 23.0 |
| spread+refine (R=6) | - | - | - | - | - | 9.2 / -3.3 | 28.6 / 13.8 | 23.9 / 9.8 | 29.7 / 20.6 | 30.1 / 19.8 |
| tsne | 5.6 / -7.2 | 5.5 / -7.3 | 5.4 / -8.8 | 6.9 / -7.0 | 7.7 / -6.5 | 7.0 / -7.2 | - | - | - | - |
| tsne+refine (R=2) | - | - | - | 9.5 / -3.2 | 9.0 / -3.2 | 8.6 / -3.6 | 10.7 / -2.8 | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | 17.6 / -1.0 | 15.6 / -0.3 | 17.3 / -0.0 | - | - |
| true AoA + LUT (0 sweeps) | 32.0 / 24.0 | | | | | | | | | |
| exhaustive (68 sweeps) | 32.1 / 24.2 | | | | | | | | | |
Zero-phase median SNR -69.5 dB, exhaustive median SNR -36.1 dB.
t-SNE angle error / LUT: K=1: AoA median |err| 5.2/6.9 deg, RMSE 20.4/37.3; LUT >1 beam off 51%/60%; diverged 0.0%; K=2: AoA median |err| 5.2/7.8 deg, RMSE 20.7/28.9; LUT >1 beam off 50%/66%; diverged 0.0%; K=3: AoA median |err| 5.7/7.8 deg, RMSE 20.7/26.7; LUT >1 beam off 53%/66%; diverged 0.0%; K=4: AoA median |err| 5.6/6.9 deg, RMSE 21.0/25.3; LUT >1 beam off 51%/63%; diverged 0.0%; K=5: AoA median |err| 5.4/7.8 deg, RMSE 20.5/28.1; LUT >1 beam off 51%/66%; diverged 0.0%; K=6: AoA median |err| 4.8/8.0 deg, RMSE 20.2/26.0; LUT >1 beam off 49%/67%; diverged 0.0%; K=7: AoA median |err| 5.6/7.2 deg, RMSE 20.4/26.9; LUT >1 beam off 53%/66%; diverged 0.0%; K=8: AoA median |err| 5.2/6.7 deg, RMSE 19.9/25.5; LUT >1 beam off 52%/65%; diverged 0.0%

**M60** (codebook 68, 737 users)
| scheme | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -6.1 / -16.1 | -4.4 / -13.5 | -3.2 / -12.3 | -2.7 / -11.4 | -2.1 / -10.5 | -2.0 / -9.9 | -2.0 / -9.9 | -1.9 / -9.7 | -1.5 / -9.3 | -0.7 / -8.3 |
| spread | 0.0 / 0.0 | 2.1 / -6.2 | 0.4 / 0.0 | 16.0 / 4.6 | 12.9 / 4.3 | 14.3 / -2.0 | 16.3 / 5.0 | 20.1 / 10.1 | 22.5 / 8.8 | 24.6 / 13.0 |
| hierarchical | 0.0 / 0.0 | 2.1 / -6.2 | 0.4 / 0.0 | 16.0 / 4.6 | 16.1 / 4.7 | 18.7 / 10.2 | 23.8 / 13.0 | 28.0 / 17.6 | 30.7 / 23.2 | 30.7 / 23.2 |
| spread+refine (R=6) | - | - | - | - | - | 9.1 / -3.3 | 28.7 / 14.2 | 23.9 / 10.1 | 29.7 / 21.0 | 30.0 / 20.2 |
| tsne | 6.7 / -6.8 | 5.3 / -8.5 | 5.6 / -8.2 | 5.3 / -8.9 | 6.8 / -7.9 | 5.9 / -8.5 | - | - | - | - |
| tsne+refine (R=2) | - | - | - | 12.6 / -1.5 | 9.1 / -2.3 | 10.1 / -2.7 | 10.4 / -2.7 | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | 22.5 / -0.4 | 17.6 / 0.2 | 20.0 / 0.0 | - | - |
| true AoA + LUT (0 sweeps) | 32.0 / 24.0 | | | | | | | | | |
| exhaustive (68 sweeps) | 32.2 / 24.2 | | | | | | | | | |
Zero-phase median SNR -69.5 dB, exhaustive median SNR -36.2 dB.
t-SNE angle error / LUT: K=1: AoA median |err| 5.7/6.0 deg, RMSE 20.2/31.6; LUT >1 beam off 53%/59%; diverged 0.0%; K=2: AoA median |err| 5.4/7.5 deg, RMSE 20.5/31.2; LUT >1 beam off 51%/66%; diverged 0.0%; K=3: AoA median |err| 5.8/7.3 deg, RMSE 20.0/27.0; LUT >1 beam off 55%/65%; diverged 0.0%; K=4: AoA median |err| 5.8/7.3 deg, RMSE 21.2/23.2; LUT >1 beam off 54%/66%; diverged 0.0%; K=5: AoA median |err| 5.5/6.8 deg, RMSE 20.3/27.5; LUT >1 beam off 54%/64%; diverged 0.0%; K=6: AoA median |err| 5.1/7.7 deg, RMSE 19.9/27.5; LUT >1 beam off 50%/67%; diverged 0.0%; K=7: AoA median |err| 5.6/7.1 deg, RMSE 20.3/27.0; LUT >1 beam off 53%/66%; diverged 0.0%; K=8: AoA median |err| 5.6/7.8 deg, RMSE 19.9/25.8; LUT >1 beam off 55%/67%; diverged 0.0%

**M120** (codebook 136, 737 users)
| scheme | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -6.3 / -15.6 | -4.8 / -13.6 | -3.5 / -12.0 | -2.8 / -11.0 | -2.1 / -10.0 | -2.0 / -9.6 | -2.0 / -9.6 | -2.0 / -9.6 | -1.9 / -9.5 | -1.7 / -9.5 |
| spread | 0.0 / 0.0 | 2.1 / -5.7 | 0.2 / 0.0 | 15.1 / 6.6 | 11.9 / 4.0 | 13.9 / 6.1 | 17.3 / 6.9 | 19.1 / 10.6 | 22.1 / 6.9 | 23.2 / 14.3 |
| hierarchical | 0.0 / 0.0 | 2.1 / -5.7 | 0.2 / 0.0 | 15.1 / 6.6 | 15.6 / 7.2 | 20.2 / 12.1 | 24.2 / 15.3 | 29.1 / 19.8 | 36.5 / 30.1 | 36.5 / 30.1 |
| spread+refine (R=6) | - | - | - | - | - | 8.2 / -3.2 | 29.8 / 11.1 | 19.9 / 8.7 | 27.9 / 10.2 | 34.4 / 19.5 |
| tsne | 6.3 / -7.9 | 7.7 / -9.3 | 5.4 / -7.6 | 5.8 / -7.7 | 6.2 / -7.6 | 6.8 / -6.9 | - | - | - | - |
| tsne+refine (R=2) | - | - | - | 10.5 / -2.9 | 9.0 / -1.6 | 10.5 / -2.4 | 8.9 / -2.7 | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | 14.8 / -0.8 | 13.0 / 0.3 | 14.4 / -0.5 | - | - |
| true AoA + LUT (0 sweeps) | 36.8 / 30.1 | | | | | | | | | |
| exhaustive (136 sweeps) | 37.5 / 31.0 | | | | | | | | | |
Zero-phase median SNR -66.3 dB, exhaustive median SNR -26.6 dB.
t-SNE angle error / LUT: K=1: AoA median |err| 5.4/6.4 deg, RMSE 21.7/30.6; LUT >1 beam off 70%/74%; diverged 0.0%; K=2: AoA median |err| 5.5/7.2 deg, RMSE 21.3/29.2; LUT >1 beam off 71%/77%; diverged 0.0%; K=3: AoA median |err| 5.8/7.3 deg, RMSE 20.6/24.1; LUT >1 beam off 72%/79%; diverged 0.0%; K=4: AoA median |err| 6.5/7.1 deg, RMSE 20.3/22.4; LUT >1 beam off 74%/79%; diverged 0.0%; K=5: AoA median |err| 5.1/7.4 deg, RMSE 20.3/25.9; LUT >1 beam off 72%/79%; diverged 0.0%; K=6: AoA median |err| 5.3/7.9 deg, RMSE 20.0/25.2; LUT >1 beam off 71%/79%; diverged 0.0%; K=7: AoA median |err| 5.7/7.4 deg, RMSE 20.3/24.3; LUT >1 beam off 72%/79%; diverged 0.0%; K=8: AoA median |err| 5.4/7.2 deg, RMSE 20.1/24.5; LUT >1 beam off 72%/79%; diverged 0.0%

**M180** (codebook 204, 737 users)
| scheme | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -5.9 / -16.0 | -4.4 / -13.8 | -3.5 / -12.2 | -2.8 / -11.1 | -2.1 / -9.9 | -2.1 / -9.7 | -2.1 / -9.7 | -2.1 / -9.7 | -2.0 / -9.7 | -2.0 / -9.7 |
| spread | 0.0 / 0.0 | 2.3 / -6.1 | 0.0 / 0.0 | 15.8 / 5.9 | 11.8 / 6.9 | 15.7 / -0.7 | 16.5 / 5.9 | 19.5 / 9.3 | 22.0 / 12.2 | 24.8 / 14.7 |
| hierarchical | 0.0 / 0.0 | 2.3 / -6.1 | 0.0 / 0.0 | 15.8 / 5.9 | 16.2 / 6.5 | 21.7 / 11.3 | 27.0 / 15.5 | 32.2 / 21.1 | 39.3 / 31.7 | 39.3 / 31.7 |
| spread+refine (R=6) | - | - | - | - | - | 7.5 / -3.2 | 25.8 / 10.8 | 18.5 / 9.3 | 26.1 / 6.6 | 34.0 / 18.3 |
| tsne | 9.4 / -6.4 | 6.7 / -7.9 | 6.4 / -6.0 | 5.9 / -8.3 | 7.0 / -7.4 | 7.8 / -7.4 | - | - | - | - |
| tsne+refine (R=2) | - | - | - | 12.7 / -2.2 | 9.6 / -1.3 | 10.0 / -1.8 | 8.4 / -2.5 | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | 16.6 / -0.4 | 12.6 / 0.5 | 13.0 / 0.1 | - | - |
| true AoA + LUT (0 sweeps) | 40.0 / 32.8 | | | | | | | | | |
| exhaustive (204 sweeps) | 40.7 / 34.0 | | | | | | | | | |
Zero-phase median SNR -64.4 dB, exhaustive median SNR -21.9 dB.
t-SNE angle error / LUT: K=1: AoA median |err| 5.6/5.8 deg, RMSE 20.2/31.4; LUT >1 beam off 80%/80%; diverged 0.0%; K=2: AoA median |err| 5.4/7.8 deg, RMSE 20.7/26.1; LUT >1 beam off 80%/84%; diverged 0.0%; K=3: AoA median |err| 5.4/7.7 deg, RMSE 20.3/25.1; LUT >1 beam off 81%/84%; diverged 0.0%; K=4: AoA median |err| 5.7/6.8 deg, RMSE 20.3/22.0; LUT >1 beam off 80%/83%; diverged 0.0%; K=5: AoA median |err| 5.2/6.7 deg, RMSE 20.2/25.1; LUT >1 beam off 79%/84%; diverged 0.0%; K=6: AoA median |err| 5.8/6.8 deg, RMSE 20.4/24.5; LUT >1 beam off 81%/84%; diverged 0.0%; K=7: AoA median |err| 5.6/7.1 deg, RMSE 20.5/23.5; LUT >1 beam off 81%/85%; diverged 0.0%; K=8: AoA median |err| 4.8/7.7 deg, RMSE 20.0/24.1; LUT >1 beam off 80%/86%; diverged 0.0%

**M240** (codebook 272, 737 users)
| scheme | 1 | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -6.2 / -14.9 | -4.8 / -13.2 | -3.6 / -12.1 | -2.9 / -11.2 | -2.2 / -9.9 | -2.1 / -9.5 | -2.0 / -9.5 | -2.0 / -9.5 | -2.0 / -9.5 | -2.0 / -9.5 |
| spread | 0.0 / 0.0 | 2.0 / -6.0 | 0.4 / 0.0 | 16.2 / 6.6 | 13.1 / 5.2 | 14.9 / 7.0 | 17.9 / 5.8 | 20.1 / 10.8 | 22.5 / 13.1 | 25.5 / 15.3 |
| hierarchical | 0.0 / 0.0 | 2.0 / -6.0 | 0.4 / 0.0 | 16.2 / 6.6 | 16.6 / 8.1 | 22.4 / 13.7 | 29.1 / 18.8 | 34.3 / 24.8 | 39.4 / 31.0 | 42.1 / 35.4 |
| spread+refine (R=6) | - | - | - | - | - | 7.7 / -3.4 | 25.1 / 10.0 | 17.7 / 9.0 | 26.2 / 9.7 | 33.2 / 15.6 |
| tsne | 10.4 / -5.8 | 7.2 / -6.9 | 7.4 / -5.9 | 7.8 / -7.5 | 9.4 / -6.4 | 7.2 / -6.4 | - | - | - | - |
| tsne+refine (R=2) | - | - | - | 12.8 / -2.6 | 10.1 / -1.3 | 9.2 / -3.1 | 9.2 / -2.7 | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | 15.2 / -1.0 | 12.1 / 1.1 | 11.2 / -0.1 | - | - |
| true AoA + LUT (0 sweeps) | 41.8 / 34.7 | | | | | | | | | |
| exhaustive (272 sweeps) | 42.3 / 35.6 | | | | | | | | | |
Zero-phase median SNR -63.0 dB, exhaustive median SNR -19.3 dB.
t-SNE angle error / LUT: K=1: AoA median |err| 6.4/6.2 deg, RMSE 27.1/33.1; LUT >1 beam off 85%/86%; diverged 0.0%; K=2: AoA median |err| 5.2/7.4 deg, RMSE 21.1/25.7; LUT >1 beam off 87%/88%; diverged 0.0%; K=3: AoA median |err| 5.8/7.5 deg, RMSE 20.6/23.4; LUT >1 beam off 86%/88%; diverged 0.0%; K=4: AoA median |err| 5.5/6.2 deg, RMSE 20.3/21.0; LUT >1 beam off 85%/89%; diverged 0.0%; K=5: AoA median |err| 5.5/7.8 deg, RMSE 20.3/23.1; LUT >1 beam off 87%/88%; diverged 0.0%; K=6: AoA median |err| 4.9/6.6 deg, RMSE 20.1/22.1; LUT >1 beam off 85%/88%; diverged 0.0%; K=7: AoA median |err| 5.7/6.8 deg, RMSE 20.3/22.2; LUT >1 beam off 87%/89%; diverged 0.0%; K=8: AoA median |err| 6.1/6.4 deg, RMSE 19.9/22.7; LUT >1 beam off 87%/88%; diverged 0.0%

Figures: `sweeps_<surface>.png` per surface, `sizes_comparison.png` across surfaces. All numbers: `summary.csv`, `validation_and_aoa.csv`.
