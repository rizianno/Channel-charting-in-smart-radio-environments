# EMS beam selection vs number of sweeps (noise `est_G60dB`)

Generated 2026-10-22 15:42 by `experiments/Main_AoA_Serving_Sweeps.m` (branch aoa-serving), runtime 31 min.

## Codebook from the beamwidth
- 60 elements along the steering axis at lambda/2 -> aperture 30 lambda (0.30 m at 30 GHz). Half-power beamwidth in sin(theta): 0.886 lambda / D = 0.0295 (1.7 deg near the normal, 3.4 deg at 60 deg).
- Covering sin(theta) in [-1, 1] with beams crossing at -3 dB: N = 68 entries per EMS, phase step 5.29 deg/element, over -180..180.
- 180 elements per side would give width 0.0098 and about 204 entries; here only its gain is emulated (+19.1 dB = (180/60)^4), the beams stay those of the 60-element data.
- Far-field distance 2 D^2 / lambda = 18 m: users closer than that are in the near field of the surface, where the steered (far-field) beams lose gain.
- Specular direction of EMS1 / EMS2 (BS mirrored): -72.8 / -49.0 deg.

## Schemes (all vs number of sweeps, codebook fixed)
- `sequential`: entries 1, 2, 3, ... (worst baseline). `spread`: s entries spread evenly. `hierarchical`: 4 spread entries, then the two neighbours at half the spacing around the current best, down to spacing 1.
- `tsne`: K spread training sweeps (K = [1 2 3 4 5 6 7 8], both EMSs on the same entry) -> t-SNE azimuth (joint chart, all 3200 users, 305 anchors, label-free local fusion) -> look-up table (steering law with the known BS direction). `tsne+refine` (R = [2 4 6]): then measure the LUT entry and R neighbours (R/2 per side): K + R + 1 sweeps. `spread+refine`: B spread sweeps + R neighbours of the best: B + R sweeps.
- Upper bounds: `exhaustive` (all 68 entries per EMS) and `exhaustive_joint` (all 4624 pairs). `true AoA + LUT`: the look-up table with the exact angle (no sweeps).
- Selection by measurement uses each EMS path's own gain (the BS separates the two surfaces spatially: each EMS->BS channel is rank one), noise-free (in favour of the sweeping baselines); t-SNE uses estimated channels. The served SNR uses the full true channel (direct + both EMSs).
- Users: unlabelled with any channel; LoS 811, reflection only 1504, no BS-UE channel 885 (all users). Link: 23 dBm, -92 dBm, SNR = P/N0 lambda_max(H'H).

## EMS gain +0.0 dB: median SNR gain vs zero phase [dB], users without LoS
| sweeps | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 16 | 20 | 30 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -0.0 | -0.0 | -0.0 | -0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| spread | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.1 |
| hierarchical | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.1 | 0.1 | 0.1 | 0.1 |
| spread+refine (R=2) | - | - | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.1 | 0.1 | 0.1 | 0.1 |
| spread+refine (R=4) | - | - | - | - | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.1 | 0.1 | 0.1 |
| spread+refine (R=6) | - | - | - | - | - | - | 0.0 | 0.1 | 0.0 | 0.1 | 0.1 | 0.1 |
| tsne | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | - | - | - | - | - |
| tsne+refine (R=2) | - | - | - | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | - | - | - | - |
| tsne+refine (R=4) | - | - | - | - | - | 0.0 | 0.0 | 0.0 | 0.0 | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | - | 0.0 | 0.0 | 0.0 | - | - | - |
| true AoA + LUT (0 sweeps) | 0.1 | | | | | | | | | | | |
| exhaustive (68 sweeps) | 0.1 | | | | | | | | | | | |
| exhaustive_joint (4624 sweeps) | 0.2 | | | | | | | | | | | |

## EMS gain +19.1 dB: median SNR gain vs zero phase [dB], users without LoS
| sweeps | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 16 | 20 | 30 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sequential | -0.0 | -0.0 | -0.0 | -0.0 | -0.0 | -0.0 | -0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.2 |
| spread | 0.0 | 0.0 | 0.0 | 0.2 | 0.0 | 0.3 | 0.2 | 0.3 | 0.6 | 0.9 | 1.4 | 1.3 |
| hierarchical | 0.0 | 0.0 | 0.0 | 0.2 | 0.2 | 0.2 | 0.4 | 1.0 | 1.8 | 2.6 | 2.6 | 2.6 |
| spread+refine (R=2) | - | - | 0.0 | 0.0 | 0.0 | 0.6 | 0.6 | 0.5 | 1.2 | 1.2 | 1.4 | 1.8 |
| spread+refine (R=4) | - | - | - | - | 0.0 | 0.0 | 1.0 | 1.0 | 1.5 | 2.1 | 1.9 | 2.0 |
| spread+refine (R=6) | - | - | - | - | - | - | 0.1 | 1.4 | 1.4 | 2.3 | 2.1 | 2.1 |
| tsne | 0.0 | 0.1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | - | - | - | - | - |
| tsne+refine (R=2) | - | - | - | 0.1 | 0.3 | 0.3 | 0.3 | 0.3 | - | - | - | - |
| tsne+refine (R=4) | - | - | - | - | - | 0.3 | 0.5 | 0.6 | 0.6 | - | - | - |
| tsne+refine (R=6) | - | - | - | - | - | - | 0.5 | 0.7 | 0.9 | - | - | - |
| true AoA + LUT (0 sweeps) | 3.0 | | | | | | | | | | | |
| exhaustive (68 sweeps) | 3.1 | | | | | | | | | | | |
| exhaustive_joint (4624 sweeps) | 3.5 | | | | | | | | | | | |

Figures: `sweeps_gain<g>dB.png` (median / mean SNR gain and SE gain vs sweeps for users without LoS, reflection-only and no-channel users). All numbers, also LoS and all users and coverage: `summary.csv`.
