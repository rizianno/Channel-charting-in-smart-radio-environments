# Serving with AoA estimated from few EMS acquisitions (noise `est_G60dB`)

Generated 2026-10-22 10:04 by `experiments/Main_AoA_Serving.m` (branch aoa-serving), runtime 1 min.

## What is simulated
- **Serving codebook:** N = 30 entries per EMS, phase step between adjacent elements from -168 to 180 deg/element in steps of 12 (entry 15 = zero phase = baseline). For lambda/2 spacing the beam is steered to sin(theta) = sin(theta_specular) + step/180; specular direction of EMS1 / EMS2 = -72.8 / -49.0 deg.
- **Training:** K entries spread evenly over the codebook (K = [1 2 3 4 7]), e.g. K = 3 -> entries [5 15 25]. `diag`: acquisition k sets both EMSs to training entry k (K acquisitions); `grid`: all pairs (K^2 acquisitions).
- **AoA estimation:** one joint t-SNE chart (azimuth from EMS1 and EMS2) of all 3200 users, 305 anchors (15% of the users in front of both surfaces); one chart per training acquisition, then local fusion weighted by a label-free quality (no ground truth of the served users). Noise: Estimated channels as est_G0dB, with a 60 dB pilot processing gain (noise variance on each estimate divided by 10^6).
- **Look-up tables (angle -> entry, per EMS):** `physics` = steering law with the known BS direction; `data` = majority best entry of the 7 anchors with the closest true angle (offline calibration with the labelled users).
- **Link:** uplink, 23 dBm, noise -92 dBm, best rank-1 beamforming: SNR = P/N0 * lambda_max(H'H) with the TRUE channel H = direct + EMS1 + EMS2. SE = log2(1 + SNR). Coverage = share of users with SNR >= 0 dB (also 10 dB in the CSV).
- **Compared schemes:** `baseline` (zero phase on both EMSs), `best_training` (best of the training configurations measured for the user, no angle estimation), `proposed_*` (t-SNE AoA + LUT), `genie_*` (true AoA + LUT: the cost of the angle error), `full_sweep` (best of all N^2 configurations, upper bound).
- **Scored users:** unlabelled users with a non-zero channel in at least one configuration (13 users have no channel at all and cannot be served by any scheme); groups: all, in front of both surfaces, without a direct path.

## Mode `diag`: median SNR gain vs baseline [dB] / mean SE [bit/s/Hz] / coverage SNR>=0 dB [%]

**all users**
| scheme | K=1 | K=2 | K=3 | K=4 | K=7 |
|---|---|---|---|---|---|
| baseline | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| best_training | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| proposed_physics | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| proposed_data | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| genie_physics | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| genie_data | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| full_sweep | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |

**in front of both surfaces**
| scheme | K=1 | K=2 | K=3 | K=4 | K=7 |
|---|---|---|---|---|---|
| baseline | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| best_training | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| proposed_physics | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| proposed_data | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| genie_physics | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| genie_data | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| full_sweep | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |

**without direct path**
| scheme | K=1 | K=2 | K=3 | K=4 | K=7 |
|---|---|---|---|---|---|
| baseline | 0.0 / 0.00 / 0 | 0.0 / 0.00 / 0 | 0.0 / 0.00 / 0 | 0.0 / 0.00 / 0 | 0.0 / 0.00 / 0 |
| best_training | 0.0 / 0.00 / 0 | 1.5 / 0.00 / 0 | 0.0 / 0.00 / 0 | 13.8 / 0.00 / 0 | 13.9 / 0.00 / 0 |
| proposed_physics | 2.9 / 0.00 / 0 | 2.4 / 0.00 / 0 | 1.8 / 0.00 / 0 | 2.2 / 0.00 / 0 | 2.2 / 0.00 / 0 |
| proposed_data | 3.1 / 0.00 / 0 | 2.1 / 0.00 / 0 | -1.8 / 0.00 / 0 | 2.0 / 0.00 / 0 | 1.9 / 0.00 / 0 |
| genie_physics | 25.7 / 0.00 / 0 | 25.7 / 0.00 / 0 | 25.7 / 0.00 / 0 | 25.7 / 0.00 / 0 | 25.7 / 0.00 / 0 |
| genie_data | 25.4 / 0.00 / 0 | 25.4 / 0.00 / 0 | 25.4 / 0.00 / 0 | 25.4 / 0.00 / 0 | 25.4 / 0.00 / 0 |
| full_sweep | 25.8 / 0.00 / 0 | 25.8 / 0.00 / 0 | 25.8 / 0.00 / 0 | 25.8 / 0.00 / 0 | 25.8 / 0.00 / 0 |

AoA error of the estimate (users in front of both, unlabelled), RMSE EMS1/EMS2 [deg]: K=1: 25.6/27.7, K=2: 23.2/26.1, K=3: 25.4/26.4, K=4: 22.8/22.8, K=7: 22.6/21.6

## Mode `grid`: median SNR gain vs baseline [dB] / mean SE [bit/s/Hz] / coverage SNR>=0 dB [%]

**all users**
| scheme | K=1 | K=2 | K=3 | K=4 | K=7 |
|---|---|---|---|---|---|
| baseline | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| best_training | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| proposed_physics | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| proposed_data | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| genie_physics | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| genie_data | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |
| full_sweep | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 | 0.0 / 1.55 / 27 |

**in front of both surfaces**
| scheme | K=1 | K=2 | K=3 | K=4 | K=7 |
|---|---|---|---|---|---|
| baseline | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| best_training | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| proposed_physics | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| proposed_data | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| genie_physics | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| genie_data | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |
| full_sweep | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 | 0.0 / 1.74 / 33 |

**without direct path**
| scheme | K=1 | K=2 | K=3 | K=4 | K=7 |
|---|---|---|---|---|---|
| baseline | 0.0 / 0.00 / 0 | 0.0 / 0.00 / 0 | 0.0 / 0.00 / 0 | 0.0 / 0.00 / 0 | 0.0 / 0.00 / 0 |
| best_training | 0.0 / 0.00 / 0 | 1.6 / 0.00 / 0 | 0.0 / 0.00 / 0 | 13.9 / 0.00 / 0 | 14.0 / 0.00 / 0 |
| proposed_physics | 2.9 / 0.00 / 0 | 2.6 / 0.00 / 0 | 2.4 / 0.00 / 0 | 1.5 / 0.00 / 0 | 1.9 / 0.00 / 0 |
| proposed_data | 3.1 / 0.00 / 0 | 2.2 / 0.00 / 0 | 1.9 / 0.00 / 0 | 1.1 / 0.00 / 0 | 1.5 / 0.00 / 0 |
| genie_physics | 25.7 / 0.00 / 0 | 25.7 / 0.00 / 0 | 25.7 / 0.00 / 0 | 25.7 / 0.00 / 0 | 25.7 / 0.00 / 0 |
| genie_data | 25.4 / 0.00 / 0 | 25.4 / 0.00 / 0 | 25.4 / 0.00 / 0 | 25.4 / 0.00 / 0 | 25.4 / 0.00 / 0 |
| full_sweep | 25.8 / 0.00 / 0 | 25.8 / 0.00 / 0 | 25.8 / 0.00 / 0 | 25.8 / 0.00 / 0 | 25.8 / 0.00 / 0 |

AoA error of the estimate (users in front of both, unlabelled), RMSE EMS1/EMS2 [deg]: K=1: 25.6/27.7, K=2: 23.6/23.6, K=3: 24.2/25.4, K=4: 22.9/20.8, K=7: 22.4/20.9

Figures: `serving_vs_acquisitions_<mode>.png` (SNR gain, SE, coverage vs number of training acquisitions), `snr_cdf_K3_<mode>.png` (SNR CDFs with 3 training entries). All numbers: `summary.csv`.
