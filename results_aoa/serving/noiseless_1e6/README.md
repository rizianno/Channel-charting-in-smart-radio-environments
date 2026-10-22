# Serving with AoA estimated from few EMS acquisitions (noise `noiseless_1e6`)

Generated 2026-10-22 10:01 by `experiments/Main_AoA_Serving.m` (branch aoa-serving), runtime 32 min.

## What is simulated
- **Serving codebook:** N = 30 entries per EMS, phase step between adjacent elements from -168 to 180 deg/element in steps of 12 (entry 15 = zero phase = baseline). For lambda/2 spacing the beam is steered to sin(theta) = sin(theta_specular) + step/180; specular direction of EMS1 / EMS2 = -72.8 / -49.0 deg.
- **Training:** K entries spread evenly over the codebook (K = [1 2 3 4 7]), e.g. K = 3 -> entries [5 15 25]. `diag`: acquisition k sets both EMSs to training entry k (K acquisitions); `grid`: all pairs (K^2 acquisitions).
- **AoA estimation:** one joint t-SNE chart (azimuth from EMS1 and EMS2) of all 3200 users, 305 anchors (15% of the users in front of both surfaces); one chart per training acquisition, then local fusion weighted by a label-free quality (no ground truth of the served users). Noise: Exact covariance (no estimation noise) with the OLD fixed 1e-6 regularization, exactly as the original code and the paper's position experiment (Main_Position_CodebookSize on main): for like-for-like comparisons with them.
- **Look-up tables (angle -> entry, per EMS):** `physics` = steering law with the known BS direction; `data` = majority best entry of the 7 anchors with the closest true angle (offline calibration with the labelled users).
- **Link:** uplink, 23 dBm, noise -92 dBm, best rank-1 beamforming: SNR = P/N0 * lambda_max(H'H) with the TRUE channel H = direct + EMS1 + EMS2. SE = log2(1 + SNR). Coverage = share of users with SNR >= 0 dB (also 10 dB in the CSV).
- **Compared schemes:** `baseline` (zero phase on both EMSs), `best_training` (best of the training configurations measured for the user, no angle estimation), `proposed_*` (t-SNE AoA + LUT), `genie_*` (true AoA + LUT: the cost of the angle error), `full_sweep` (best of all N^2 configurations, upper bound).
- **Scored users:** unlabelled users only; groups: all, in front of both surfaces, without a direct path.

## Mode `diag`: median SNR gain vs baseline [dB] / mean SE [bit/s/Hz] / coverage SNR>=0 dB [%]

**all users**
| scheme | K=1 | K=2 | K=3 | K=4 | K=7 |
|---|---|---|---|---|---|
| baseline | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| best_training | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| proposed_physics | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| proposed_data | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| genie_physics | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| genie_data | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| full_sweep | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |

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
| baseline | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| best_training | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| proposed_physics | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| proposed_data | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| genie_physics | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| genie_data | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| full_sweep | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |

AoA error of the estimate (users in front of both, unlabelled), RMSE EMS1/EMS2 [deg]: K=1: 27.0/28.0, K=2: 25.8/24.0, K=3: 24.9/23.9, K=4: 24.1/22.1, K=7: 24.3/22.2

## Mode `grid`: median SNR gain vs baseline [dB] / mean SE [bit/s/Hz] / coverage SNR>=0 dB [%]

**all users**
| scheme | K=1 | K=2 | K=3 | K=4 | K=7 |
|---|---|---|---|---|---|
| baseline | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| best_training | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| proposed_physics | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| proposed_data | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| genie_physics | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| genie_data | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |
| full_sweep | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 | NaN / 1.54 / 27 |

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
| baseline | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| best_training | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| proposed_physics | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| proposed_data | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| genie_physics | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| genie_data | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |
| full_sweep | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 | NaN / 0.00 / 0 |

AoA error of the estimate (users in front of both, unlabelled), RMSE EMS1/EMS2 [deg]: K=1: 27.0/28.0, K=2: 24.4/21.6, K=3: 24.0/22.3, K=4: 23.8/21.1, K=7: 23.8/21.6

Figures: `serving_vs_acquisitions_<mode>.png` (SNR gain, SE, coverage vs number of training acquisitions), `snr_cdf_K3_<mode>.png` (SNR CDFs with 3 training entries). All numbers: `summary.csv`.
