# Position estimation vs codebook size (1x1 ... 5x5)

Paper position experiment (3200 users, 480 anchors, per-configuration t-SNE + fusion, noiseless, 1e-6 regularization) with K x K codebooks instead of 11x11, for two phase-step ranges. Runtime 27 min. See each range folder for details.

- [`narrow`](narrow/): phase steps 0..40 deg/element (the range of the paper's 11x11 codebook)
- [`wide`](wide/): phase steps spread evenly over -180..180 deg/element (beams over all directions)

| range | K | local fusion mean error [m] | p90 [m] | typical single mean [m] | fusion gain [m] | local fusion TW | local fusion CT |
|---|---|---|---|---|---|---|---|
| narrow | 1x1 | 13.5 | 29.4 | 13.5 | 0.0 | 0.878 | 0.875 |
| narrow | 2x2 | 9.0 | 19.6 | 12.7 | 3.6 | 0.933 | 0.938 |
| narrow | 3x3 | 8.0 | 16.2 | 12.4 | 4.4 | 0.943 | 0.948 |
| narrow | 4x4 | 7.8 | 15.4 | 12.4 | 4.6 | 0.943 | 0.947 |
| narrow | 5x5 | 7.5 | 15.1 | 12.3 | 4.8 | 0.948 | 0.952 |
| wide | 1x1 | 13.5 | 29.4 | 13.5 | 0.0 | 0.878 | 0.875 |
| wide | 2x2 | 10.4 | 24.4 | 14.0 | 3.6 | 0.913 | 0.917 |
| wide | 3x3 | 9.2 | 22.0 | 13.4 | 4.3 | 0.930 | 0.934 |
| wide | 4x4 | 9.1 | 20.2 | 13.7 | 4.6 | 0.931 | 0.934 |
| wide | 5x5 | 9.0 | 20.8 | 13.9 | 4.8 | 0.934 | 0.937 |

`comparison.png`: narrow vs wide, local fusion and typical single configuration.
