# Asilomar abstract compliance

Abstract checked: `Dropbox/Apps/Overleaf/Asilomar Conference - Channel Charting/main.tex`
("Towards Real-Time RIS Configuration: A Channel Charting Approach", version of 20 Sep).
Purpose: keep the simulations, paper and slides in line with what the abstract claims.

## Mapping of the abstract's two acquisition strategies onto the current work

The abstract proposes (i) a baseline drawing the K_acq acquisition configurations uniformly at
random from the codebook, and (ii) an *informed* strategy using a prior PDF p(theta) of the user
angles to concentrate acquisitions in the likely angular regions, claiming (ii) significantly
outperforms (i).

Evidence for (ii) from our runs: a chart places users correctly only if one of the acquisitions
illuminates them. EMS2's full-NLOS users sit in three angle ranges (about -80, -30, +35 deg); the
nested spread order hits them only by chance at K = 5 and K = 6 (cluster diagnostic,
`results_aoa/serving_nlos_size_nestedFusion_TEMP/cluster_check.txt`).

Planned pipeline:

1. Prior per surface: p_e(theta) from the anchors' known azimuths (the 15 % labelled users),
   i.e. "learned offline from previously collected CSI".
2. Method (ii), informed: per surface, choose the K beams covering the most prior mass in their
   main lobe (greedy max-coverage on the sin(theta) grid, or beams at prior quantiles). The two
   EMSs get different entries, so one slot = one pair (EMS1 entry, EMS2 entry): K acquisitions,
   not K^2.
3. Method (i), random: K entries drawn uniformly per surface, averaged over ~20 draws; nested
   spread kept as a deterministic reference.
4. Unchanged afterwards: per-user chart selection (max label-free trustworthiness) -> LUT ->
   optional refinement / hierarchical search; compared with hierarchical search and exhaustive.

## Conflicts between the abstract and the actual work

| Abstract says | What we actually do | Solution |
|---|---|---|
| Embeddings fused by weighted averaging in the latent space | Weighted averaging was worse (near-uniform weights); we use per-user selection of the chart with the highest trustworthiness | Present selection as weighted averaging with one-hot weights (softmax over trustworthiness, temperature -> 0); show the soft average as the variant that fails, with the near-uniform-weights explanation |
| Baseline (i): acquisitions drawn uniformly at random | Never run (we used same-entry, grid, nested) | Run the random baseline (about 20 draws) |
| Informed strategy (ii) significantly outperforms (i) | Not yet shown | Run (ii); if the gap is small, soften to "reduces the number of acquisitions needed" |
| LUT built offline from the labelled portion of the training set | LUT from the steering law (physics) | Add an anchor-based LUT (anchors' offline exhaustive best codeword mapped by angle: nearest anchor or binned majority) and report it next to the physics LUT |
| Codebook constructed following the procedure of [OurPaperA] | New full-range phase-gradient codebook (-180..180 deg steps, N from the beamwidth: 68 to 272); the old one steered only a ~16-28 deg sector | Describe our codebook explicitly as a generalisation; fix lambda/2 (data) vs lambda/4 (OurPaperA Remark 1) |
| "A RIS", "a single-antenna UE" | Two EMSs; dataset UE has 2 antennas, used in the covariance features | Say "one or more surfaces, each configured independently with its own chart"; declare 2 UE antennas or rerun with one |
| Baselines: random codeword (no acquisitions) and GFlowNet adapted as a reference | Neither implemented; we have exhaustive, sequential, spread and hierarchical search | Add random codeword (trivial); replace GFlowNet by the closed-form per-element optimum (EMS-to-BS link is rank one, so co-phasing is optimal and upper-bounds any per-element method), plus the complexity argument |
| Results on received power, signalling overhead and spectral efficiency | SNR gain (median, 10th percentile); absolute SNR via the EMS is very low (about -36 dB at 60 x 60 with exhaustive), so SE is about 0 | Decide the link budget (Tx power / bandwidth normalisation) or report SE at 240 x 240; add overhead-aware effective SE (1 - K_acq T_slot / T_frame) log2(1 + SNR) vs K_acq; add received-power CDFs |
| Implicitly all users served | Charts and scoring on full-NLOS users only | State the scope: the RIS serves users without a usable direct link, detected at the BS from direct-link power |
| Azimuth and/or elevation | Azimuth only, no vertical steering | Fine under "and/or"; state that elevation is fixed |
| Approach exhaustive search; saturation point small relative to K | Consistent: selection saturates at K about 6 of 68 on the real surface | Keep; add top-1 codeword accuracy (`beam_agreement`) |

Caveat for paper and slides: hierarchical search beats t-SNE from about 10 sweeps and in the
worst 10 % of users. Frame the contribution as the few-acquisition regime, or as t-SNE warm-starting
hierarchical search.

## Suggested order

1. Informed and random acquisition modes (core of the abstract).
2. Anchor-based LUT and random-codeword baseline.
3. Per-element upper bound.
4. SE and CDF plots, after the link-budget decision.
