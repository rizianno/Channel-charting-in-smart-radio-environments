function abl_readme(scn, CAT, noise, U, charts, allScn, outDir, T)
% ABL_README  README.md of one scenario folder: what is simulated, how it differs from the
% reference and from the scenarios that change exactly one factor, outputs and definitions.
    st = CAT.structures(strcmp({CAT.structures.id}, scn.structure));
    ft = CAT.features(strcmp({CAT.features.id}, scn.feature));
    dtxt = {'**off**: the direct BS-UE channel is removed (ideal: the BS rejects it completely); only the EMS-reflected paths remain.', ...
            '**on**: the direct BS-UE channel (LoS + one bounce, `H_LOS_1bounce`) is included, as in the reference.'};
    L = {};
    L{end+1} = sprintf('# %s', scn.name);
    L{end+1} = '';
    L{end+1} = sprintf('Noise model: `%s`. Codebook sizes: %s. Generated %s by `experiments/Main_AoA_Ablations.m`.', ...
                       noise.name, strjoin(arrayfun(@(k) sprintf('%dx%d',k,k), scn.Ks, 'UniformOutput', false), ', '), ...
                       char(datetime('now','Format','yyyy-MM-dd HH:mm')));
    L{end+1} = '';
    L{end+1} = '## What this simulation does';
    L{end+1} = sprintf('- **Users and charts (%s).** %s', scn.structure, st.text);
    for c = charts
        L{end+1} = sprintf('  - chart `%s`: %d users, %d anchors, label = azimuth from %s', c.id, ...
                           numel(c.users), c.nA, strjoin(arrayfun(@(e) sprintf('EMS%d',e), c.cols, 'UniformOutput', false), ' and ')); %#ok<AGROW>
    end
    L{end+1} = sprintf('- **Input feature (%s).** %s', scn.feature, ft.text);
    L{end+1} = sprintf('- **Direct path.** %s', dtxt{scn.direct+1});
    L{end+1} = sprintf('- **Noise (%s).** %s', noise.name, noise.text);
    L{end+1} = '- **Unchanged everywhere:** dataset, surface geometry, t-SNE (`tsne_sem`, perplexity 150), Log-Euclidean distance, 15% supervision, fusion rule (where there is fusion).';
    L{end+1} = '';
    L{end+1} = '| K | phase gradients per EMS [deg] | acquisitions |';
    L{end+1} = '|---|---|---|';
    for K = scn.Ks
        a = unique(T.acquisitions(T.K == K));
        L{end+1} = sprintf('| %dx%d | %s | %d |', K, K, mat2str(round(CAT.levels{K},2)), a(1)); %#ok<AGROW>
    end
    L{end+1} = '';
    L{end+1} = '## How it differs from the other scenarios';
    ref = struct('structure','S1','feature','F0','direct',true);
    d = diffs(scn, ref, CAT);
    if isempty(d)
        L{end+1} = '- This **is the reference** (joint chart, users in front of both surfaces, per-sweep covariance + fusion, direct path on).';
    else
        L{end+1} = sprintf('- **vs the reference** (S1 joint, F0 per-sweep covariance + fusion, direct on): %s.', strjoin(d, '; '));
    end
    L{end+1} = '- **Scenarios that differ in exactly one factor** (compare with these to isolate it):';
    for o = allScn
        dd = diffs(o, scn, CAT);
        if numel(dd) == 1
            L{end+1} = sprintf('  - `%s`: %s', o.name, dd{1}); %#ok<AGROW>
        end
    end
    L{end+1} = '- Older reference (no noise, old fixed 1e-6 floor): `results_aoa/coarse_3x3` (3x3) and `results_aoa/full_11x11`.';
    L{end+1} = '';
    L{end+1} = '## Files';
    L{end+1} = '- `K<k>/summary.csv`: median / mean / p90 / RMSE of |error| [deg] per method, user subset and EMS (EMS1, EMS2, pooled = both sets of errors together, joint = sqrt(e1^2+e2^2) for joint charts).';
    L{end+1} = '- `K<k>/per_user_*.csv`: per-user truth, estimate and error of the main method, with flags (front of each surface, direct path present, reachable from each surface, anchor).';
    L{end+1} = '- `K<k>/cdf.png`: CDFs per EMS, pooled and (joint charts) joint. `K<k>/map_per_ems*.png`: error maps per EMS (absolute 0-30 deg and signed +-30 deg; S3 also all users, 0-90 deg).';
    L{end+1} = '- `trend.csv`, `trend.png`: error vs codebook size (number of acquisitions).';
    L{end+1} = '';
    L{end+1} = '## Definitions';
    L{end+1} = '- **Error:** estimated minus true azimuth, wrapped to [-180, 180) deg; azimuth measured from each surface normal, positive toward +y.';
    L{end+1} = '- **Subsets:** `common` = in front of both surfaces and not an anchor in any scenario (identical users everywhere: use it to compare scenarios); `natural` = every user this scenario estimates in front of that surface; `reach` = natural with a non-zero channel to that surface; `all` (S3) = every user, also behind the surface (worst case).';
    if strcmp(scn.feature,'F0')
        L{end+1} = '- **Methods:** `typical_single` = average over the configurations of a single-configuration chart; `best_single` = best configuration (oracle, chosen on the scored users); `global_*` / `local_*` = fusion of all configurations. `*_orig` weights each chart by its quality measured against the **true** labels of all users, exactly as the original code (so it uses the unlabeled users'' ground truth); `*_free` measures the same quality against the input dissimilarities instead (no labels needed). Fusion gain = typical_single minus the fused error.';
        L{end+1} = '- Note: the original code calls this quality "TW" but it is column 2 of `drscore_pointwise`, i.e. the continuity; kept as is for comparability.';
    else
        L{end+1} = '- **Method:** `chart` = the single chart built from all acquisitions (no fusion); its improvement with K comes only from the richer combined feature.';
    end
    txt = strjoin(L, newline);
    fid = fopen(fullfile(outDir,'README.md'), 'w'); fprintf(fid, '%s\n', txt); fclose(fid);
end

function d = diffs(a, b, CAT)
% factors in which scenario a differs from b, described for a
    d = {};
    if ~strcmp(a.structure, b.structure)
        d{end+1} = sprintf('users/charts %s instead of %s', ...
            CAT.structures(strcmp({CAT.structures.id}, a.structure)).tok, ...
            CAT.structures(strcmp({CAT.structures.id}, b.structure)).tok);
    end
    if ~strcmp(a.feature, b.feature)
        d{end+1} = sprintf('feature %s instead of %s', ...
            CAT.features(strcmp({CAT.features.id}, a.feature)).tok, ...
            CAT.features(strcmp({CAT.features.id}, b.feature)).tok);
    end
    if a.direct ~= b.direct
        onoff = {'off','on'};
        d{end+1} = sprintf('direct path %s instead of %s', onoff{a.direct+1}, onoff{b.direct+1});
    end
end
