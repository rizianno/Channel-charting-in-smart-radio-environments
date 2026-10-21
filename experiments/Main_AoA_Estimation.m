%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AoA ESTIMATION WITH FEW EMS SWEEPS
%
% Same pipeline as Main_Multi_Slot_Fusion.m (per-codebook semi-supervised
% t-SNE on Log-Euclidean covariance dissimilarities, then global / local
% fusion), but the chart target is the azimuth of each UE as seen from
% EMS1 and EMS2 instead of its (x,y) position. Azimuths are measured from
% each surface normal, in (-90, 90) deg, so there is no wrap-around.
%
% Optional environment variables (quick vs full runs, no code edits):
%   PHI_OPTIONS_DEG  phase gradients per EMS, e.g. "0:20:40" or "0,40"
%                    (default: params.phi_options from init_parameters)
%   AOA_MAX_USERS    cap on the number of users, e.g. "400" (smoke tests)
%   RUN_TAG          name of the output folder, e.g. "coarse_3x3"
%
% Outputs in results_aoa/<RUN_TAG>/ (light files, meant to be synced back):
%   summary.csv, per_codebook.csv, notes.txt, *.png
%   results.mat  (full arrays; *.mat is git-ignored)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
set(0,'DefaultFigureVisible','off');

%% ------------------------------------------------------------------------
% INIT
%% ------------------------------------------------------------------------
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));

params = init_parameters();
params.target = 'aoa';
rng(params.random_seed);

s = getenv('PHI_OPTIONS_DEG');
if ~isempty(s)
    params.phi_options   = deg2rad(parse_numbers(s));
    params.nPhasesPerEMS = numel(params.phi_options);
end
s = getenv('AOA_MAX_USERS');
if ~isempty(s), params.max_users = str2double(s); end

tag = getenv('RUN_TAG');
if isempty(tag)
    tag = sprintf('aoa_%dx%d_%s', params.nPhasesPerEMS, params.nPhasesPerEMS, ...
                  char(datetime('now','Format','yyyyMMdd_HHmm')));
end
outDir = fullfile(root,'results_aoa',tag);
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ------------------------------------------------------------------------
% DATA (only users in front of both surfaces have a defined AoA)
%% ------------------------------------------------------------------------
Data  = load_dataset(params);
state = prepare_data(Data, params);
clear Data

params.supervised_Len = round(state.Len * params.Supervision_percentage);
lab = 1:params.supervised_Len;
unl = params.supervised_Len+1:state.Len;

fprintf('Users with a defined AoA (in front of both EMSs): %d\n', state.Len);
fprintf('Labeled anchors: %d, evaluated users: %d\n', numel(lab), numel(unl));

%% ------------------------------------------------------------------------
% CODEBOOK
%% ------------------------------------------------------------------------
codebook = generate_codebook(params);
fprintf('EMS sweeps (codebook configurations): %d  [phi per EMS: %s deg]\n', ...
    codebook.nCodebooks, mat2str(round(rad2deg(params.phi_options),2)));

%% ------------------------------------------------------------------------
% PHASE 1 + 2: same functions as the position experiment
%% ------------------------------------------------------------------------
t0 = tic;
results_phase1 = run_channel_charting(state, params, codebook);
results_global = global_weighting(results_phase1, state, params);
results_local  = local_weighting(results_phase1, state, params);
best           = evaluate_best_codebook(results_phase1, state, params);
metrics        = compute_all_metrics(results_global, results_local, best, state, params);
runtime_min    = toc(t0)/60;

%% ------------------------------------------------------------------------
% PER-EMS ANGLE ERRORS
%% ------------------------------------------------------------------------
Y = state.Y;
est.best   = best.xy;
est.global = results_global.xy;
est.local  = results_local.xy;
% Reference: predict every user's AoA as the mean AoA of the labeled users
est.anchor_mean = repmat(mean(Y(lab,:),1), state.Len, 1);

methods = {'best_codebook','global_fusion','local_fusion','anchor_mean_baseline'};
fields  = {'best','global','local','anchor_mean'};
TW = {metrics.TW_best, metrics.TW_global, metrics.TW_local, NaN};
CT = {metrics.CT_best, metrics.CT_global, metrics.CT_local, NaN};

rows = cell(numel(methods),1);
err  = struct();
for k = 1:numel(methods)
    E = est.(fields{k})(unl,:) - Y(unl,:);
    e1 = abs(E(:,1)); e2 = abs(E(:,2)); ej = sqrt(sum(E.^2,2));
    err.(fields{k}) = ej;
    rows{k} = {methods{k}, mean(ej), median(ej), prctile(ej,90), ...
               mean(e1), median(e1), prctile(e1,90), ...
               mean(e2), median(e2), prctile(e2,90), ...
               mean(TW{k}), mean(CT{k})};
end
T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'method','joint_mean_deg','joint_median_deg','joint_p90_deg', ...
     'EMS1_mean_deg','EMS1_median_deg','EMS1_p90_deg', ...
     'EMS2_mean_deg','EMS2_median_deg','EMS2_p90_deg','TW_mean','CT_mean'});
writetable(T, fullfile(outDir,'summary.csv'));

Tcb = table((1:codebook.nCodebooks).', rad2deg(codebook.all(:,1)), rad2deg(codebook.all(:,2)), ...
    results_phase1.ErrorPerCodebook, results_phase1.TW_PerCodebook_mean, ...
    results_phase1.CT_PerCodebook_mean, results_global.weights(:), ...
    'VariableNames', {'codebook','phi1_deg','phi2_deg','joint_mean_err_deg', ...
                      'TW_mean','CT_mean','global_weight'});
writetable(Tcb, fullfile(outDir,'per_codebook.csv'));

%% ------------------------------------------------------------------------
% FIGURES
%% ------------------------------------------------------------------------
col = struct('best',[0 0.45 0.74],'global',[0.85 0.33 0.10], ...
             'local',[0.47 0.67 0.19],'anchor_mean',[0.5 0.5 0.5]);

% 1) CDF of the joint angular error
f = figure('Position',[100 100 700 500]); hold on; grid on; box on;
for cb = 1:codebook.nCodebooks
    ecb = results_phase1.ErrorPerUser(unl,cb);
    plot(sort(ecb), linspace(0,1,numel(ecb)), 'Color',[0.85 0.85 0.85], ...
        'HandleVisibility','off');
end
for k = 1:numel(fields)
    ej = err.(fields{k});
    plot(sort(ej), linspace(0,1,numel(ej)), 'LineWidth',2, ...
        'Color', col.(fields{k}), 'DisplayName', strrep(methods{k},'_',' '));
end
xlabel('Joint AoA error  \surd(\Delta\theta_1^2+\Delta\theta_2^2)  [deg]');
ylabel('Empirical CDF'); legend('Location','southeast');
title(sprintf('AoA error, %d EMS sweeps, %d users', codebook.nCodebooks, numel(unl)));
exportgraphics(f, fullfile(outDir,'aoa_error_cdf.png'), 'Resolution',150); close(f);

% 2) Estimated vs true azimuth, per EMS (global fusion)
f = figure('Position',[100 100 1000 450]);
for e = 1:2
    subplot(1,2,e); hold on; grid on; box on; axis equal;
    scatter(Y(unl,e), est.global(unl,e), 8, 'filled', 'MarkerFaceAlpha',0.5);
    plot([-90 90],[-90 90],'k--');
    xlim([-90 90]); ylim([-90 90]);
    xlabel(sprintf('True azimuth from EMS%d [deg]', e));
    ylabel(sprintf('Estimated azimuth from EMS%d [deg]', e));
    title(sprintf('EMS%d (global fusion)', e));
end
exportgraphics(f, fullfile(outDir,'aoa_true_vs_est.png'), 'Resolution',150); close(f);

% 3) Where the errors are (local fusion), with EMS and BS positions
%    color scale 0-30 deg (larger errors saturate) so the typical errors are visible
geo = struct('R_xy', state.R_xy(unl,:), 'EMS_position', state.EMS_position, ...
             'BS_position', state.BS_position);
plot_aoa_error_map(geo, err.local, fullfile(outDir,'aoa_error_map.png'), 30);
% light per-user CSV: lets aoa_error_map.m redraw it with any color range
write_per_user_errors(geo, Y(unl,:), est.local(unl,:), fullfile(outDir,'per_user_errors.csv'));
% same, per EMS (two maps, own colorbars)
plot_aoa_error_map_per_ems(geo, Y(unl,:), est.local(unl,:), ...
    fullfile(outDir,'aoa_error_map_per_ems.png'), 30, 'abs');

% 4) Error per codebook configuration
f = figure('Position',[100 100 700 400]); hold on; grid on; box on;
bar(results_phase1.ErrorPerCodebook, 'FaceColor',[0.7 0.7 0.7]);
yline(mean(err.global),'--','global fusion','Color',col.global,'LineWidth',1.5);
yline(mean(err.local), '--','local fusion', 'Color',col.local, 'LineWidth',1.5);
xlabel('Codebook configuration'); ylabel('Mean joint AoA error [deg]');
title('AoA error per EMS configuration');
exportgraphics(f, fullfile(outDir,'aoa_error_per_codebook.png'), 'Resolution',150); close(f);

%% ------------------------------------------------------------------------
% NOTES + FULL RESULTS
%% ------------------------------------------------------------------------
[~, gitHash] = system(sprintf('git -C "%s" describe --always --dirty', root));
fid = fopen(fullfile(outDir,'notes.txt'),'w');
fprintf(fid, 'Run tag         : %s\n', tag);
fprintf(fid, 'Date            : %s\n', char(datetime('now')));
fprintf(fid, 'Code commit     : %s\n', strtrim(gitHash));
fprintf(fid, 'Users (AoA def.): %d  (in front of both EMSs)\n', state.Len);
fprintf(fid, 'Labeled anchors : %d (%.0f%%), evaluated: %d\n', numel(lab), ...
        100*params.Supervision_percentage, numel(unl));
fprintf(fid, 'EMS sweeps      : %d  (phi per EMS [deg]: %s)\n', codebook.nCodebooks, ...
        mat2str(round(rad2deg(params.phi_options),2)));
fprintf(fid, 'Fusion weighting: %s (%s)\n', params.WEIGHTING_METRIC, params.SCENARIO_AGGREGATION);
fprintf(fid, 'Runtime         : %.1f min\n\n', runtime_min);
fprintf(fid, 'Azimuth: angle of the UE seen from each EMS, 0 deg = surface normal,\n');
fprintf(fid, '         positive toward +y (codebook-gradient axis); range (-90,90).\n');
fprintf(fid, 'Joint error: sqrt(dtheta1^2 + dtheta2^2) per user, unlabeled users only.\n');
fprintf(fid, 'best_codebook is an oracle (picked by its error on the evaluated users).\n');
fprintf(fid, 'anchor_mean_baseline: every user predicted at the mean labeled AoA.\n\n');
fprintf(fid, '%s\n', evalc('disp(T)'));
fclose(fid);

save(fullfile(outDir,'results.mat'), 'T','Tcb','est','err','Y','params','geo', ...
     'results_phase1','results_global','results_local','best','metrics','-v7.3');

disp(T);
fprintf('\nOutputs written to %s\n=== DONE ===\n', outDir);


function v = parse_numbers(s)
% Parse "0:20:40" or "0,20,40" (digits and : , . - and spaces only).
    s = strtrim(s);
    assert(~isempty(regexp(s,'^[0-9:,\.\-\s]+$','once')), ...
        'PHI_OPTIONS_DEG may only contain numbers, ":", "," and "-".');
    v = str2num(['[' s ']']); %#ok<ST2NM>
end
