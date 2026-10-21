function aoa_error_map_per_ems(tag, cmax1, cmax2, mode)
% AOA_ERROR_MAP_PER_EMS  AoA error seen from EMS1 and from EMS2, as two maps
% side by side, each with its own colorbar and range.
%
%   aoa_error_map_per_ems                                   asks for everything
%   aoa_error_map_per_ems('full_11x11', 20)                 same range (20 deg) for both
%   aoa_error_map_per_ems('full_11x11', 20, 15)             EMS1 0-20, EMS2 0-15 deg
%   aoa_error_map_per_ems('full_11x11', 20, 20, 'signed')   estimate - true, -20..20 deg
%                                                           (blue under-, red over-estimate)
%
% Reads results_aoa/<tag>/per_user_errors.csv (no dataset, no re-simulation),
% shows the figure and saves results_aoa/<tag>/aoa_error_map_per_ems[_signed].png.

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root,'src'));

    if nargin < 1 || isempty(tag)
        d = dir(fullfile(root,'results_aoa','*','per_user_errors.csv'));
        [~, runs] = cellfun(@fileparts, {d.folder}, 'UniformOutput', false);
        fprintf('Runs: %s\n', strjoin(runs, ', '));
        tag = input(sprintf('Run [%s]: ', runs{end}), 's');
        if isempty(tag), tag = runs{end}; end
    end
    if nargin < 2 || isempty(cmax1)
        cmax1 = input('Top of the color scale for EMS1 [deg] [30]: ');
        if isempty(cmax1), cmax1 = 30; end
        if nargin < 3
            cmax2 = input(sprintf('Top of the color scale for EMS2 [deg] [%g]: ', cmax1));
        end
    end
    if nargin < 3 || isempty(cmax2), cmax2 = cmax1; end
    if nargin < 4 || isempty(mode), mode = 'abs'; end

    runDir = fullfile(root,'results_aoa',tag);
    T = readtable(fullfile(runDir,'per_user_errors.csv'));
    G = readtable(fullfile(runDir,'per_user_errors_geometry.csv'), 'ReadRowNames', true);
    geo = struct('R_xy', [T.x_m T.y_m], ...
                 'EMS_position', G{{'EMS1','EMS2'},:}, ...
                 'BS_position',  G{'BS',:});
    Y_true = [T.az1_true_deg T.az2_true_deg];
    Y_est  = [T.az1_est_deg  T.az2_est_deg];

    suffix = '';
    if strcmpi(mode,'signed'), suffix = '_signed'; end
    outFile = fullfile(runDir, ['aoa_error_map_per_ems' suffix '.png']);
    plot_aoa_error_map_per_ems(geo, Y_true, Y_est, outFile, [cmax1 cmax2], mode, true);

    e = abs(Y_est - Y_true);
    fprintf('%s: within range -> EMS1 %.0f%% (<%g deg), EMS2 %.0f%% (<%g deg). Saved %s\n', ...
            tag, 100*mean(e(:,1) < cmax1), cmax1, 100*mean(e(:,2) < cmax2), cmax2, outFile);
end
