function aoa_error_map(tag, cmax)
% AOA_ERROR_MAP  Redraw the AoA error map of a finished run with your own color range.
%
%   aoa_error_map                     asks which run and the top of the color scale
%   aoa_error_map('full_11x11', 20)   errors >= 20 deg take the top color
%
% Works on any PC with the repo (no dataset, no re-simulation): it reads
% results_aoa/<tag>/per_user_errors.csv, shows the figure and overwrites
% results_aoa/<tag>/aoa_error_map.png.

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root,'src'));

    if nargin < 1 || isempty(tag)
        d = dir(fullfile(root,'results_aoa','*','per_user_errors.csv'));
        [~, runs] = cellfun(@fileparts, {d.folder}, 'UniformOutput', false);
        fprintf('Runs: %s\n', strjoin(runs, ', '));
        tag = input(sprintf('Run [%s]: ', runs{end}), 's');
        if isempty(tag), tag = runs{end}; end
    end
    if nargin < 2 || isempty(cmax)
        cmax = input('Top of the color scale [deg] [30]: ');
        if isempty(cmax), cmax = 30; end
    end

    runDir = fullfile(root,'results_aoa',tag);
    T = readtable(fullfile(runDir,'per_user_errors.csv'));
    G = readtable(fullfile(runDir,'per_user_errors_geometry.csv'), 'ReadRowNames', true);
    geo = struct('R_xy', [T.x_m T.y_m], ...
                 'EMS_position', G{{'EMS1','EMS2'},:}, ...
                 'BS_position',  G{'BS',:});

    outFile = fullfile(runDir,'aoa_error_map.png');
    plot_aoa_error_map(geo, T.joint_err_deg, outFile, cmax, true);
    fprintf('%s: %.0f%% of users below %g deg. Saved %s\n', tag, ...
            100*mean(T.joint_err_deg < cmax), cmax, outFile);
end
