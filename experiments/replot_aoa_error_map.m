%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REDRAW aoa_error_map.png OF FINISHED RUNS (no re-simulation)
%   RUN_TAGS="coarse_3x3,full_11x11"  (default: every folder in results_aoa)
%   AOA_MAP_CMAX="40"                 top of the color scale [deg]
% Uses results_aoa/<tag>/results.mat; runs saved before 'geo' was stored
% get the user positions back from the dataset with the run's own params.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
set(0,'DefaultFigureVisible','off');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));

cmax = str2double(getenv('AOA_MAP_CMAX'));
if isnan(cmax), cmax = 40; end

tags = strtrim(strsplit(getenv('RUN_TAGS'), ','));
if isempty(tags{1})
    d = dir(fullfile(root,'results_aoa','*'));
    tags = {d([d.isdir] & ~startsWith({d.name},'.')).name};
end

Data = [];
for k = 1:numel(tags)
    runDir = fullfile(root,'results_aoa',tags{k});
    S = load(fullfile(runDir,'results.mat'));
    if isfield(S,'geo')
        geo = S.geo;
    else
        if isempty(Data), Data = load_dataset(S.params); end
        state = prepare_data(Data, S.params);
        unl   = S.params.supervised_Len+1:state.Len;
        assert(isequal(state.Y, S.Y), 'User set differs from the saved run.');
        geo = struct('R_xy', state.R_xy(unl,:), 'EMS_position', state.EMS_position, ...
                     'BS_position', state.BS_position);
    end
    plot_aoa_error_map(geo, S.err.local, fullfile(runDir,'aoa_error_map.png'), cmax);
    fprintf('redrawn: %s\n', fullfile(runDir,'aoa_error_map.png'));
end
