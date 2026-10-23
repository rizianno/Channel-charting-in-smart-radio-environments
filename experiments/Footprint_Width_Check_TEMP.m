%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMPORARY DIAGNOSTIC: useful footprint of one training beam for the chart.
% For every cached single-entry chart (nested run, both EMSs on the same entry), the AoA error of
% each scored user vs its distance from the beam direction, |sin(theta_user) - sin(theta_beam)|,
% in units of the beamwidth (0.886 * 2 / M in sin). The footprint w = the distance up to which
% the chart still places users well; it sets the beam spacing of the PDF-informed training.
% Output: printed + results_aoa/serving_nlos_size_informed_TEMP/footprint_width.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
params = init_parameters();  params.target = 'aoa';
Data = load_dataset(params);  U = abl_users(Data, params);
nlos = find(~U.has_direct);
reachE = [any(Data.H_I_RIS1(nlos,:,:) ~= 0, [2 3]), any(Data.H_I_RIS2(nlos,:,:) ~= 0, [2 3])] & U.front(nlos,:);
azBS = compute_ems_azimuth(U.BS_position(1:2), U.EMS_position, params);  spec = -azBS;
edges = [0 0.5 1 2 3 4 6 8 12 16 24 32 inf];                      % distance bins [beamwidths]
out = {sprintf('AoA error of single-entry charts vs distance of the user from the beam (sin domain, in beamwidths)'), ...
       sprintf('bins: %s', mat2str(edges))};
for tag = {'M60_raytraced','M60','M120','M180','M240'}
    M = sscanf(tag{1}, 'M%d');  wu = 0.886*2/M;  N = ceil(2/wu);  N = N + mod(N,2);
    step = ((1:N) - N/2) * 360/N;
    f = [1/2 1];  for lev = 2:12, f = [f, (1:2:2^lev-1)/2^lev]; if numel(unique(round(f*N))) >= N, break; end, end %#ok<AGROW>
    ord = unique(min(max(round(f*N),1),N), 'stable');  ord = ord(1:10);
    cf = fullfile(root,'cache_ablations', sprintf('nested_est_G60dB_%s.mat', tag{1}));
    if ~isfile(cf), out{end+1} = sprintf('%s: no cache', tag{1}); continue; end %#ok<SAGROW>
    C = load(cf);  R = containers.Map(C.keys, C.vals);
    for e = 1:2
        m = find(reachE(:,e));  nA = round(params.Supervision_percentage * numel(m));
        az = U.az(nlos(m), e);  unl = (nA+1:numel(m)).';
        D = [];  E = [];
        for i = ord
            o = R(sprintf('est_G60dB|%s|nest|e%d|%d', tag{1}, e, i));
            ub = mod(sind(spec(e)) + step(i)/180 + 1, 2) - 1;          % beam direction (sin)
            dd = abs(mod(sind(az(unl)) - ub + 1, 2) - 1) / wu;          % distance [beamwidths]
            D = [D; dd]; E = [E; abs(mod(o.est(unl) - az(unl) + 180, 360) - 180)]; %#ok<AGROW>
        end
        [~, ~, b] = histcounts(D, edges);
        med = accumarray(b, E, [numel(edges)-1 1], @median, NaN);
        p10 = accumarray(b, E < 10, [numel(edges)-1 1], @mean, NaN);
        cnt = accumarray(b, 1, [numel(edges)-1 1]);
        out{end+1} = sprintf('%s EMS%d (beamwidth %.4f in sin): median |err| per bin %s deg', tag{1}, e, wu, mat2str(round(med.',1))); %#ok<SAGROW>
        out{end+1} = sprintf('%s EMS%d: share |err| < 10 deg per bin %s; users per bin %s', tag{1}, e, mat2str(round(100*p10.')), mat2str(cnt.')); %#ok<SAGROW>
    end
end
fprintf('%s\n', out{:});
od = fullfile(root,'results_aoa','serving_nlos_size_informed_TEMP');  if ~exist(od,'dir'), mkdir(od); end
fid = fopen(fullfile(od,'footprint_width.txt'),'w');  fprintf(fid, '%s\n', out{:});  fclose(fid);
