%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WHICH USERS NEED THE EMSs, AND WHICH CAN EACH EMS SERVE?
%   1) users in line of sight of the BS (blue), reached only through a reflection (orange), or
%      with no BS-UE channel at all (red). H_LOS_1bounce holds LoS + single-bounce paths summed,
%      so LoS is identified from its power relative to free space: the ratio is bimodal (LoS
%      around -25 dB, reflection-only around -56 dB, gap near -36 dB) -> threshold NLOS_THR_DB
%   2,3) users WITHOUT LoS (reflection-only or no BS-UE channel) that EMS1 / EMS2 can reach (in front of it, non-zero
%        incident channel), coloured by the best uplink SNR that EMS alone gives them over a
%        fine codebook (N phase steps over -180..180 deg/element), and the angular sector they
%        span seen from that EMS (what a focused codebook would have to cover).
% Environment: NLOS_N (default 60 = 2D/lambda resolvable beams of a 30-lambda aperture)
% Outputs: results_aoa/serving/nlos_maps/ (PNG, summary.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
set(0,'DefaultFigureVisible','off');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
params = init_parameters();  params.target = 'aoa';
Data = load_dataset(params);
U = abl_users(Data, params);
N = str2double(getenv('NLOS_N'));  if isnan(N), N = 60; end
thr = str2double(getenv('NLOS_THR_DB'));  if isnan(thr), thr = -36; end
lam = 3e8/30e9;
Pd  = squeeze(sum(sum(abs(Data.H_LOS_1bounce).^2,2),3)) / numel(Data.H_LOS_1bounce(1,:,:));
dd  = sqrt(sum(([U.xy, params.UE_height*ones(U.n,1)] - U.BS_position).^2, 2));
rfs = 10*log10(Pd ./ (lam./(4*pi*dd)).^2);              % direct-channel power vs free space [dB]
los = U.has_direct & rfs > thr;                         % line of sight
refl = U.has_direct & ~los;                             % BS reaches the user only by a reflection
step = ((1:N) - N/2) * 360/N;
rho  = 10^((23 - (-92))/10);
Nb = params.Par.Tx.ArrSize;  Nr = size(Data.H_LOS_1bounce,3);
outDir = fullfile(root,'results_aoa','serving','nlos_maps');
if ~exist(outDir,'dir'), mkdir(outDir); end

best = zeros(U.n,2);  bestStep = zeros(U.n,2);
for e = 1:2
    Ho = double(Data.(sprintf('H_o_RIS%d',e)));
    X  = reshape(permute(Data.(sprintf('H_I_RIS%d',e)),[2 3 1]), params.M*params.N, []);
    for i = 1:N
        th = exp(-1j*kron(deg2rad(step(i))*(1:params.M), ones(1,params.N))).';
        H  = reshape(Ho*(th.*X), Nb, Nr, U.n);
        h1 = squeeze(H(:,1,:)); h2 = squeeze(H(:,2,:));
        g11 = sum(abs(h1).^2,1); g22 = sum(abs(h2).^2,1); g12 = sum(conj(h1).*h2,1);
        l = ((g11+g22)/2 + sqrt(((g11-g22)/2).^2 + abs(g12).^2)).';
        b = l > best(:,e);  best(b,e) = l(b);  bestStep(b,e) = step(i);
    end
end
snrE = 10*log10(rho*best);                  % best SNR through EMS e alone [dB]
nlos = ~los;                                 % no line of sight (reflection-only or no channel)
serv = nlos & U.front & U.reach & best > 0;  % NLoS users EMS e can reach

L = {sprintf('Users: %d; LoS: %d; reflection-only: %d; no BS-UE channel: %d (LoS threshold %g dB vs free space)', ...
             U.n, nnz(los), nnz(refl), nnz(~U.has_direct), thr), ...
     sprintf('Without LoS and reachable by EMS1: %d, by EMS2: %d, by both: %d, by neither: %d', ...
             nnz(serv(:,1)), nnz(serv(:,2)), nnz(all(serv,2)), nnz(nlos & ~any(serv,2))), ...
     sprintf('Codebook used for the best SNR: %d phase steps, %.0f deg/element apart', N, 360/N)};

% 1) LoS / reflection-only / no channel
f = figure('Position',[100 100 760 720]); hold on; grid on; box on; axis equal;
scatter(U.xy(los,1),  U.xy(los,2),  9, [0 0.35 0.85], 'filled', 'DisplayName', sprintf('line of sight (%d)', nnz(los)));
scatter(U.xy(refl,1), U.xy(refl,2), 9, [1 0.55 0], 'filled', 'DisplayName', sprintf('reflection only (%d)', nnz(refl)));
scatter(U.xy(~U.has_direct,1), U.xy(~U.has_direct,2), 9, [0.85 0.1 0.1], 'filled', ...
        'DisplayName', sprintf('no BS-UE channel (%d)', nnz(~U.has_direct)));
markers(U);
legend('Location','southoutside','Orientation','horizontal');
xlabel('x [m]'); ylabel('y [m]'); title('Users by BS-UE link: LoS, reflection only, none');
exportgraphics(f, fullfile(outDir,'map_los_reflection_none.png'), 'Resolution',130); close(f);
histogram_fig(rfs(U.has_direct), thr, fullfile(outDir,'los_power_histogram.png'));

% 2,3) NLoS users each EMS can serve
for e = 1:2
    s  = serv(:,e);
    az = U.az(s,e);  sn = sind(az);
    p  = prctile(az, [5 95]);  ps = prctile(sn, [5 95]);
    bw = 2/(params.M);                       % beam spacing in sin(theta) for a 60-element lambda/2 aperture
    L{end+1} = sprintf(['EMS%d: %d users without LoS, azimuth 5-95%% range %.0f..%.0f deg (sin %.2f..%.2f) -> about %d of %d resolvable beams; ' ...
                        'best SNR via EMS%d: median %.1f dB, %d users >= 0 dB; best phase steps used: %.0f..%.0f deg/element'], ...
                       e, nnz(s), p(1), p(2), ps(1), ps(2), ceil(diff(ps)/bw), 2*params.M/2, e, median(snrE(s,e)), ...
                       nnz(snrE(s,e) >= 0), min(bestStep(s,e)), max(bestStep(s,e))); %#ok<SAGROW>
    f = figure('Position',[100 100 800 700]); hold on; grid on; box on; axis equal;
    scatter(U.xy(los,1), U.xy(los,2), 6, [0.82 0.82 0.82], 'filled', 'DisplayName','line of sight (EMS not needed)');
    o = nlos & ~s;
    scatter(U.xy(o,1), U.xy(o,2), 8, [1 0.7 0.7], 'filled', 'DisplayName', sprintf('no LoS, not reachable by EMS%d', e));
    scatter(U.xy(s,1), U.xy(s,2), 14, snrE(s,e), 'filled', 'DisplayName', sprintf('no LoS, reachable by EMS%d', e));
    colormap(turbo); cb = colorbar; ylabel(cb, sprintf('best uplink SNR through EMS%d alone [dB]', e));
    clim(prctile(snrE(s,e), [2 98]));
    markers(U, e);
    legend('Location','southoutside','Orientation','vertical');
    xlabel('x [m]'); ylabel('y [m]');
    title(sprintf('EMS%d: %d users without LoS it can reach; azimuth %.0f..%.0f deg (5-95%%)', e, nnz(s), p(1), p(2)));
    exportgraphics(f, fullfile(outDir, sprintf('map_nlos_served_by_EMS%d.png', e)), 'Resolution',130); close(f);
end
fid = fopen(fullfile(outDir,'summary.txt'),'w'); fprintf(fid, '%s\n', L{:}); fclose(fid);
fprintf('%s\n', L{:});

function markers(U, e)
    if nargin < 2, e = 0; end
    for k = 1:2
        fc = [1 1 0]; if e && k ~= e, fc = [0.7 0.7 0.7]; end
        plot(U.EMS_position(k,1), U.EMS_position(k,2), 'ks', 'MarkerSize',11, 'LineWidth',1.5, ...
             'MarkerFaceColor',fc, 'HandleVisibility','off');
        text(U.EMS_position(k,1), U.EMS_position(k,2)-3.5, sprintf('EMS%d',k), 'HorizontalAlignment','center');
    end
    plot(U.BS_position(1), U.BS_position(2), 'k^', 'MarkerSize',11, 'MarkerFaceColor','r', 'HandleVisibility','off');
    text(U.BS_position(1), U.BS_position(2)+3.5, 'BS', 'HorizontalAlignment','center');
end

function histogram_fig(r, thr, file)
    f = figure('Position',[100 100 650 380]); hold on; grid on; box on;
    histogram(r, -80:2:0, 'FaceColor',[0.4 0.4 0.4]);
    xline(thr, 'r--', sprintf('LoS threshold %g dB', thr), 'LineWidth',1.5);
    xlabel('BS-UE channel power relative to free space [dB]'); ylabel('users');
    title('Bimodal: line of sight (right) vs reflection only (left)');
    exportgraphics(f, file, 'Resolution',130); close(f);
end
