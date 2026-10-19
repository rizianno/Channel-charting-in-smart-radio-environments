function state = prepare_data(Data, params)
% PREPARE_DATA
% Extracts and prepares dataset components for the EMS-based channel model.
%
% This function:
%   - selects a subset of users
%   - extracts direct and EMS-reflected channel components
%   - organizes data into a consistent state structure
%
% IMPORTANT:
%   The dataset .mat file uses legacy field names (H_I_RIS1, H_o_RIS1, etc.).
%   Internally all variables use "EMS" to match the paper terminology.
%
% Outputs:
%   state struct containing:
%       H_dir      : Direct (LOS + bounce) channel
%       H_i_EMS1   : Incident channel to EMS 1
%       H_o_EMS1   : Outgoing channel from EMS 1
%       H_i_EMS2   : Incident channel to EMS 2
%       H_o_EMS2   : Outgoing channel from EMS 2
%       R_xy       : Ground-truth user positions
%       Len        : Number of users

    %% --------------------------------------------------------------------
    % USER SUBSAMPLING
    % ---------------------------------------------------------------------
    idxUsers = 1:params.Total_samples_used;

    % AoA target: the azimuth is only defined for users in front of BOTH
    % surfaces (a surface cannot serve a user behind its wall).
    isAoA = isfield(params,'target') && strcmpi(params.target,'aoa');
    if isAoA
        emsPos = [Data.RIS1_position; Data.RIS2_position];
        [az_all, front_all] = compute_ems_azimuth( ...
            Data.R_xy(idxUsers,:), emsPos, params);
        keep = all(front_all,2);
        idxUsers = idxUsers(keep);
        az_all   = az_all(keep,:);
        if isfield(params,'max_users') && ~isempty(params.max_users)
            nKeep    = min(numel(idxUsers), params.max_users);
            idxUsers = idxUsers(1:nKeep);
            az_all   = az_all(1:nKeep,:);
        end
        state.EMS_position = emsPos;
        state.BS_position  = Data.BS_position;
    end

    state.Len = length(idxUsers);

    %% --------------------------------------------------------------------
    % DIRECT CHANNEL (LOS + SCATTERING)
    % ---------------------------------------------------------------------
    state.H_dir = Data.H_LOS_1bounce(idxUsers,:,:);

    %% --------------------------------------------------------------------
    % USER POSITIONS (GROUND TRUTH)
    % ---------------------------------------------------------------------
    state.R_xy = Data.R_xy(idxUsers,:);

    %% --------------------------------------------------------------------
    % CHART TARGET (what t-SNE anchors and all errors refer to)
    % ---------------------------------------------------------------------
    if isAoA
        state.AoA_deg = az_all;          % [Len x 2]: azimuth from EMS1, EMS2
        state.Y       = state.AoA_deg;
    else
        state.Y       = state.R_xy;      % original behaviour
    end

    %% --------------------------------------------------------------------
    % EMS CHANNEL COMPONENTS
    % ---------------------------------------------------------------------
    % NOTE:
    % The dataset field names use legacy notation (H_I_RIS1, H_o_RIS1, etc.).
    % We map them to EMS naming here and nowhere else in the codebase.

    % EMS 1
    state.H_i_EMS1 = Data.H_I_RIS1(idxUsers,:,:);   % incident
    state.H_o_EMS1 = Data.H_o_RIS1;                 % outgoing

    % EMS 2
    state.H_i_EMS2 = Data.H_I_RIS2(idxUsers,:,:);   % incident
    state.H_o_EMS2 = Data.H_o_RIS2;                 % outgoing

    %% --------------------------------------------------------------------
    % SANITY CHECKS (CRITICAL FOR DEBUGGING)
    % ---------------------------------------------------------------------
    assert(size(state.H_dir,1) == state.Len, 'Mismatch in H_dir size');
    assert(size(state.R_xy,1)  == state.Len, 'Mismatch in R_xy size');

end