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