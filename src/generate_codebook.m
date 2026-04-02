function codebook = generate_codebook(params)
% GENERATE_CODEBOOK
% Constructs the EMS phase-gradient codebook.
%
% Each codebook entry corresponds to a pair of horizontal phase gradients:
%   [phi1, phi2] → applied to EMS 1 and EMS 2 respectively
%
% The full codebook is obtained by exhaustively combining:
%   params.phi_options × params.phi_options
%
% Output:
%   codebook.all         : [nCodebooks × 2] matrix of phase pairs
%   codebook.nCodebooks  : total number of configurations

    %% --------------------------------------------------------------------
    % Generate all combinations of phase gradients (EMS1, EMS2)
    % ---------------------------------------------------------------------
    [phi1_grid, phi2_grid] = meshgrid(params.phi_options, params.phi_options);

    codebook.all = [phi1_grid(:), phi2_grid(:)];
    codebook.nCodebooks = size(codebook.all,1);

    %% --------------------------------------------------------------------
    % Sanity check
    % ---------------------------------------------------------------------
    expected = params.nPhasesPerEMS^2;
    assert(codebook.nCodebooks == expected, ...
        'Codebook size mismatch');

end