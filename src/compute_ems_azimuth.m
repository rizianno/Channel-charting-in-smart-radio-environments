function [az_deg, in_front, dircos_m] = compute_ems_azimuth(R_xy, ems_pos, params)
% COMPUTE_EMS_AZIMUTH
% Azimuth of each user as seen from each EMS, in the EMS's own frame:
%   0 deg  = along the EMS outward normal,
%   > 0    = toward the EMS codebook-gradient axis (+m).
% Users in front of a surface get azimuths in (-90, 90) deg, so the label
% never wraps around. Users behind a surface cannot be served by it and
% are flagged with in_front = false.
%
% Inputs:
%   R_xy    : [nU x 2] user (x,y) positions [m]
%   ems_pos : [nE x 3] EMS centre positions [m]
%   params  : UE_height, EMS_normal [nE x 3], EMS_axis_m [nE x 3]
%
% Outputs:
%   az_deg   : [nU x nE] azimuth [deg]
%   in_front : [nU x nE] logical, user on the illuminated side of the EMS
%   dircos_m : [nU x nE] direction cosine along the gradient axis, i.e. the
%              quantity the phase-gradient codebook acts on (for a later
%              lookup table)

    nU = size(R_xy,1);
    nE = size(ems_pos,1);
    P  = [R_xy, params.UE_height*ones(nU,1)];

    az_deg   = zeros(nU,nE);
    in_front = false(nU,nE);
    dircos_m = zeros(nU,nE);

    for e = 1:nE
        d  = P - ems_pos(e,:);                 % EMS -> user
        dn = d * params.EMS_normal(e,:).';     % component along the normal
        dm = d * params.EMS_axis_m(e,:).';     % component along the gradient axis

        az_deg(:,e)   = atan2d(dm, dn);
        in_front(:,e) = dn > 0;
        dircos_m(:,e) = dm ./ vecnorm(d,2,2);
    end

end
