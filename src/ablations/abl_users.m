function U = abl_users(Data, params)
% ABL_USERS  All dataset users with AoA labels, visibility flags and the anchor sets.
%   az(:,e)      : azimuth seen from EMS e [deg], from its normal, positive toward +y;
%                  for users behind the surface it lies beyond +-90 (geometric, up to +-180)
%   front(:,e)   : user in front of EMS e;  both = in front of both
%   has_direct   : direct channel not identically zero
%   reach(:,e)   : incident channel UE -> EMS e not identically zero
%   anchorsJoint : 15% of the users in front of both (dataset order), as the reference runs
%   anchorsEMS{e}: 15% of the users in front of EMS e (dataset order)
    U.n  = size(Data.R_xy,1);
    U.xy = Data.R_xy;
    U.EMS_position = [Data.RIS1_position; Data.RIS2_position];
    U.BS_position  = Data.BS_position;
    [U.az, U.front] = compute_ems_azimuth(U.xy, U.EMS_position, params);
    U.both = all(U.front, 2);

    U.has_direct = squeeze(sum(sum(abs(Data.H_LOS_1bounce).^2, 2), 3)) > 0;
    U.reach = [squeeze(sum(sum(abs(Data.H_I_RIS1).^2, 2), 3)) > 0, ...
               squeeze(sum(sum(abs(Data.H_I_RIS2).^2, 2), 3)) > 0];

    sup  = params.Supervision_percentage;
    idxB = find(U.both);
    U.anchorsJoint = idxB(1:round(sup*numel(idxB)));
    for e = 1:2
        idxF = find(U.front(:,e));
        U.anchorsEMS{e} = idxF(1:round(sup*numel(idxF)));
    end
    U.anyAnchor = false(U.n,1);
    U.anyAnchor([U.anchorsJoint; U.anchorsEMS{1}; U.anchorsEMS{2}]) = true;
end
