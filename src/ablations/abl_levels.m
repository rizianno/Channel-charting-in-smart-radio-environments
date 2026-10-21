function lev = abl_levels(CAT, scn, K)
% ABL_LEVELS  Phase steps per EMS [deg/element] of scenario scn at codebook size K x K.
    cb = 'narrow';
    if isfield(scn, 'codebook'), cb = scn.codebook; end
    lev = CAT.codebooks(strcmp({CAT.codebooks.id}, cb)).levels{K};
end
