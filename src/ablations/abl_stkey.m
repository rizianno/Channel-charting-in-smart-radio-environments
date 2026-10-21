function k = abl_stkey(direct, phi)
% ABL_STKEY  Name of one acquisition: direct path flag + phase gradient of each EMS (or off).
    p = cell(1,2);
    for e = 1:2
        if isnan(phi(e)), p{e} = 'off'; else, p{e} = sprintf('%.3f', phi(e)); end
    end
    k = sprintf('d%d_%s_%s', direct, p{1}, p{2});
end
