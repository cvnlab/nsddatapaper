function nsd_et_savefig(h, filename, outpath)
% ----------------------------------------------------------------------- %
%  save figures in favorite formats. MN, July 2021
% ----------------------------------------------------------------------- %
saveas(h,fullfile(outpath, strcat(filename, '.jpg')));
saveas(h,fullfile(outpath, strcat(filename, '.svg')));
end