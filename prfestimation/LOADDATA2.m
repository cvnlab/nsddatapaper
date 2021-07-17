% Expects hh and cc to be defined.
% Creates <somedata> in the workspace.

% h5read parameters
arg1 = [(cc-1)*numinchunk+1 1];
if cc==ceil(datasizes(hh,1)/numinchunk)
  arg2 = [datasizes(hh,1)-(cc-1)*numinchunk Inf];
else
  arg2 = [numinchunk Inf];
end

% load the data
tic;
somedata = h5read(sprintf('~/ext/figurefiles/nsd/datab3nativesurface_subj%02d_%s_betas.hdf5', ...
                          subjix,datanames{hh}),'/betas',arg1,arg2);
toc;
