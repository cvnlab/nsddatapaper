% Creates <alldata> in the workspace.
%
% This may be too expensive now that we are doing func1pt8mm...  we load only 1:3

% load the data
tic;
alldata = {};
for hh=1:3   %length(datanames)
  alldata{hh} = h5read(sprintf('~/ext/figurefiles/nsd/datab3nativesurface_subj%02d_%s_betas.hdf5',subjix,datanames{hh}),'/betas');
end
toc;
