%% Based on files from alex white, refine.

zz = 8;

subjid = sprintf('subj%02d',zz);   % which subject
cmap   = jet(256);   % colormap for ROIs
rng    = [0 7];      % should be [0 N] where N is the max ROI index
roilabels = {'OWFA' 'VWFA-1' 'VWFA-1b' 'VWFA-2' 'MFS-WORDS' 'ATL-words' 'unlabeled'};  % 1 x N cell vector of strings
mgznames = {'flocwordtval'};%{'curvature'};           % quantities of interest (1 x Q)
crngs = {[-10 10]};%{[-1 1]};          % ranges for the quantities (1 x Q)
cmaps = {cmapsign4(256)};%{copper};   % colormaps for the quantities (1 x Q)
  %roivals = [];                                                                    % start with a blank slate?
  %roivals = cvnloadmgz(sprintf('%s/%s/label/?h.visualsulc.mgz',cvnpath('freesurfer'),subjid));  % load in an existing file?
roivals = cvnloadmgz(sprintf('/Users/kendrick/Dropbox/nsdwords/floc-words/ROIs/labeledWordROIs/subj%02d/?h.floc-words.mgz',zz));

% do it
cvndefinerois;

% final things to do:
% -throw away VWFA1-b and unlabeled
% -erase speckles
% -fill holes
% -save into official directory freesurfer/subj01/label/lh.floc-words.mgz
% -note that the final rng is 0 5!

%% continue

% make .ctab file with final names!

% NOTE: Apr 19 2020, regen from try2words

% whittle based on t>0
hemis = {'lh' 'rh'};
for zz=1:8
  mkdirquiet(sprintf('~/Desktop/%d/',zz));
  for hemi=1:2
    data0 = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.flocwordtval.mgz',zz,hemis{hemi}));
    file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.floc-words.mgz',zz,hemis{hemi});
    roi0 = cvnloadmgz(file0);
    for roi=1:5
      roi0(roi0==roi) = (data0(roi0==roi)>0) * roi;
    end
    copyfile(file0,sprintf('~/Desktop/%d/',zz));
    nsd_savemgz(roi0(:),file0,sprintf('~/nsd/nsddata/freesurfer/subj%02d/',zz));
  end
end

%% create ctab file

for subjix=1:8
  copyfile(sprintf('~/Dropbox/KKTEMP/%s.mgz.ctab','floc-words'), ...
           sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/',subjix));
end

%%%%% proceed to analysis_drawrois_facesandwords.m
