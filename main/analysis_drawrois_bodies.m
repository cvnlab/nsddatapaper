%%%%%% THIS IS JUST VERY SIMILAR TO PLACES, should write more general code




setenv('SUBJECTS_DIR','/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer');

%% %%%%% transfer all current ROIs to fsaverage

files = {'floc-bodies'};
hemis = {'lh' 'rh'};

vals = zeros(163842,2,8,1);  % V x hemis x 8 x quantities
for ff=1:length(files)
  for subjix=1:8
    for hh=1:length(hemis)
      a5 = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.%s.mgz',subjix,hemis{hh},files{ff}));
      if ~isempty(a5)
        vals(:,hh,subjix,ff) = nsd_mapdata(subjix,sprintf('%s.white',hemis{hh}),'fsaverage',squeeze(a5));
      end
    end
  end
end

%% %%%%% map the group average to the individual subject

subjix = 1;
groupavg = {};
for hh=1:2
  mx = max(vals(:));
  for vv=1:mx
    groupavg{vv}{hh} = nsd_mapdata(subjix,'fsaverage',sprintf('%s.white',hemis{hh}),mean(vals(:,hh,:,1)==vv,3));
  end
end

%% %%%%% DEFINE THE ROIS

subjid = sprintf('subj%02d',subjix);   % which subject
cmap   = jet(256);   % colormap for ROIs
rng    = [0 4];      % should be [0 N] where N is the max ROI index
roilabels = {'EBA' 'FBA-1' 'FBA-2' 'mTL-bodies'};  % 1 x N cell vector of strings
mgznames = {'flocbodiestval' 'flocbodiesanglemetric' 'corticalsulc' 'Kastner2015'};%{'curvature'};           % quantities of interest (1 x Q)
mgznames = [mgznames groupavg];
crngs = {[-10 10] [-90 90] [0 28] [0 25] [0 1] [0 1] [0 1] [0 1]};%{[-1 1]};          % ranges for the quantities (1 x Q)
cmaps = {cmapsign4(256) cmapsign4(256) jet(256) jet(256) jet jet jet jet};%{copper};   % colormaps for the quantities (1 x Q)
threshs = {[] [] 0.5 0.5 1/16 1/16 1/16 1/16};
roivals = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/%s/label/?h.floc-bodies.mgz',subjid));  % load in an existing file?

% do it
cvndefinerois;

% use sphere-occip, occipA2 (more lateral view)

%% %%%%% create ctab file

for subjix=1:8
  copyfile(sprintf('~/Dropbox/KKTEMP/%s.mgz.ctab','floc-bodies'), ...
           sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/',subjix));
end

%% %%%%% whittle the ROIs down

% whittle based on t>0
hemis = {'lh' 'rh'};
for zz=1:8
  mkdirquiet(sprintf('~/Desktop/%d/',zz));
  for hemi=1:2
    data0 = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.flocbodiestval.mgz',zz,hemis{hemi}));
    file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.floc-bodies.mgz',zz,hemis{hemi});
    roi0 = cvnloadmgz(file0);
    for roi=1:4
      roi0(roi0==roi) = (data0(roi0==roi)>0) * roi;
    end
    copyfile(file0,sprintf('~/Desktop/%d/',zz));
    nsd_savemgz(roi0(:),file0,sprintf('~/nsd/nsddata/freesurfer/subj%02d/',zz));
  end
end

%% %%%%% visualize results (inspections and quality assessment)

outputdir = '~/Dropbox/KKTEMP/floc-bodies';
mkdirquiet(outputdir);

% show native inflated: corticalsulc, ROIs (with border), tvalues (with border)
todo = {{'floc-bodies' 'flocbodiestval'}};
views = [4 5 6 7];
viewnames = {'inflated-parietal' 'inflated-medial' 'inflated-lateral' 'inflated-medial-ventral'};  %  %{'medial-ventral' 'parietal' 'medial'};
for zz=1:length(todo)
  strs = todo{zz};

  for subjix=1:8

    % load corticalsulc
    file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/*.corticalsulc.mgz',subjix);
    sulc0 = cvnloadmgz(file0);

    % load ROI
    file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/*.%s.mgz',subjix,strs{1});
    roi0 = cvnloadmgz(file0);
    
    % load tval
    file0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/*.%s.mgz',subjix,strs{2});
    tval0 = cvnloadmgz(file0);

    % write map
    for vv=1:length(views)

      % write corticalsulc
      cvnlookup(sprintf('subj%02d',subjix),views(vv),sulc0,[0 28],jet(256),0.5,[],0);
      imwrite(rgbimg,sprintf('%s/subj%02d-%s_corticalsulc.png',outputdir,subjix,viewnames{vv}));

      % write curvature
      cvnlookup(sprintf('subj%02d',subjix),views(vv),sulc0,[0 28],jet(256),1000,[],0);
      imwrite(rgbimg,sprintf('%s/subj%02d-%s_curvature.png',outputdir,subjix,viewnames{vv}));

      % write ROIs
      extraopts = {'roiname',strs{1},'roiwidth',0,'drawroinames',true};
      cvnlookup(sprintf('subj%02d',subjix),views(vv),roi0,[0 3],jet,0.5,[],0,extraopts);
      imwrite(rgbimg,sprintf('%s/subj%02d-%s_%s.png',outputdir,subjix,viewnames{vv},strs{1}));

      % write tvalue
      extraopts = {'roiname',strs{1},'roiwidth',1,'roicolor','w','drawroinames',true};
      cvnlookup(sprintf('subj%02d',subjix),views(vv),tval0,[-10 10],cmapsign4(256),[],[],0,extraopts);
      imwrite(rgbimg,sprintf('%s/subj%02d-%s_%s.png',outputdir,subjix,viewnames{vv},strs{2}));

    end

  end

end

% show fsaverage flat: tvalues (with border name)
hemis = {'lh' 'rh'};
vals = [];
for subjix=1:8

  for hh=1:2
    tval0 = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.%s.mgz',subjix,hemis{hh},'flocbodiestval'));
    roi0 = cvnloadmgz(sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/%s.%s.mgz',subjix,hemis{hh},'floc-bodies'));
    vals(:,:,hh,subjix) = nsd_mapdata(subjix,sprintf('%s.white',hemis{hh}),'fsaverage',[tval0(:) roi0(:)]);
  end

  extraopts = {'roimask',{vflatten(vals(:,2,:,subjix))==1 vflatten(vals(:,2,:,subjix))==2 vflatten(vals(:,2,:,subjix))==3 vflatten(vals(:,2,:,subjix))==4}, ...
               'roiwidth',1,'roicolor','w','drawroinames',{'EBA' 'FBA-1' 'FBA-2' 'mTL-bodies'}};
  cvnlookup('fsaverage',13,vflatten(vals(:,1,:,subjix)),[-10 10],cmapsign4(256),[],[],0,extraopts);
  imwrite(rgbimg,sprintf('%s/fsaverage_subj%02d_%s.png',outputdir,subjix,'floc-bodies'));

end

% show fsaverage flat: probability maps
roilabels = {'EBA' 'FBA-1' 'FBA-2' 'mTL-bodies'};
for rr=1:4
  cvnlookup('fsaverage',13,vflatten(mean(vals(:,2,:,:)==rr,4)),[0 1],copper(256),1/8/2,[],0);
  imwrite(rgbimg,sprintf('%s/fsaverage_probmap_%s.png',outputdir,roilabels{rr}));
end

%% then proceed to convert surface to volume (analysis_surfaceroistovolume.m)
