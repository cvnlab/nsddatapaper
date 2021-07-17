%% This document shows how fsaverage-defined ROIs are converted to individual-subject
%% surface ROIs, and then how those individual-subject surface ROIs are converted
%% to the anat0pt8, func1mm, and func1pt8mm spaces.


% First, we will take various ROIs defined on fsaverage and transfer these ROIs
% via nearest-neighbor to the individual subject surfaces.
hemis = {'lh' 'rh'};
rois = {'Kastner2015' 'HCP_MMP1' 'gVTC' 'hVTC' 'gEVC' 'visualsulc' 'nsdgeneral' 'corticalsulc' 'streams' 'fVTC'};
for subjix=1:8
  for r=1:length(rois)
    for p=1:length(hemis)
      fsdir = sprintf([nsd_datalocation '/freesurfer/subj%02d'],subjix);
      sourcedata = sprintf([nsd_datalocation '/freesurfer/fsaverage/label/%s.%s.mgz'],hemis{p},rois{r});
      outputfile = sprintf([nsd_datalocation '/freesurfer/subj%02d/label/%s.%s.mgz'],subjix,hemis{p},rois{r});
      nsd_mapdata(subjix,'fsaverage',sprintf('%s.white',hemis{p}),sourcedata,[],[],outputfile,[],fsdir);
        % copy the .ctab and .txt file too:
      assert(copyfile(sprintf([nsd_datalocation '/freesurfer/fsaverage/label/%s.mgz.*'],rois{r}), ...
                      sprintf([nsd_datalocation '/freesurfer/subj%02d/label/'],subjix)));
    end
  end
end


% Next, we will take each ROI file and convert its labels to a volume in the
% anat0pt8 space, and then convert that volume to the func1mm and func1pt8mm spaces.
% Note that we have to be careful about handling the hemispheres: if the
% hemispheres are processed separately, there could be a voxel that "belongs"
% to both hemispheres, which would be weird.
hemis = {'lh' 'rh'};
rois = {'Kastner2015' 'HCP_MMP1' 'visualsulc' 'nsdgeneral' 'corticalsulc' 'streams' 'fVTC' 'floc-faces' 'floc-words' 'prf-visualrois' 'prf-eccrois' 'floc-places' 'floc-bodies'};
Ns = [25 180 14 1 28 7 1 5 5 7 5 3 4];
for subjix=1:8
  for r=1:length(rois)
    N = Ns(r);
  
    % load labels from both hemispheres
    data = {};
    for p=1:length(hemis)
      
      % load labeling
      data0 = cvnloadmgz(sprintf('%s/freesurfer/subj%02d/label/%s.%s.mgz',nsd_datalocation,subjix,hemis{p},rois{r}));

      % sanity check that all labels are integers between 0 and N
      assert(all(isfinite(data0)));
      temp = flatten(union(data0(:),[]));
      assert(all(ismember(temp,0:N)));
      
      % record
      data{p} = data0;  % column vector

    end
    
    % apply an offset such that the rh values live from N+1 onwards
    data{2} = data{2} + (N+1);
    
    % perform the mapping. after this step, non-cortical voxels are set to -1,
    % and cortical voxels range between 0 and N (LH) and N+1 and N+1+N (RH).
    mapdata = nsd_mapdata(subjix,{'lh.layerB1' 'lh.layerB2' 'lh.layerB3' ...
                                  'rh.layerB1' 'rh.layerB2' 'rh.layerB3'}, ...
                          'anat0pt8',data([1 1 1 2 2 2]),'surfacewta',-1);

    % create some other versions: use a winner-take-all scheme to map the 
    % labels from the anat0pt8 space to the func1mm and func1pt8mm spaces.
    mapdata2 = nsd_mapdata(subjix,'anat0pt8','func1pt0',mapdata,'wta',-1);
    mapdata3 = nsd_mapdata(subjix,'anat0pt8','func1pt8',mapdata,'wta',-1);

    % now we have to separate out the left and right hemispheres, taking
    % care to remove the offset from the rh volume values. after this step,
    % both the lh and rh contain values between -1 and N. a value of -1
    % means "not in this hemisphere's cortical voxels", and a value of 0
    % means "in this hemisphere's cortical voxels but not explicitly labeled".
    alldata = {mapdata mapdata2 mapdata3};
    dirnames = {'anat' 'func1mm' 'func1pt8mm'};
    targetres = [0.8 1 1.8];
    for p=1:length(alldata)
        % left
      hemi = 'lh';
      outputfile = sprintf([nsd_datalocation '/ppdata/subj%02d/%s/roi/%s.%s.nii.gz'],subjix,dirnames{p},hemi,rois{r});
      vollh = copymatrix(alldata{p},alldata{p} >= (N+1),-1);                 % set RH values to -1
      nsd_savenifti(vollh,repmat(targetres(p),[1 3]),outputfile);
        % right
      hemi = 'rh';
      outputfile = sprintf([nsd_datalocation '/ppdata/subj%02d/%s/roi/%s.%s.nii.gz'],subjix,dirnames{p},hemi,rois{r});
      volrh = copymatrix(alldata{p},alldata{p} >= 0 & alldata{p} <= N,-1);   % set LH values to -1
      volrh(alldata{p} >= (N+1)) = volrh(alldata{p} >= (N+1)) - (N+1);       % remove offset from RH values
      nsd_savenifti(volrh,repmat(targetres(p),[1 3]),outputfile);
        % combined
      outputfile = sprintf([nsd_datalocation '/ppdata/subj%02d/%s/roi/%s.nii.gz'],subjix,dirnames{p},rois{r});
      volboth = vollh;
      volboth(vollh==-1) = volrh(vollh==-1);
      nsd_savenifti(volboth,repmat(targetres(p),[1 3]),outputfile);
    end

  end
end
