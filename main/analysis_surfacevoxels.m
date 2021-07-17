function analysis_surfacevoxels(subjix)

% function analysis_surfacevoxels(subjix)
%
% <subjix> is 1-8
%
% Save .mgz files with surface voxels at 1, 2, and 3 mm resolution
% onto each of the following surfaces: white, pial, layerB1-B3.

% define
hemis = {'lh' 'rh'};
surfs = {'white' 'pial' 'layerB1' 'layerB2' 'layerB3'};
mms = [1 2 3];
n = 256;

% do it
for subjix=1:8

  % make directory
  dir0 = sprintf('~/nsd/nsddata/freesurfer/subj%02d/label/surfacevoxels/',subjix);
  mkdirquiet(dir0);

  % create alternating volumes
  vols = [];
  for rr=1:length(mms)
    vol = zeros(n,n,n);
    ix = mod(ceil((1:n)/mms(rr)),2)==1;
    vol(ix,:,:) = vol(ix,:,:) + 1;
    vol(:,ix,:) = vol(:,ix,:) + 2;
    vol(:,:,ix) = vol(:,:,ix) + 4;
    vols(:,:,:,rr) = vol;
  end

  % process each hemisphere and surface
  for hh=1:length(hemis)
    for ss=1:length(surfs)
      nsd_mapdata(subjix,'anat1pt0',sprintf('%s.%s',hemis{hh},surfs{ss}),vols,'nearest',[], ...
        sprintf('%s/%s.surfacevoxels_%s.mgz',dir0,hemis{hh},surfs{ss}),[],'~/nsd/nsddata/freesurfer/subj%02d');
    end
  end
  
end
