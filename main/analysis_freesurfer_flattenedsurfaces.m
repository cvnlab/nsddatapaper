%% PREPARE GUIDELINES

% ensure SUBJECTS_DIR is correct in MATLAB environment

% load in what the fsaverage flat surface looks like
surfL = cvnreadsurface('fsaverage','lh','full.flat.patch.3d','orig');
surfR = cvnreadsurface('fsaverage','rh','full.flat.patch.3d','orig');

% map to individual subject's surface (LH)
subjix = 1;
nsd_mapdata(subjix,'fsaverage','lh.white',vflatten(double(surfL.patchmask)),[],[],'lh.temp.mgz',[], ...
  '/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj%02d');

% map to individual subject's surface (RH)
subjix = 1;
nsd_mapdata(subjix,'fsaverage','rh.white',vflatten(double(surfR.patchmask)),[],[],'rh.temp.mgz',[], ...
  '/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj%02d');

%% CUT THE SURFACES

% ensure SUBJECTS_DIR is correct in shell

% read instructions in FreeSurferOccipitalFlattenedPatch.pdf

tksurfer subj01 lh inflated
tksurfer subj02 lh inflated
tksurfer subj03 lh inflated
tksurfer subj04 lh inflated
tksurfer subj05 lh inflated
tksurfer subj06 lh inflated
tksurfer subj07 lh inflated
tksurfer subj08 lh inflated

tksurfer subj01 rh inflated
tksurfer subj02 rh inflated
tksurfer subj03 rh inflated
tksurfer subj04 rh inflated
tksurfer subj05 rh inflated
tksurfer subj06 rh inflated
tksurfer subj07 rh inflated
tksurfer subj08 rh inflated

steps:
- load curvature
- make binary display
- load overlay (the temp.mgz file) [select 'no registration']
- configure overlay, green-red, [0 1]
- five cuts (Cut Line) [generally straight line, pretty close to the folding of the lateral surface]
- medial wall (Cut Closed Line)
- rotate to a nice view and take snapshot
- click outside the big circle but on cortex (Fill Uncut Area)
- save patch as [lh,rh].full.patch.3d

%% FLATTEN THE PATCHES

for p=1:8
  subjid = sprintf('subj%02d',p);

  % prep
  cd(sprintf([cvnpath('freesurfer') '/%s/surf'],subjid));

  % perform the flattening
  unix('mris_flatten lh.full.patch.3d lh.full.flat.patch.3d');  %-w 10 
  unix('mris_flatten rh.full.patch.3d rh.full.flat.patch.3d');

end

%% ROTATE THE RESULTS

% go see analysis_freesurfer_flattenedsurfaces_rotation.m

%% CHECK THE RESULTS

% check the results
for p=1:8
  data = sprintf([nsd_datalocation '/freesurfer/subj%02d/label/lh.Kastner2015.mgz'],p);
  cvnlookup(sprintf('subj%02d',p),13,data,[0 25],jet(256),0.5,[],0);  % {'reset' true}
  imwrite(rgbimg,sprintf('~/Dropbox/KKTEMP/subj%02d.png',p));
  cvnlookup(sprintf('subj%02d',p),10,data,[0 25],jet(256),0.5,[],0);  % {'reset' true}
  imwrite(rgbimg,sprintf('~/Dropbox/KKTEMP/subj%02d_fsaverage.png',p));
end
