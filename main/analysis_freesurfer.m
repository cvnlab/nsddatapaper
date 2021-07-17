%% %%%%%%%% Manually define parameters

nsdsetup;

%% %%%%%%%% Run FreeSurfer

% install the runtime into the FS60 folder:
%  curl "https://surfer.nmr.mgh.harvard.edu/fswiki/MatlabRuntime?action=AttachFile&do=get&target=runtime2012bLinux.tar.gz" -o "runtime.tar.gz"
%  tar xvf runtime.tar.gz

% prepare an expert options file
echo "mris_inflate -n 50" > ~/nsddata/templates/expert.opts

% run on stone!
for zz=1:length(allnsdids)
  % build-stamp: freesurfer-x86_64-unknown-linux-gnu-stable6-20170118:
  cvnrunfreesurfer(sprintf('FS6EXPERT_subj%02d',zz), ...
    sprintf('%s/ppdata/%s-structurals/T1_0pt8_masked.nii.gz',nsddir,allnsdids{zz}), ...
    sprintf('-hires -openmp 32 -norandomness -expert /home/surly-raid4/kendrick-data/nsd/nsddata/templates/expert.opts -brainstem-structures -hippocampal-subfields-T1T2 %s/ppdata/%s-structurals/T2_0pt8_masked.nii.gz HST1T2',nsddir,allnsdids{zz}));
end

% at this point, the FSIDs are like "FS6EXPERT_subj01".

%% %%%%%%%% Run FreeSurfer for the individual reps [SIDE EFFORT]

% run on stone! [it seems to be actually: freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0-2beb96c]
%   
% edit .tcshrc:
% setenv SUBJECTS_DIR /home/surly-raid4/kendrick-data/nsd/nsddata_other/freesurferoriginals
%
% matlab:
% setenv('SUBJECTS_DIR','/home/surly-raid4/kendrick-data/nsd/nsddata_other/freesurferoriginals');

% process all individual T1s from all subjects
for zz=1:8
  for rep=1:10

    % check for the file
    file0 = sprintf('%s/ppdata/%s-structurals/T1_rep%02d_0pt8_masked.nii.gz',nsddir,allnsdids{zz},rep);
    if ~exist(file0,'file')
      continue;
    end

    % build-stamp: freesurfer-x86_64-unknown-linux-gnu-stable6-20170118:
    cvnrunfreesurfer(sprintf('subj%02d_rep%02d',zz,rep), ...
      file0, ...
      sprintf('-hires -openmp 32 -norandomness -expert /home/surly-raid4/kendrick-data/nsd/nsddata/templates/expert.opts -brainstem-structures -hippocampal-subfields-T1'));  % -hippocampal-subfields-T1T2

  end
end  

% Notice that the processing is exactly the same except for the T1T2 hippocampus stuff.
%
% For the individual-reps analysis, we stop the processing here. No posthoc commands.

%% %%%%%%%% Visual inspection of surface reconstructions (white and pial on T1)

% create image directory
mkdir /home/stone/kendrick/Dropbox/KKTEMP/FS6EXPERT_subj08/
cd /home/stone-ext1/freesurfer/subjects/FS6EXPERT_subj08
cd ~/nsd/nsddata/freesurfer/subj08

% FOR INITIAL PASS VISUALIZATION
freeview -v mri/brain.finalsurfs.mgz \
  -f surf/lh.white:edgecolor=blue:edgethickness=1 \
     surf/lh.pial:edgecolor=cyan:edgethickness=1 \
     surf/rh.white:edgecolor=red:edgethickness=1 \
     surf/rh.pial:edgecolor=yellow:edgethickness=1

% FOR AFTER-FS-EDITS VISUALIZATION
freeview -v mri/T1.mgz \
  -v mri/surfaceimperfections.mgz \
  -f surf/lh.white:edgecolor=blue:edgethickness=1 \
     surf/lh.pial:edgecolor=cyan:edgethickness=1 \
     surf/rh.white:edgecolor=red:edgethickness=1 \
     surf/rh.pial:edgecolor=yellow:edgethickness=1

% For surfaceimperfections, use NIH 1.5 3.5 and opacity 0.5.

% save movie frames
make the visible area a large square by adjusting size of panels
for three orthogonal views, save movie frames using slices 0 through 340,
  saving the images to something like: subj01/axial

% consolidate into a single .mp4 file
views = {'sagittal' 'coronal' 'axial'};
for p=1:8
  for q=1:length(views)
    dir0 = sprintf('/Users/kendrick/Desktop/subj%02d/%s',p,views{q});
%     imagesequencetomovie(dir0, ...
%       sprintf('%s.mov',dir0),15);
    unix_wrapper(sprintf('ffmpeg -framerate 15 -pattern_type glob -i ''%s/*.png'' -crf 18 -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p %s.mp4',dir0,dir0));
  end
end

% post-process .mp4 with Handbrake:
% Fast 1080p:
%  sagittal crop 700x4
%  coronal crop 1000 for left/right   700 for top/down
%  axial crop 1000 for left/right   400 for top/down

%% %%%%%%%% Manual editing of FreeSurfer surfaces

% Follow procedure described by Luca and Ruyuan (see wiki).
% This procedure generates a new directory named something like:
%   /stone/ext1/freesurfer/subjects/C1051
% and has a special file named mri/binaryvolume.mgz

% Generate visualization of the FreeSurfer white and pial surfaces
% for the new FS subject ID (see above).

% Manually inspect to ensure sanity and quality!!!

% Get new FS directories in the right location
rsync -av /stone/ext1/freesurfer/subjects/FS6EXPERT_subj01_ver10/* /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj01/
rsync -av /stone/ext1/freesurfer/subjects/FS6EXPERT_subj02_ver7/*  /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj02/
rsync -av /stone/ext1/freesurfer/subjects/FS6EXPERT_subj03_ver5/*  /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj03/
rsync -av /stone/ext1/freesurfer/subjects/FS6EXPERT_subj04_ver5/*  /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj04/
rsync -av /stone/ext1/freesurfer/subjects/FS6EXPERT_subj05_ver6/*  /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj05/
rsync -av /stone/ext1/freesurfer/subjects/FS6EXPERT_subj06_ver5/*  /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj06/
rsync -av /stone/ext1/freesurfer/subjects/FS6EXPERT_subj07_ver5/*  /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj07/
rsync -av /stone/ext1/freesurfer/subjects/FS6EXPERT_subj08_ver5/*  /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj08/

% MORE! [do we have to do anything more?]
rsync -av /home/stone-ext1/freesurfer/subjects/fsaverage3 /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/
rsync -av /home/stone-ext1/freesurfer/subjects/fsaverage4 /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/
rsync -av /home/stone-ext1/freesurfer/subjects/fsaverage5 /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/
rsync -av /home/stone-ext1/freesurfer/subjects/fsaverage6 /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/
rsync -av /home/stone-ext1/freesurfer/subjects/fsaverage_sym /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/

% Make group-writable
sudo chmod -R g+w /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj??

% Remove stale imglookup?
ls -lad /home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj??/surf/imglookup

% do some binary volume preparation
for p=1:8
  cvnpreparebinaryvolume(sprintf('subj%02d',p));
end

% use freeview to edit the volume
cd ~/nsd/nsddata/freesurfer/subj01/
cd /Volumes/SeagateEXfa/NSD_SEGMENTATION_EA_2019/binaryvolume/subj01/
freeview \
  -v mri/T2.mgz \
  -v mri/T1.mgz \
  -v mri/binaryvolume.mgz \
  -f surf/lh.white:edgecolor=blue:edgethickness=1 \
     surf/lh.pial:edgecolor=cyan:edgethickness=1 \
     surf/rh.white:edgecolor=red:edgethickness=1 \
     surf/rh.pial:edgecolor=yellow:edgethickness=1

% The binaryvolume is initially 1 in brain tissue.
% The goal is to label some of these binaryvolume voxels as "2".
% Don't restrict the labeling of 2 and it's okay to be somewhat sloppy
%   (can extend beyond gray matter a bit).  We just want to mark gray matter
%   that might be suspect in terms of bad FreeSurfer surfaces.
% Save the drawn volume back to binaryvolume.mgz.

% Rename:
% mv mri/binaryvolume.mgz mri/surfaceimperfections.mgz

%% %%%%%%%% Perform posthoc FreeSurfer commands

  % build-stamp: freesurfer-x86_64-unknown-linux-gnu-stable6-20170118:
cvnrunfreesurfer2('fsaverage','-hires -openmp 32 -norandomness');  % THIS IS A BIT EXPERIMENTAL. IS THE MEANING OF THIS VALID??   [DID NOT RE-RUN FOR FS6]
for p=1:8
  cvnrunfreesurfer2(sprintf('subj%02d',p),'-hires -openmp 32 -norandomness');
end

%% %%%%%%%% Generate dense layer truncated surfaces

% The dense and truncated are probably not of interest, but let's just bring it along!

parpool;

  % build-stamp: freesurfer-x86_64-unknown-linux-gnu-stable6-20170118:
cvnmakelayers for 'fsaverage' --> go run the hack in "02 general fsaverage pre-processing.m"   [WELL, DID NOT RE-RUN FOR FS6. JUST INHERIT.]
for p=1:8
  cvnmakelayers(sprintf('subj%02d',p),linspace(.25,.75,3),'B','pt');
end

%% %%%%%%%% Get fsaverage atlases onto individual subjects' surfaces

analysis_freesurfer_atlas;

%% %%%%%%%% Calculate lookups for transferring to/from fsaverage

analysis_freesurfer_fsaveragetransfer;

%% %%%%%%%% Prepare flattened patches

analysis_freesurfer_flatgrids;

%% %%%%%%%% Make fsaverage flat surface

% Note: we already ran "05e fsaverage flat.m" one-time-only for fsaverage.

%% %%%%%%%% Calculate transformations between spaces

% Note: June 5 2019: only needed to run modes 3-9.
% Note: June 6 2019: fixed bug (invwarp) and re-ran modes 6-9.
% Note: June 23 2019, re-run all for edited surfaces.
% Note: April 19, 2020, re-run mode [3 8] (fix bug in FS interpretation)

for subjix=1:8
  for mode=1:9
    preprocess_nsd_calculatetransformations(subjix,mode);
  end
end

% make sure the function nsd_mapdata.m is written and correct!

%% %%%%%%%% Deal with flattened surfaces

% go do analysis_freesurfer_flattenedsurfaces.m

%% %%%%%%%% Semi-inflated brain

% generate semi-inflated brain
hemis = {'lh' 'rh'};
for p=1:8
  for hh=1:length(hemis)
    hemi = hemis{hh};
    dir0 = sprintf('~/nsddata/freesurfer/subj%02d/surf/',p);
    assert(movefile(sprintf('%s/%s.sulc',dir0,hemi),sprintf('%s/%s.sulc.bak',dir0,hemi)));  % KEEP A BACKUP
    unix_wrapper(sprintf('mris_inflate -n 5 %s/%s.smoothwm %s/%s.semiinflated',dir0,hemi,dir0,hemi));
    assert(movefile(sprintf('%s/%s.sulc',dir0,hemi),sprintf('%s/%s.sulcsemiinflated',dir0,hemi)));
    assert(movefile(sprintf('%s/%s.sulc.bak',dir0,hemi),sprintf('%s/%s.sulc',dir0,hemi)));
    %mris_inflate -w 1 -n 50 ~/nsddata/freesurfer/subj08/surf/lh.smoothwm ~/nsddata/freesurfer/subj08/surf/lh.semiinflated
  end
end
