%% %%%%%%%% Organize the data

% manually organize dicom directories

% quick check of directory sizes
du -s ~/nsd/rawdata/*stru*/dicom/*

%% %%%%%%%% Prepare some templates

% get MNI152_T1_1mm templates
cd /opt/local/fsl-5.0.7/fsl/data/standard
rsync -av MNI152_T1_1mm* /home/stone/generic/Dropbox/nsddata/templates/

% create a version of the brain mask with extra dilation
fslmaths MNI152_T1_1mm_brain_mask_dil.nii.gz -kernel box 9 -dilM MNI152_T1_1mm_brain_mask_dil_dilM.nii.gz

% manually create MNI152_T1_1mm_LPIdot.gz using ITK-SNAP.
% this is a simple binary mark with a big dot located at left, posterior, inferior.

%% %%%%%%%% Manually define parameters

nsdsetup;

%% %%%%%%%% Convert to NIFTI and possibly perform gradunwarp

% strategy:
% - T1w,T2w,dMRI,SpinEchoFieldMap,TOF -> NIFTI + gradunwarp
% - HRT2 (at 7T) -> NIFTI
% - SWI -> NIFTI

% make NIFTIs for each subject
parpool;
for zz=1:length(allnsdids)    % [1 2 3 4 5 6 7 8] for HRT2
  outputdir = sprintf('%s/ppdata/%s-structurals',nsddir,allnsdids{zz});
  cvncollectT1s({outputdir},structuraldirs{zz},'prisma','T1w',1);
  cvncollectT1s({outputdir},structuraldirs{zz},'prisma','T2w',1);
  cvncollectT1s({outputdir},structuraldirs{zz},'prisma','dMRI',0);
  cvncollectT1s({outputdir},structuraldirs{zz},'prisma','SpinEchoFieldMap',0);
  cvncollectT1s({outputdir},structuraldirs{zz},'prisma','TOF',0);

  dir0 = stripfile(stripfile(subscript(matchfiles(sprintf('~/nsd/rawdata/*%s*/dicom/*t2_tse*',allnsdids{zz})),1,1)))
  cvncollectT1s({outputdir},dir0,[],'t2_tse',0);  % no gradunwarp.

  % SWI
  dir0 = stripfile(stripfile(subscript(matchfiles(sprintf('~/nsd/rawdata/*%s*/dicom/*Mag_Images',allnsdids{zz})),1,1)))
  cvncollectT1s({outputdir},dir0,'7tps','Mag_Images',0);  % gradunwarp just in case we need it
  cvncollectT1s({outputdir},dir0,'7tps','Pha_Images',0);  % gradunwarp just in case we need it
end
% Note: Data might reside in multiple directories!
% Note: Skip BIAS because we don't care.

% For HRT2, go ahead and export into nsd files
for p=1:8
  outputdir = sprintf('%s/ppdata/%s-structurals',nsddir,allnsdids{p});
  file0 = matchfiles(sprintf('%s/*t2tsecor*.nii.gz',outputdir));
  dir0 = sprintf('~/nsddata/ppdata/subj%02d/anat/HRT2/',p);
  rmdirquiet(dir0);
  mkdirquiet(dir0);
  copyfile(file0{1},sprintf('%s/HRT2_raw.nii.gz',dir0));
end

%% %%%%%%%% Prepare T1s and T2s

% for each subject, prepare the T1s and T2s (coregistration, gradunwarp, multiple resolutions)
% note that this step also brings along the individual anatomical "rep" volumes.
for zz=1:8

  % prepare T1s
  preprocess_nsd_structuralalignment(zz,1);
  
  % prepare T2s
  if ismember(zz,[2])
    preprocess_nsd_structuralalignment(zz,2,1);  % pause to get better initial seed for T2 to T1 for split-sessions
  else
    preprocess_nsd_structuralalignment(zz,2);
  end

end

% check the results by looking at the orthMID figures.
% reject any volumes and re-run if necessary!

%% %%%%%%%% Prepare SWI and TOF

for zz=1:8
  preprocess_nsd_SWI(zz);
  preprocess_nsd_TOF(zz);
end

%% %%%%%%%% Determine MNI alignment and obtain brain mask

% Based on the prepared 1.0-mm T1, determine nonlinear warp to MNI.
% We obtain a subject-native brain mask using the determined warp.
for zz=1:8
  preprocess_nsd_MNIandbrainmask(zz);
end

% Note: June 6 2019: re-generated the inverse warp (standard2str). this fixed a bug. but we did not re-generate the brain mask!!

% Check the sanity of the T1_to_MNI_nonlin.nii.gz [across all subjects].
% Check the sanity of the brainmask_1pt0.nii.gz.
% Check the sanity of the inverse warp.

%% %%%%%%%% Apply brain mask in order to de-identify the structural volumes

% apply brainmask to the T1s, T2s, SWIs, TOFs (at all resolutions).
% note that this step also prepares the individual anatomical "rep" volumes.
for zz=1:8
  preprocess_nsd_applybrainmask(zz);
end

% the resulting files are now the official anatomical volumes.
% check that these files do not omit any brain!

%% %%%%%%%% Determine alignment between EPI and T2

% use ANTS to match the EPI (mean of first 5 sessions at 1.0-mm resolution)
% to the T2 (1.0-mm version). the nonlinear warp is really what we want.
for zz=1:8
  preprocess_nsd_functionaltostructuralalignment(zz);
end

%% %%%%%%%% Check results!!

% write out some orthMID inspections of all modalities
for zz=1:8
  preprocess_nsd_structuralinspection(zz);
end

% make a version showing all subjects
autoqc_nsd_grandfinale('structurals');

% visually inspect to ensure sanity.
%
% thus far, note that we quality-controlled things like:
% - no imaging artifacts, good alignment, good liberal brain mask.

%% %%%%%%%% Add surface voxels

for zz=1:8
  analysis_surfacevoxels(zz);
end
