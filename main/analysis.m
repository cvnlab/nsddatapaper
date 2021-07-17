%% %%%%%%%% Manual data transfer, organization, and naming

% - check log entries on the google sheets
% - make official data directory with official name
% - rsync -av the files from naxos into "dicom"
%   - stick extra scans into "extra_scans"
%   - remember to check file/dir sizes
% - rsync -av .mat files and "data" and "scripts" folders into "mat_files_from_scan"
%   - check chronological order of the intended experiments and delete any unwanted files!
% - get slice.png (convert bmp to png)

% - for phantom test data:
%   - transfer data and organize and run autoqc_fmri only.
%   - hm, for DATECENSORED-ST001-NSD_phantom, just transferred only.

% hints below...

transcode eyetracking video [very fast 480] and store it away /research/nsd/eyetracking/

DATECENSORED-NSD139-structural4
0216
NSD139
subj139
139*

  # TYPE
mkdir /home/surly-raid1/kendrick-data/nsd/rawdata/DATECENSORED-NSD139-structural4/

  # TYPE
rsync -av ~/dicom/*0216*NSD139*session*/* /home/surly-raid1/kendrick-data/nsd/rawdata/DATECENSORED-NSD139-structural4/dicom/

get flicker .mat specially?

stimmac
cd ~/Desktop/cvnlab/kendrick/
  # TYPE
rsync -av 2019*subj139*     kendrick@stone.cmrr.umn.edu:"~/fmridata/nsd/rawdata/DATECENSORED-NSD139-structural4/mat_files_from_scan/"
rsync -av eye_2019*subj139* kendrick@stone.cmrr.umn.edu:"~/fmridata/nsd/rawdata/DATECENSORED-NSD139-structural4/eyetracking/"
mv *2019*subj139* \ old\ files/
  rsync -av fLoc/data/139*    kendrick@stone.cmrr.umn.edu:"~/fmridata/nsd/rawdata/DATECENSORED-NSD139-structural4/mat_files_from_scan/data/"
  rsync -av fLoc/scripts/139* kendrick@stone.cmrr.umn.edu:"~/fmridata/nsd/rawdata/DATECENSORED-NSD139-structural4/mat_files_from_scan/scripts/"
  rm -rf fLoc/data/139*
  rm -rf fLoc/scripts/139*

  # TYPE
convert /home/stone/tmp/*NSD139*T2*.bmp /home/surly-raid1/kendrick-data/nsd/rawdata/DATECENSORED-NSD139-structural4/slice_HRT2.png
rm -rf /home/stone/tmp/*NSD139*T2*.bmp
convert /home/stone/tmp/*NSD139*.bmp /home/surly-raid1/kendrick-data/nsd/rawdata/DATECENSORED-NSD139-structural4/slice.png
rm -rf /home/stone/tmp/*NSD139*.bmp

  # PHYSIO?
rsync -av /home/stone/tmp/Physio*2019* /home/surly-raid1/kendrick-data/nsd/rawdata/DATECENSORED-NSD139-structural4/physio/
rm -rf /home/stone/tmp/Physio*2019*
ls -lad ~/nsd/rawdata/DATECENSORED-NSD139-structural4/physio/* | wc -l
% 4x60 = 240; 300; 360; 420; 480; 540; 600; 660; 720; 780; 840; 900; 960; 1020; 1080; 1200; 1320; 1380; 1440; 1500; 1560; 1620; 1680
%   1740; 60√√
% do we need to remove excess files? (or to a new scan directory)

cd /home/surly-raid1/kendrick-data/nsd/rawdata/DATECENSORED-NSD139-structural4/mat_files_from_scan/
ls -la
rm *run99*
ls -latr * data/*/*.mat

cd /home/surly-raid1/kendrick-data/nsd/rawdata/DATECENSORED-NSD139-structural4/dicom/
mkdir extra_scans
mv *local* *SBRef extra_scans/
ls -la
mv MR-SE003* extra_scans/
du -hs *
du -s *

% consider editing nsdsetup.m for special resting-state extra runs??

% check physio files (exactly 15 .resp files) [ran on July 4 2019]
a1 = matchfiles('/home/surly-raid1/kendrick-data/nsd/rawdata/*/physio');
for p=1:length(a1)
  a1{p}
  b = unix_wrapper(sprintf('ls -1 %s/*.resp | wc -l',a1{p}));
  assert(str2double(b)==15);
end

% check that pre-processing is chugging away correctly [periodically]
./checknsdfinal

STATUS:
303 134syn√
304 258syn√
305 149syn√
306 929syn√
307 168√
308 814√
309 258imag√
310 134imag√
311 400syn√
312 139syn√
313 400imag√
314 139imag√
315 149imag√
316 814syn√
317 168syn√
318 168imag√
319 814imag√
320 929imag√

nsdsyntheticpilot: 134 exists [some fixation responses 3/1]

%% %%%%%%%% IF QUESTIONABLE SCAN SESSION, STOP HERE.

%% %%%%%%%% PERFORM ANY DATA FIXES IF NECESSARY!

% history:
% - datafix_script1.m
% - datafix_script2.m
% - datafix_script4.m

%% %%%%%%%% Manually define some parameters and upload/update

nsdsetup;

%% %%%%%%%% AUTOMATIC: start it

%% %%%%%%%% Quick check of AcquisitionDate and AcquisitionTime

% manually inspect the AcquisitionDate and AcquisitionTime saved in the DICOM
for zz=length(datadirs):length(datadirs)
  files0 = matchfiles(sprintf('%s/dicom/extra_scans/*',datadirs{zz}));
  files1 = matchfiles(sprintf('%s/*.dcm',files0{1}));
  a0 = dicominfo(files1{1});
  fprintf('%s: %s %s\n',datadirs{zz},a0.AcquisitionDate,a0.AcquisitionTime);
end

% assert to check that scan directory name has a date identical to the AcquisitionDate
for zz=1:length(datadirs), zz
  files0 = matchfiles(sprintf('%s/dicom/MR*bold*',datadirs{zz}));
  files1 = matchfiles(sprintf('%s/*.dcm',files0{1}));
  f = regexp(datadirs{zz},'.+?rawdata/(.+?)-(.+?)-.+$','tokens');
  a0 = dicominfo(files1{1});
  assert(isequal(f{1}{1},a0.AcquisitionDate));
  
  tmp = a0.PatientName.FamilyName(1:length(f{1}{2}));
  if ~isequal(f{1}{2},tmp)
    fprintf('WARNING: %d: %s\n',zz,tmp);
  end
end

%% %%%%%%%% Quick check of SliceLocation

  % what about
  % a.ImageOrientationPatient
  % a.ImagePositionPatient

for zz=length(datadirs):length(datadirs)
  autoqc_acquisition(zz);
end

%% %%%%%%%% Automatic QC of physio data

for zz=length(datadirs):length(datadirs)
  autoqc_physio(zz,'/home/stone/generic/Dropbox/nsd/autoqc_physio/');
end

%% %%%%%%%% Automatic QC and analysis of behavioral data

% prf
for zz=6+(1:14)
  analyzebehavior_prf(datadirs{zz});
end
analyzebehavior_prf_visualize(ppdirs(6+(1:14)),'/home/stone/generic/Dropbox/nsd/analyzebehavior_prf');

% floc
for zz=6+(1:14)
  analyzebehavior_floc(datadirs{zz},312,12);
end
analyzebehavior_floc_visualize(ppdirs(6+(1:14)),'/home/stone/generic/Dropbox/nsd/analyzebehavior_floc');

% nsd basic file check [this does not include nsdsynthetic or nsdimagery]
analyzebehavior_nsd_filecheck;

% nsd
for zz=21:length(datadirs)
  if ~ismember(zz,[nsdsyntheticsession nsdimagerysession])
    analyzebehavior_nsd(datadirs{zz},nsdexpfile);
  end
end
analyzebehavior_nsd_visualize(ppdirs(setdiff(21:length(ppdirs),[nsdsyntheticsession nsdimagerysession])),'/home/stone/generic/Dropbox/nsd/analyzebehavior_nsd');

% nsdsynthetic [ran it√]
for zz=21:length(datadirs)
  if ismember(zz,[nsdsyntheticsession])
    analyzebehavior_nsdsynthetic(datadirs{zz},nsdsyntheticexpfile);
  end
end

% nsdimagery [ran it√]
for zz=21:length(datadirs)
  if ismember(zz,[nsdimagerysession])
    analyzebehavior_nsdimagery(datadirs{zz},nsdimageryexpfile);
  end
end

% OBSOLETE (caused by synthetic)
% % nsd by subject
% mkdirquiet('/home/stone/generic/Dropbox/nsd/analyzebehavior_nsd_bysubject');
% for zz=1:length(allnsdids)
%   tmp0 = tempdir;
%   ix = cellfun(@(x) ~isempty(regexp(x,sprintf('%s-nsd',allnsdids{zz}))),dirs);
%   if any(ix)
%     analyzebehavior_nsd_visualize(ppdirs(ix),tmp0);
%     movefile([tmp0 '/accuracyrt.png'],sprintf('/home/stone/generic/Dropbox/nsd/analyzebehavior_nsd_bysubject/%s.png',allnsdids{zz}));
%   end
% end

%% %%%%%%%% Automatic QC of raw MRI data

for zz=length(datadirs):length(datadirs)
  autoqc_fmri(datadirs{zz},expecteddim{sessiontypes(zz)});
end
autoqc_fmri_visualize(ppdirs,'/home/stone/generic/Dropbox/nsd/autoqc_fmri',length(ppdirs):length(ppdirs));

% % THIS WAS EXPERIMENTAL. DEPRECATED.
% for zz=290:302
%   autoqc_fmri2(datadirs{zz});
% end

%% %%%%%%%% Perform volume-based pre-processing

% [hm, keep in mind that we don't yet inherit the correct NIFTI headers]
%   [also, field of view change, apparent voxel resolution change]
%
% [pipeline institute the distortion effects quantification?]
%   [epi and fieldmap might be interesting for large data set stability inspections?]
%
% note that the nsd sessions now include a fieldmap regularization based on magnitude (Nov 20 2018).
% 
% NOTES:
% - /(hard-coded change to fieldmap - regularization of the fieldmap implemented! 
%  this means my copy of the preprocess will do this always. it seems to run nicely.)
%  [circa Nov 20 2018. will need to check in eventually. is the solution robust?]
%  [methods notes: 1 after threshold; squaring before the threshold]

% Define motion-correction mask.
for zz=length(datadirs):length(datadirs)
  definemcmask(datadirs{zz});
end

% Perform volume-based pre-processing.
for zz=length(datadirs):length(datadirs)
  feval(ppscripts{pptypes(zz)},zz);
end

% Need to do any fixes to the pp data?
%
% history:
% - datafix_script3a.m

% Automatic QC of pre-processing results.
for zz=length(datadirs):length(datadirs)
  autoqc_ppfmri(ppdirs{zz});
end
autoqc_ppfmri_visualize(ppdirs,'/home/stone/generic/Dropbox/nsd/autoqc_ppfmri');

% Manual QC of pre-processing results.
% - EPIoriginal: any major artifacts?
% - EPIoriginal->EPIundistort->braincropped: basic sanity of the undistortion
% - EPIfinal: is overall stability okay?
% - meanvol, stdvol: are we missing brain in slices?
% - temporalsnr: is relatively stable and good?
% - fieldmapbraincropped, fieldmapsmoothedALT: any massive movement or problems?
% - randomvolumes: quick sanity check
% - motionallruns: anything crazy?
%
% Record observations on Google Sheets.
% Re-do any processing if necessary.

%% %%%%%%%% Analyze first-pass volume-based pp for the nsd screening sessions (prf and floc) and for sanity-checking nsdXX

% NOTE: we want to use the vanilla analyzePRF, so we do this:
rmpath(genpath('/home/stone/kendrick/code/analyzePRF/'));

% deal with prffloc
for zz=6+(1:14)
  glm_prfINITIALSCREENING(ppdirs{zz}, [1 2 5 6 9  10]);
  glm_flocINITIALSCREENING(ppdirs{zz},[3 4 7 8 11 12]);
end

% deal with nsdXX
for zz=21:length(ppdirs)
  if ~ismember(zz,[nsdsyntheticsession nsdimagerysession])
    glm_nsdSIMPLE(ppdirs{zz});  % one-regressor GLM
  end
end

% download and inspect figures

%% %%%%%%%% Select subjects based on nsd screening sessions

% figure: BOLDSNR
% figure: PARTICIPANTSUMMARY
% figure: FINALPARTICIPANTSUMMARY

%% %%%%%%%% Create leaderboard for main nsd experiment

autoqc_nsd(setdiff(21:length(ppdirs),[nsdsyntheticsession nsdimagerysession]),'/home/stone/generic/Dropbox/nsd/autoqc_nsd');

%% %%%%%%%% Perform gradunwarp on the mean volume from first-stage volume-based pre-processing.

for zz=7:length(datadirs), zz
  preprocess_nsd_epigradunwarp(datadirs{zz});
end

%% %%%%%%%% Based on gradunwarp results, determine crop range and insert into nsdsetup.m

% NOTE: only gets run for first nsd session!

for zz=length(datadirs):length(datadirs), zz
  if isempty(nsdcropranges{nsdsubjixfun(zz)})
    keyboard;
    %% manually go through preprocess_nsd_crop.m
    %% rehash path; clear nsdsetup;
    %% close nomachine!
  end
end

%% %%%%%%%% Align each session to the reference session [order matters and re-dos are tricky!].

  % do nsd01 first
for zz=[21 22 23 24 25 27 28 31]
  preprocess_nsd_epialignment(zz);
end
  % do screening sessions next
for zz=[7 10 18 16 12 17 8 20]
  preprocess_nsd_epialignment(zz);
end
  % then:
for zz=[26 29 30 32:length(datadirs)]
  preprocess_nsd_epialignment(zz);
end

%% %%%%%%%% Perform second-stage volume-based pre-processing.

for zz=[7 10 18 16 12 17 8 20 21:length(datadirs)]
  preprocess_nsd_secondstage(zz);
  if ismember(zz,[nsdimagerysession])
    preprocess_nsd_secondstageLOW1s(zz);   % IMAGERY LOW RESOLUTION is at 1s
  else
    preprocess_nsd_secondstageLOW(zz);
  end
end

% NOTES ON NAN/INVALID DATA:   [CHECK THIS]
% - After pre-processing, the bad voxels are marked
%   in valid.mat (basically, when valid.mat is 0, this means it is
%   a bad vertex). In the data files (e.g. run01.mat), the time-series
%   data for all bad vertices are set to all 0 (in all runs). Thus,
%   the data files only have finite values.

% Do we need to perform any fixes to the pp data?
%
% history:
% - datafix_script3b.m

%% %%%%%%%% Do some basic GLM analyses.

% Perform GLM.
for zz=21:length(datadirs)
  if ~ismember(zz,[nsdsyntheticsession nsdimagerysession])

    glm_nsd(ppdirs{zz},'preprocessVER1SECONDSTAGE',   'designmatrixSINGLETRIAL_nsd.mat',           ...
      '~/nsd/ppdata/hrfbasis.mat',          1,  'glmdata'   ,'~/nsd/ppdata/hrfmanifold.mat');
    glm_nsd(ppdirs{zz},'preprocessVER1SECONDSTAGELOW','designmatrixSINGLETRIAL_3pertrial_nsd.mat', ...
      '~/nsd/ppdata/hrfbasis_3pertrial.mat',4/3,'glmdataLOW','~/nsd/ppdata/hrfmanifold_3pertrial.mat');

    % Resting-state
    if sessiontypes(zz)==8
      glm_nsdrestingstate(ppdirs{zz},'preprocessVER1SECONDSTAGE',   'designmatrixSINGLETRIAL_nsd.mat',           ...
        '~/nsd/ppdata/hrfbasis.mat',          1,  'glmdata'   ,'~/nsd/ppdata/hrfmanifold.mat');
      glm_nsdrestingstate(ppdirs{zz},'preprocessVER1SECONDSTAGELOW','designmatrixSINGLETRIAL_3pertrial_nsd.mat', ...
        '~/nsd/ppdata/hrfbasis_3pertrial.mat',4/3,'glmdataLOW','~/nsd/ppdata/hrfmanifold_3pertrial.mat');
    end
  
  end

  % deal with nsdsynthetic
  if ismember(zz,[nsdsyntheticsession])
    glm_nsdsynthetic(ppdirs{zz},'preprocessVER1SECONDSTAGE',   'designmatrixSINGLETRIAL_nsdsynthetic.mat',           ...
      '~/nsd/ppdata/hrfbasis.mat',          1,  'glmdata'   ,'~/nsd/ppdata/hrfmanifold.mat');
    glm_nsdsynthetic(ppdirs{zz},'preprocessVER1SECONDSTAGELOW','designmatrixSINGLETRIAL_3pertrial_nsdsynthetic.mat', ...
      '~/nsd/ppdata/hrfbasis_3pertrial.mat',4/3,'glmdataLOW','~/nsd/ppdata/hrfmanifold_3pertrial.mat');
  end

  % deal with nsdimagery
  if ismember(zz,[nsdimagerysession])
    glm_nsdimagery(ppdirs{zz},'preprocessVER1SECONDSTAGE',   'designmatrixSINGLETRIAL_nsdimagery.mat',           ...
      '~/nsd/ppdata/hrfbasis.mat',          1,'glmdata'   ,'~/nsd/ppdata/hrfmanifold.mat');
    glm_nsdimagery(ppdirs{zz},'preprocessVER1SECONDSTAGELOW','designmatrixSINGLETRIAL_nsdimagery.mat', ...
      '~/nsd/ppdata/hrfbasis.mat',          1,'glmdataLOW','~/nsd/ppdata/hrfmanifold.mat');
  end

end

% Automatic QC of GLM results.
for zz=[7 10 18 16 12 17 8 20 21:length(datadirs)]
  autoqc_glm_nsd(zz,'glmdata');
  autoqc_glm_nsd(zz,'glmdataLOW');
end

%% %%%%%%%% Do some grand finales.

for ss=1:8
  autoqc_nsd_grand(ss,'glmdata');
  autoqc_nsd_grand(ss,'glmdataLOW');
end
autoqc_nsd_grandfinale('glmdata',42);
autoqc_nsd_grandfinale('glmdataLOW',42);

%% %%%%%%%% AUTOMATIC END:
%%          1. look at PP QC figures in sessions/nsd (see above questions)
%%             exclude MOVIEfinal, MOVIEoriginal, randomvolumes, contours* [should we delete from server?]
%%          2. look at GLMdenoise_nsdSIMPLE_figures, sessioncoregistration, glmdata-glmBASIC 
%%          3. look at nsd Dropbox figures
%%          4. e-mail subject with the results.
%%          5. write bonus on spreadsheet
