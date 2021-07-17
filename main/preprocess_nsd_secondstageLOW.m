function preprocess_nsd_secondstageLOW(sessix)


nsdsetup;
datadir = datadirs{sessix};


% start parallel MATLAB to speed up execution.
numworkers = 6;  % keep low to prevent memory swapping!?
if isempty(gcp('nocreate'))
  parpool(numworkers);
end

% *** EDIT ME ***
% what directory do the data live in?
%%%%%datadir = '/home/stone-ext1/fmridata/20180601-ST001-CVNS001ODC/';
ppdir = regexprep(datadir,'rawdata','ppdata');
mkdirquiet(ppdir);


% *** EDIT ME ***
% where should i save figures to?
figuredir = [ppdir '/preprocessVER1SECONDSTAGELOWfigures'];

  %%% CASE 2: separate DICOM directories [if you use this case, please comment out CASE 1 and CASE 3]

  % the first set of DICOMs are the fieldmaps; the second set of DICOMs are the magnitude brains.
  % after loading in the second set of DICOMs, we consider only the first N slices, where N is
  % the number of slices in the first set of DICOMs (the reason is that there may be two volumes
  % in the second set of DICOMs).
%   fieldmapfiles = cellfun(@(x,y) {x y}, ...
%     matchfiles([datadir '/dicom/PH*field_mapping*']), ...
%     matchfiles([datadir '/dicom/MR*field_mapping*']), ...
%     'UniformOutput',0);
  fieldmapfiles = matchfiles([datadir '/dicom/*field_mapping*']);
  fieldmapfiles = cellfun(@(x,y) {x y}, ...
    fieldmapfiles(2:2:end), ...
    fieldmapfiles(1:2:end), ...
    'UniformOutput',0);

% %%%% OOPS.  THE FIELDMAPS HAD TOO MANY SLICES IN THIS SESSION (ESSA MISTAKE)
% %%%% SO WE USE THIS HACK:
% fieldmapslicerange = 8:42-7;

% A NEW TWEAK/ADJUSTMENT:
fieldmapslicerangeALT = fieldmapslicerangeALTs{sessix};




% what DICOM directories should we interpret as EPI runs?
% it is okay if you also match things that are not DICOM directories; we'll just ignore them.
% *** CHECK ME ***
epifilenames = matchfiles([datadir '/dicom/*bold*']);
isok = @(x) isempty(regexp(x,'SBRef'));
epifilenames = epifilenames(cellfun(isok,epifilenames));
%%epifilenames = epifilenames([1 4 5 6]);  % had to ignore some crap runs



if isempty(fixepifuns{sessix})
  fixepifun = [];
else
  fixepifun = fixepifuns{sessix};
end


  
% if you didn't acquire fieldmaps with the same slice thickness as the 
% functionals, we can work around this problem if your fieldmaps are
% a positive integer multiple of the slice thickness of the functionals,
% and if the total field-of-view in the slice dimension is the same.
% all we do is upsample the fieldmaps using nearest neighbor interpolation.
% this is done immediately and we then act as if the fieldmaps were acquired
% at the correct slice thickness.  (of course, we could be more flexible
% and fix other circumstances, but we'll do this as the need arises.)
% if you want the work around, supply the appropriate positive integer 
% for <fieldmapslicefactor>.  if [], do nothing special.
% *** EDIT ME ***
fieldmapslicefactor = 2;

% what are the time values to associate with the fieldmaps?
% if [], default to 1:N where N is the number of fieldmaps.
% *** EDIT ME ***
if isempty(fieldmaptweaks{sessix})
  fieldmaptimes = [1 5 9 13];
else
  fieldmaptimes = fieldmaptweaks{sessix}{1};
end

% what is the difference in TE (in milliseconds) for the two volumes in the fieldmaps?
% *** CHECK ME ***
fieldmapdeltate = 1.02;

% should we attempt to unwrap the fieldmaps? (note that 1 defaults to a fast, 2D-based strategy; 
% see preprocessfmri.m for details.)  if accuracy is really important to you and the 2D strategy 
% does not produce good results, consider switching to a full 3D strategy like 
% fieldmapunwrap = '-f -t 0' (however, execution time may be very long).
fieldmapunwrap = 1;

% how much smoothing (in millimeters) along each dimension should we use for the fieldmaps?
% the optimal amount will depend on what part of the brain you care about.
% I have found that 7.5 mm may be a good general setting.
% The extra 1 indicates extra regularization.
fieldmapsmoothing = [5 5 5 1];  %[7.5 7.5 7.5];






% what DICOM directories should we interpret as in-plane runs?
% it is okay if you also match things that are not DICOM directories; we'll just ignore them.
inplanefilenames = [];   %matchfiles([datadir '/dicom/*Inplane*']);






% do you want the special "push alternative data" mode?  if so, specify the 'record.mat'
% file from the previous call and <mcmask> must not be {}.  otherwise, leave as [].
%wantpushalt = [];
%wantpushalt = [ppdir '/preprocessVER1figures/record.mat'];

wantpushalt = {[ppdir '/preprocessVER1figures/record.mat'] 1 0 1};  % since we change temporal rate, we need to re-do motion



% this input is important only if you acquired something other than just the magnitude images
% for the EPI data.  the format is documented in dicomloaddir.m (the input <phasemode>).
% for the purposes of this script, the second element of <epiphasemode> must have exactly one
% element (and that element should probably be either 1 or 5).  if you acquired just the 
% magnitude images, just leave epiphasemode set to [].
epiphasemode = [];

% what is the desired in-plane matrix size for the EPI data?
% this is useful for downsampling your data (in order to save memory) 
% in the case that the data were reconstructed at too high a resolution.  
% for example, if your original in-plane matrix size was 70 x 70, the 
% images might be reconstructed at 128 x 128, in which case you could 
% pass in [70 70].  what we do is to immediately downsample each slice
% using lanczos3 interpolation.  if [] or not supplied, we do nothing special.
epidesiredinplanesize = [];

% what is the slice order for the EPI runs?
% special case is [] which means to omit slice time correction.

  %episliceorder = 'interleavedalt';

  % HOW TO FIGURE OUT THE SLICE PARAMETER
  % a = dicominfo('MR-ST001-SE011-0041.dcm');
  % mbfactor = 3;
  % [d,ix] = sort(a.Private_0019_1029(1:end/mbfactor));
  % ord = [];
  % for p=1:length(ix)
  %   ord(ix(p)) = p;
  % end
  % ord = repmat(ord,[1 mbfactor]);
  % mat2str(ord)
  % 
  % % KEITH'S WAY:
  % dcmfile = 'MR-ST001-SE011-0041.dcm';
  % D=dicominfo(dcmfile);
  % slicetimes=D.Private_0019_1029;
  % slicediff=diff(sort(slicetimes));
  % slicediff=median(slicediff(slicediff>0));
  % sliceindex=(round(slicetimes/slicediff)+1)';
  % mat2str(sliceindex)


% *** EDIT ME ***

% 84 slices, MB3

% get a new sampling rate (cubic interpolation)
episliceorder = {};
episliceorder{1} = [1 14 27 12 25 10 23 8 21 6 19 4 17 2 15 28 13 26 11 24 9 22 7 20 5 18 3 16 1 14 27 12 25 10 23 8 21 6 19 4 17 2 15 28 13 26 11 24 9 22 7 20 5 18 3 16 1 14 27 12 25 10 23 8 21 6 19 4 17 2 15 28 13 26 11 24 9 22 7 20 5 18 3 16];
episliceorder{2} = epitrtweaks{sessix} * (4/3);




% what is the phase-encode direction for the EPI runs? (see preprocessfmri.m for details.)
% (note that 'InPlanePhaseEncodingDirection' in the dicominfo of the EPI files gives either COL 
% which means the phase-encode direction is up-down in the images (which is 1 or -1 in our
% convention) or ROW which means the direction is left-right in the images (which is 2 or -2 in
% our convention).  the problem is that we don't know if the phase direction has been flipped,
% which you can do manually via CV vars.  it appears that if you don't explicitly flip, you should
% use the positive version (i.e. 1 or 2), meaning that the direction is towards the bottom
% (in the case of 1) or to the right (in the case of 2).  but in any case you should always
% check the sanity of the results!
%epiphasedir = repmat([-1 1],[1 4]);
%epiphasedir = repmat([-1 -1 -1  1 1 1],[1 2]);
% *** EDIT ME ***
epiphasedir = -1;  % this is for coronal-looking slices, indicating that the phase-encode goes upwards (F >> H)
epiphasedir = 1;  % axial, towards the bottom

% what is the total readout time in milliseconds for an EPI slice?
% (note that 'Private_0043_102c' in the dicominfo of the EPI files gives the time per phase-encode line in microseconds.
% I confirmed that this time is correct by checking against the waveforms displayed by plotter.)
% *** EDIT ME ***
epireadouttime = 0.67 * (104/2);  % divide by 2 if you are using 2x acceleration
epireadouttime = 0.999974 * (162/3);  % divide by 2 if you are using 2x acceleration
epireadouttime = 0.66 * (120/2);  % divide by 2 if you are using 2x acceleration
  % *** consider using Keith's dicom_readout_msec.m to check the readout time value?

% what fieldmap should be used for each EPI run? ([] indicates default behavior, which is to attempt
% to match fieldmaps to EPI runs 1-to-1, or if there is only one fieldmap, apply that fieldmap
% to all EPI runs, or if there is one more fieldmap than EPI runs, interpolate each successive
% pair of fieldmaps; see preprocessfmri.m for details.)
  %epifieldmapasst = {1 [1 2] [2 3] 3};    %%[splitmatrix([1:12; 2:13]',1) {4 4 4 4}];
  %epifieldmapasst = {1 1 [1 2] [2 3]  [3 4] [4 5]  [5 6] [6 7]};
% *** EDIT ME ***
if isempty(fieldmaptweaks{sessix})
  epifieldmapasst = splitmatrix([1:length(epifilenames); 2:length(epifilenames)+1]',1);
else
  epifieldmapasst = fieldmaptweaks{sessix}{2};
end

% how many volumes should we ignore at the beginning of each EPI run?
numepiignore = 0;

% what volume should we use as reference in motion correction? ([] indicates default behavior which is
% to use the first volume of the first run; see preprocessfmri.m for details.  set to NaN if you
% want to omit motion correction.)
motionreference = [1 1];  %[1 137];

% for which volumes should we ignore the motion parameter estimates?  this should be a cell vector
% of the same length as the number of runs.  each element should be a vector of indices, referring
% to the volumes (after dropping volumes according to <numepiignore>).  can also be a single vector
% of indices, in which case we use that for all runs.  for volumes for which we ignore the motion
% parameter estimates, we automatically inherit the motion parameter estimates of the closest
% volumes (if there is a tie, we just take the mean).  [] indicates default behavior which is to 
% do nothing special.
epiignoremcvol = [];

% by default, we tend to use double format for computation.  but if memory is an issue,
% you can try setting <dformat> to 'single', and this may reduce memory usage.
dformat = 'single';

% apply Gaussian spatial smoothing? 3-element vector indicating FWHM in mm.
% [] means do nothing.
epismoothfwhm = [];  % [3 3 3]

% what cut-off frequency should we use for filtering motion parameter estimates? ([] indicates default behavior
% which is to low-pass filter at 1/90 Hz; see preprocessfmri.m for details.)
  %motioncutoff = [];
motioncutoff = Inf;




% what extra transformation should we use in the final resampling step? ([] indicates do not perform an extra transformation.)
% tr = maketransformation([0 0 0],[1 2 3],[117.040207900659 115.623739214828 43.7254395779036],[1 2 3],[-0.670461404521411 0.704835952293662 -0.249037864439899],[116 116 42],[232 232 84],[1 1 1],[0 0 0],[0 0 0],[0 0 0]);
% extratrans = transformationtomatrix(tr,0,[2 2 2]);





% construct coordinates desired in the final space (original voxel-size separation)
crop0 = nsdcropranges{nsdsubjixfun(sessix)};
  %[xx,yy,zz] = ndgrid(1:1/1.8:120,1:1/1.8:120,1:1/1.8:84);
[xx,yy,zz] = ndgrid(crop0{1}(1):crop0{1}(2), ...
                    crop0{2}(1):crop0{2}(2), ...
                    crop0{3}(1):crop0{3}(2));
C = [flatten(xx); flatten(yy); flatten(zz)];
C(4,:) = 1;

% load T (this matrix maps from session gradunwarped EPI to reference gradunwarped EPI).
% for the first session, T should be the identity matrix.
T = loadmulti([ppdir '/preprocessVER1_sessioncoregistration/alignment.mat'],'T');

% map coordinates from reference EPI space to session EPI space
C = inv(T)*C;

% load warpcoords (this tells us where in original to sample to fix it)
warpcoords = loadmulti([ppdir '/preprocessVER1_gradunwarp/warpcoords.mat'],'warpcoords');

% determine new warp coordinates
Cnew = [];
for p=1:3
  Cnew(p,:) = ba_interp3_wrapper(reshape(warpcoords(p,:),[120 120 84]),C(1:3,:),'cubic');
end

% prepare input
extratrans = {reshape(Cnew(1,:),size(xx)) ...
              reshape(Cnew(2,:),size(yy)) ...
              reshape(Cnew(3,:),size(zz))};

% final setting
finalepisizeOVERRIDE = [1.8 1.8 1.8];







% what is the desired resolution for the resampled volumes? ([] indicates to just use the original EPI resolution.)
targetres = [];

% should we perform slice shifting?  if so, specify band-pass filtering cutoffs in Hz, like [1/360 1/20].
% probably should be left as [] which means to do nothing special.
sliceshiftband = [];

% these are constants that are used in fmriquality.m.  it is probably 
% fine to leave this as [], which means to use default values.
% NaN means to skip the fmriquality calculations.
fmriqualityparams = [];

% what kind of time interpolation should we use on the fieldmaps (if applicable)?
% ([] indicates to use the default, which is cubic interpolation.)
fieldmaptimeinterp = 'linear';

% should we use a binary 3D ellipse mask in the motion parameter estimation?
% if [], do nothing special (i.e. do not use a mask).
% if {}, then we will prompt the user to interactively determine the
%   3D ellipse mask (see defineellipse3d.m for details).  upon completion,
%   the parameters will be reported to the command window so that you can
%   simply supply those parameters if you run again (so as to avoid user interaction).
% if {MN SD}, then these will be the parameters that determine the mask to be used.
% *** EDIT ME ***
%mcmask = {};
%mcmask = { [0.46 0.5 0.14] [0.3 0.38 0.76];};
mcmask = loadmulti([ppdir '/mcmask.mat'],'mcmask');

% how should we handle voxels that have NaN values after preprocessing?
% if [], we use the default behavior which is to zero out all voxels that have a NaN
% value at any point in the EPI data.  see preprocessfmri.m for other options.
maskoutnans = [];

% savefile:  what .nii files (accepting a 1-indexed integer) should we save the final EPI data to?
%            (in the special EPI flattening case, we save the data to raw binary files (time x voxels) instead of .nii files.)
% savefileB: what .nii file should we save the valid voxels (binary mask) to?  ([] means do not save.)
% savefileC: what .nii file should we save the mean volume to?  ([] means do not save.)
% savefileD: what .nii file should we save the mad volume to?  ([] means do not save.)
% savefileE: what .nii file should we save the tsnr volume to?  ([] means do not save.)
% (we automatically make parent directories if necessary.)
savefile = [ppdir '/preprocessVER1SECONDSTAGELOW/run%02d.nii'];
savefileB = [ppdir '/preprocessVER1SECONDSTAGELOW/valid.nii'];
savefileC = [ppdir '/preprocessVER1SECONDSTAGELOW/mean.nii'];
savefileD = [ppdir '/preprocessVER1SECONDSTAGELOW/mad.nii'];
savefileE = [ppdir '/preprocessVER1SECONDSTAGELOW/tsnr.nii'];

% what .txt file should we keep a diary in?
diaryfile = [ppdir '/preprocessVER1SECONDSTAGELOW/diary.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW:

  mkdirquiet(stripfile(diaryfile));
  diary(diaryfile);
preprocessfmri_standard;
  diary off;
