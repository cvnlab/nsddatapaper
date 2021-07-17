function preprocess_nsd_MNIandbrainmask(subjix)

% function preprocess_nsd_MNIandbrainmask(subjix)
%
% <subjix> is 1-8
%
% Starting with the prepared 1.0-mm T1, determine MNI alignment using fnirt.
% On basis of the determined alignment, we backproject to the native subject space
% a liberal 1.0-mm brain mask. We also create other resolutions for the
% brain mask.

% # Based on
% # https://github.com/Washington-University/HCPpipelines/blob/master/PreFreeSurfer/scripts/BrainExtraction_FNIRTbased.sh

% setup
nsdsetup;
outputdir = sprintf('%s/ppdata/%s-structurals',nsddir,allnsdids{subjix});

% define
Input = sprintf('%s/T1_1pt0.nii.gz',outputdir);                         % master T1 (1mm)
Ref = sprintf('%s/templates/MNI152_T1_1mm.nii.gz',nsddatadir);          % MNI template
RefMask = sprintf('%s/templates/MNI152_T1_1mm_brain_mask_dil.nii.gz',nsddatadir);  % MNI liberal brain mask
RefMask2 = sprintf('%s/templates/MNI152_T1_1mm_brain_mask_dil_dilM.nii.gz',nsddatadir);
WD = sprintf('%s/MNI',outputdir);
BaseName = 'T1';
OutputBrainMask = sprintf('%s/brainmask_1pt0.nii.gz',outputdir);        % output file
ConfigFile = sprintf('%s/templates/T1_2_MNI152_2mm.cnf',nsddatadir);    % fnirt config

% make dir
mkdirquiet(WD);

% affine to MNI
unix_wrapper(sprintf('flirt -interp spline -dof 12 -in "%s" -ref "%s" -omat "%s"/roughlin.mat -out "%s"/"%s"_to_MNI_roughlin.nii.gz -nosearch',Input,Ref,WD,WD,BaseName));

% nonlinear to MNI
unix_wrapper(sprintf('fnirt --in="%s" --ref="%s" --aff="%s"/roughlin.mat --refmask="%s" --fout="%s"/str2standard.nii.gz --jout="%s"/NonlinearRegJacobians.nii.gz --refout="%s"/IntensityModulatedT1.nii.gz --iout="%s"/"%s"_to_MNI_nonlin.nii.gz --logout="%s"/NonlinearReg.txt --intout="%s"/NonlinearIntensities.nii.gz --cout="%s"/NonlinearReg.nii.gz --config="%s"',Input,Ref,WD,RefMask,WD,WD,WD,WD,BaseName,WD,WD,WD,ConfigFile));
  % str2standard.nii.gz is the warpfield?
  % T1_to_MNI_nonlin.nii.gz - resampled T1 (linear interp) [but, alternatively, could use applywarp...]

% invert the warp
unix_wrapper(sprintf('invwarp --ref="%s" -w "%s"/str2standard.nii.gz -o "%s"/standard2str.nii.gz',Input,WD,WD));  % NOTE: Should be Input, not Ref!

% use nearest-neighbor to get brain mask in subject native space
unix_wrapper(sprintf('applywarp --rel --interp=nn --in="%s" --ref="%s" -w "%s"/standard2str.nii.gz -o "%s"',RefMask2,Input,WD,OutputBrainMask));

% convert to other resolutions
a1 = load_untouch_nii(gunziptemp(OutputBrainMask));
  % 0.8
[newvol,newvolsize] = changevolumeres(double(a1.img),[1 1 1],[320 320 320]);
assert(isequal(newvolsize,[.8 .8 .8]));
sz = sizefull(newvol,3);
file1 = sprintf('%s/brainmask_0pt8.nii',outputdir);
save_nii(make_nii(int16(newvol>0.5),[.8 .8 .8],([1 1 1]+sz)/2),file1);  % notice int16 format
gzip(file1); delete(file1);
  % 0.5
[newvol,newvolsize] = changevolumeres(double(a1.img),[1 1 1],[512 512 512]);
assert(isequal(newvolsize,[.5 .5 .5]));
sz = sizefull(newvol,3);
file1 = sprintf('%s/brainmask_0pt5.nii',outputdir);
save_nii(make_nii(int16(newvol>0.5),[.5 .5 .5],([1 1 1]+sz)/2),file1);  % notice int16 format
gzip(file1); delete(file1);




% OLD:
% 
% % convert to 0.8
% file1 = sprintf('%s/brainmask_0pt8.nii.gz',outputdir);
% unix_wrapper(sprintf('mri_convert --voxsize 0.8 0.8 0.8 %s %s',OutputBrainMask,file1));
% 
% % convert to 0.5
% file1 = sprintf('%s/brainmask_0pt5.nii.gz',outputdir);
% unix_wrapper(sprintf('mri_convert --voxsize 0.5 0.5 0.5 %s %s',OutputBrainMask,file1));
