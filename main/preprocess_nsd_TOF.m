function preprocess_nsd_TOF(subjix)

% function preprocess_nsd_TOF(subjix)
% 
% <subjix> is 1-8
%
% Use ANTS to use a BSplineSyN to match the TOF image
% to the T1 1.0-mm resolution. We use some carefully crafted parameters.
%
% We then create a final 0.5-mm TOF image and some other resolutions (0.8 mm, 1.0 mm).
%
% Note: We use the non-gradunwarped version of the TOF.

% setup
nsdsetup;
outputdir = sprintf('%s/ppdata/%s-structurals',nsddir,allnsdids{subjix});

% find file
files0 = matchfiles(sprintf('%s/*TOF*1.nii.gz',outputdir));  % get the non-gradunwarped one!
assert(length(files0)==1);

% dampen the most bright voxels (this causes confusion in the registration)
a1 = load_untouch_nii(files0{1});
temp = double(a1.img);
temp(temp>80) = 50;  % NOTE!!!
a1.img = cast(temp,class(a1.img));
b = tempname;
file1b = [b '.nii.gz'];
save_untouch_nii(a1,file1b);
  %%figure;imagesc(makeimagestack(double(temp(:,:,1:30:end)))); drawnow;

% define
file0 = sprintf('%s/T1_1pt0_masked.nii.gz',outputdir);          % target
file1 = files0{1};                                              % source (moving)
file2 = sprintf('%s/brainmask_1pt0.nii.gz',outputdir);          % mask (on the target)
file3 = sprintf('%s/T1_0pt5_masked.nii.gz',outputdir);          % for the final transformed image

% determine registration using syn (NOTICE: 200)
a = 'antsRegistration --dimensionality 3 --float 0 --output [%s,%sWarped.nii.gz] --interpolation Linear --winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 --initial-moving-transform [%s,%s,1] --transform Rigid[0.1] --metric MI[%s,%s,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[0.1] --metric MI[%s,%s,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform BSplineSyN[0.1,200,0,3] --metric CC[%s,%s,1,4] --convergence [100x70x50x20,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --masks [%s]';
pre0 = sprintf('%s/TOFsyn_',outputdir);
b = tempname;
savetext([b '.sh'],sprintf(a,pre0,pre0,file0,file1b,file0,file1b,file0,file1b,file0,file1b,file2));
unix_wrapper(sprintf('sh %s',[b '.sh']));

% use the warp to create final version (0.5-mm)
a = 'antsApplyTransforms --dimensionality 3 --input %s --reference-image %s -t %s/%s -t %s/%s --output %s/%s --interpolation BSpline';
b = tempname;
savetext([b '.sh'],sprintf(a,file1,file3,outputdir,'TOFsyn_1Warp.nii.gz',outputdir,'TOFsyn_0GenericAffine.mat',outputdir,'TOF_0pt5.nii.gz'));
unix_wrapper(sprintf('sh %s',[b '.sh']));

% change double to int16 format (mangle in place), multiply by 10!
file4 = sprintf('%s/TOF_0pt5.nii.gz',outputdir);
a = load_untouch_nii(file4);
nsd_savenifti(int16(10*double(a.img)),[.5 .5 .5],file4);

% convert to other resolutions
a1 = load_untouch_nii(file4);
  % 0.8
[newvol,newvolsize] = changevolumeres(double(a1.img),[.5 .5 .5],[320 320 320]);
nsd_savenifti(int16(newvol),[.8 .8 .8],sprintf('%s/TOF_0pt8.nii.gz',outputdir));
  % 1.0
[newvol,newvolsize] = changevolumeres(double(a1.img),[.5 .5 .5],[256 256 256]);
nsd_savenifti(int16(newvol),[1 1 1],   sprintf('%s/TOF_1pt0.nii.gz',outputdir));

% clean up
delete(file1b);
