function preprocess_nsd_functionaltostructuralalignment(subjix)

% function preprocess_nsd_functionaltostructuralalignment(subjix)
% 
% <subjix> is 1-8
%
% Use ANTS to determine an EPIaffine and then an EPIsyn
% to match the EPI (mean of first 5 sessions 1.0-mm resolution)
% to the T2 (1.0-mm version). We use some carefully crafted parameters.

% setup
nsdsetup;
outputdir = sprintf('%s/ppdata/%s-structurals',nsddir,allnsdids{subjix});

% bet
unix_wrapper(sprintf('bet %s/T2_1pt0.nii.gz %s/T2_1pt0_betmasked.nii.gz -m',outputdir,outputdir));

% define
file0 = sprintf('%s/T2_1pt0_masked.nii.gz',outputdir);
file1 = sprintf('%s/ppdata/subj%02d/func1mm/meanFIRST5.nii.gz',nsddatadir,subjix);
file2 = sprintf('%s/T2_1pt0_betmasked_mask.nii.gz',outputdir);

% affine
a = 'antsRegistration --dimensionality 3 --float 0 --output [%s,%sWarped.nii.gz] --interpolation Linear --winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 --initial-moving-transform [%s,%s,1] --transform Rigid[0.1] --metric MI[%s,%s,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[0.1] --metric MI[%s,%s,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --masks [%s]';
pre0 = sprintf('%s/EPIaffine_',outputdir);
b = tempname;
savetext([b '.sh'],sprintf(a,pre0,pre0,file0,file1,file0,file1,file0,file1,file2));
unix_wrapper(sprintf('sh %s',[b '.sh']));

% syn
a = 'antsRegistration --dimensionality 3 --float 0 --output [%s,%sWarped.nii.gz] --interpolation Linear --winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 --initial-moving-transform [%s,%s,1] --transform Rigid[0.1] --metric MI[%s,%s,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[0.1] --metric MI[%s,%s,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform BSplineSyN[0.1,400,0,3] --metric CC[%s,%s,1,4] --convergence [100x70x50x20,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --masks [%s]';
pre0 = sprintf('%s/EPIsyn_',outputdir);
b = tempname;
savetext([b '.sh'],sprintf(a,pre0,pre0,file0,file1,file0,file1,file0,file1,file0,file1,file2));
unix_wrapper(sprintf('sh %s',[b '.sh']));
