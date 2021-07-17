% divide T2 by mean EPI, and divide by an automatically determined threshold
% (within brain voxels) such that 1 is meaningful. this "signaldropout"
% volume is useful, showing what voxels have high signal intensities in
% the T2 but are corrupted by signal dropout in the EPI.

% define
reses = {[1 1 1] [1.8 1.8 1.8]};
trs = [1 4/3];
vers = {'func1mm' 'func1pt8mm'};

% do it
for ver=1:length(vers)
  for p=1:8

    % load
    a1 = load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj%02d/%s/mean.nii.gz',p,vers{ver}));
    a2 = load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj%02d/%s/T2_to_%s.nii.gz',p,vers{ver},vers{ver}));
    a3 = load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj%02d/%s/aseg.nii.gz',p,vers{ver}));
    
    % divide (T2/EPI)
    data = zerodiv(double(a2.img),double(a1.img));  % if divide by 0, leave at 0

    % compute mask (aseg ~= 0)
    mask = double(a3.img)~=0;
    thresh = findtailthreshold(data(mask))
    
    % normalize units
    data = data ./ thresh;  % now, 1 means the threshold
    
    % save signaldropout volume
    outputfile = sprintf('~/nsd/nsddata/ppdata/subj%02d/%s/signaldropout.nii.gz',p,vers{ver});
    nsd_savenifti(data,reses{ver},outputfile,trs(ver));

    % save masked volume (to help visualization)
    outputfile = sprintf('~/nsd/nsddata/ppdata/subj%02d/%s/signaldropout_masked.nii.gz',p,vers{ver});
    nsd_savenifti(copymatrix(data,~mask,0),reses{ver},outputfile,trs(ver));
    
  end
end

% notes:
% - we tried and VERY similar to manually drawn dropout regions.. but this shouldn't be surprising. everything is sane.
