function preprocess_nsd_applybrainmask(subjix)

% function preprocess_nsd_applybrainmask(subjix)
%
% <subjix> is 1-8
%
% Apply brainmask to the T1s, T2s, SWIs, TOFs (at all resolutions).
% We save out files named like _masked.nii.gz

% setup
nsdsetup;
outputdir = sprintf('%s/ppdata/%s-structurals',nsddir,allnsdids{subjix});

% change directory
olddir = pwd;
cd(outputdir);

% do it
res = {'0pt5' '0pt8' '1pt0'};
typs = {'T1' 'T2' 'SWI' 'TOF'};
nums = 0:10;  % 0 is fully averaged; 1-10 are the individual reps
for p=1:length(res)
  for q=1:length(typs)
    for r=1:length(nums)
      if nums(r)==0
        typ0 = [typs{q}];
      else
        typ0 = [typs{q} sprintf('_rep%02d',nums(r))];
      end
      file0 = sprintf('%s_%s.nii.gz',typ0,res{p});
      if ~exist(file0,'file')
        continue;
      end
      unix_wrapper(sprintf('fslmaths %s_%s.nii.gz -mas brainmask_%s.nii.gz %s_%s_masked.nii.gz -odt input', ...
                           typ0,res{p},res{p},typ0,res{p}));
    end
  end
end

% revert
cd(olddir);




% % do it
% res = {'0pt5' '0pt8' '1pt0'};
% typs = {'T1'};
% for p=1:length(res)
%   unix_wrapper(sprintf('fslmaths brainmask_%s.nii.gz -thr 0.5 -bin brainmask_%s_thresh.nii.gz',res{p},res{p}));
%   for q=1:length(typs)
%     unix_wrapper(sprintf('fslmaths %s_%s.nii.gz -mas brainmask_%s_thresh.nii.gz %s_%s_masked.nii.gz', ...
%                          typs{q},res{p},res{p},typs{q},res{p}));
%   end
% end
% 

