function preprocess_nsd_structuralinspection(subjix)

% function preprocess_nsd_structuralinspection(subjix)
%
% <subjix> is 1-8
% 
% write out some orthMID inspections of the T1, T2, EPI, SWI, TOF, DWI.

% setup
nsdsetup;
outputdir = sprintf('%s/ppdata/%s-structurals',nsddir,allnsdids{subjix});
outputdir2 = sprintf('%s/inspection',outputdir);
mkdirquiet(outputdir2);

% write out orthMID figure for inspection
files = {'T1_1pt0_masked.nii.gz' 'T2_1pt0_masked.nii.gz' 'EPIaffine_Warped.nii.gz' 'EPIsyn_Warped.nii.gz' 'SWI_1pt0_masked.nii.gz' 'TOF_1pt0_masked.nii.gz' 'DWI_1pt0.nii.gz'};
prefixes = {'T1' 'T2' 'EPIaffine' 'EPIsyn' 'SWI' 'TOF' 'DWI'};
for zz=1:length(files)
  file0 = sprintf('%s/%s',outputdir,files{zz});
  if ~exist(file0,'file')
    warning('%s does not exist',file0);
    continue;
  end
  a1 = load_untouch_nii(file0);
  vol = double(a1.img);
  mn = nanmean(vol(:));
  mx = max(size(vol));
  im1 = placematrix(mn*ones(mx,mx),rotatematrix(vol(:,:,round(end/2)),1,2,1),[1 1]);
  im2 = placematrix(mn*ones(mx,mx),rotatematrix(squish(vol(:,round(end/2),:),2),1,2,1),[1 1]);
  im3 = placematrix(mn*ones(mx,mx),rotatematrix(squish(vol(round(end/2),:,:),2),1,2,1),[1 1]);
  im = cat(3,im1,im2,im3);
  imwrite(uint8(255*makeimagestack(im,1)),gray(256),sprintf('%s/%s.png',outputdir2,prefixes{zz}));
end
