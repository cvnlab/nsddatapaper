function cvnlookupwrite(rgbimg,rawimg,common,ext,wantomit)

% function cvnlookupwrite(rgbimg,rawimg,common,ext,wantomit)
%
% <rgbimg>,<rawimg> are from cvnlookup
% <common> is a filename path
% <ext> is the specific file name
% <wantomit> (optional) is whether to omit the .nii.gz
%
% imwrite the <rgbimg> as .png and save <rawimg> as a surface .nii.gz.
% files are like <common>_<ext>.png.

% inputs
if ~exist('wantomit','var') || isempty(wantomit)
  wantomit = 0;
end

imwrite(rgbimg,sprintf('%s_%s.png',common,ext));
if ~wantomit
  nsd_savenifti(flipdim(rotatematrix(rawimg,1,2,1),2),[1 1 1],sprintf('%s_%s.nii.gz',common,ext));
end
