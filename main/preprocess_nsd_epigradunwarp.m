function preprocess_nsd_epigradunwarp(datadir)

% must be run on stone only!

%%%%%%%%%% SETUP

% internal constants
gradfile = '7tps';

% calc
ppdir = regexprep(datadir,'rawdata','ppdata');

%%%%%%%%%% MAKE A NIFTI FROM FIRST VOLUME

% find first EPI volume
epidirs = matchfiles([datadir '/dicom/*bold*']);
epifiles = matchfiles([epidirs{1} '/*']);

% copy first EPI volume to a temp DICOM file
dir0 = maketempdir;  % like /tmp/FMG/
file0 = [dir0 'temp.dcm'];
assert(copyfile(epifiles{1},file0));

% convert dicom to NIFTI
result = unix_wrapper(sprintf('dcm2nii -r N -x N %s',file0));
temp = regexp(result,'GZip\.\.\.(.+)\n','tokens');
file1 = [dir0 temp{1}{1}];  % this is the NIFTI file

% run Keith's fslreorient2std
unix_wrapper(sprintf('fslreorient2std_inplace %s',file1));

% rename the file to something easier
file2 = [dir0 'temp.nii.gz'];
assert(movefile(file1,file2));

%%%%%%%%%% CALL GRADUNWARP

% call gradunwarp
prefix0 = [dir0 'temp'];
unix_wrapper(sprintf('gradunwarp -w %s_warp.nii.gz -m %s_mask.nii.gz %s.nii.gz %s_gradunwarped.nii.gz %s',prefix0,prefix0,prefix0,prefix0,gradfile));
  % temp.nii.gz is the original volume
  % temp_mask.nii.gz is some weird valid mask
  % temp_gradunwarped.nii.gz is the unwarped volume
  % temp_warp.nii.gz is three volumes with the x- y- and z- warp coordinates

%%%%%%%%%% LOAD AND PROCESS WARPCOORDS

% load warp coordinates (this tells us where to sample in the original volume to get the unwarped volume)
a1 = load_untouch_nii(gunziptemp(sprintf('%s_warp.nii.gz',prefix0)));
warpcoords = a1.img;  % X x Y x Z x 3; where to sample from ... seems to be mm units
for p=1:3
  warpcoords(:,:,:,p) = warpcoords(:,:,:,p) / a1.hdr.dime.pixdim(1+p) + ...
                          a1.hdr.dime.pixdim(1+p)/2;  % convert from mm to matrix units
end

% NOTE: this has to be done to go from our internal mean.nii space to the .nii.gz space:
%   a1.img = flipdim(permute(a2.img,[2 1 3]),2);
%     OR, EQUIVALENTLY:
%   a1.img = permute(flipdim(a2.img,1),[2 1 3]);
% let's construct a T matrix that does that.
T1 = [-1  0 0 a1.hdr.dime.dim(2)+1;
       0  1 0 0;
       0  0 1 0;
       0  0 0 1];
T2 = [0 1 0 0;
      1 0 0 0;
      0 0 1 0;
      0 0 0 1];
T = T2*T1;

% first, we need to adjust the orientation of the warpcoords
% volume itself to go from .nii.gz space to internal mean.nii space.
warpcoords = flipdim(permute(warpcoords,[2 1 3 4]),1);

% next, we need to adjust the coordinates in warpcoords so that
% it tells us where to sample in our internal mean.nii space
warpcoords = squish(warpcoords,3)';   % 3 x XYZ
warpcoords(4,:) = 1;                  % 4 x XYZ
warpcoords = inv(T)*warpcoords;       % 4 x XYZ; now it is relative to internal mean.nii space

%%%%%%%%%% UNWARP MEAN.NII

% load mean pp volume
a1 = load_untouch_nii([ppdir '/preprocessVER1/mean.nii']);

% go ahead and unwarp the volume
newvol = reshape(ba_interp3_wrapper(double(a1.img),warpcoords(1:3,:),'cubic'),size(a1.img));

%%%%%%%%%% SAVE RESULTS

% create output directory
mkdirquiet([ppdir '/preprocessVER1_gradunwarp']);

% save warp coordinates
save([ppdir '/preprocessVER1_gradunwarp/warpcoords.mat'],'warpcoords');

% save unwarped volume
a1.img = cast(newvol,class(a1.img));
save_untouch_nii(a1,[ppdir '/preprocessVER1_gradunwarp/mean.nii']);

% clean up
rmdirquiet(dir0);
