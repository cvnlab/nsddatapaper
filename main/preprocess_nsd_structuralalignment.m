function preprocess_nsd_structuralalignment(subjix,mode,wantpause)

% function preprocess_nsd_structuralalignment(subjix,mode,wantpause)
%
% <subjix> is 1-8
% <mode> is 1 (T1) or 2 (T2)
% <wantpause> is whether to pause for getting a better initial seed
%
% we do substantial work regarding preparing anatomical volumes:
% - reject any bad anatomical volumes
% - coregister volumes together (T1 and T2)
% - deal with gradunwarp coordinates
% - upsample when preparing the final volumes
% - create multiple resolutions of the results
%
% NOTE: there is an error in the gradunwarp handling; see below.
 
%%%%%%%%%% SETUP

% input
if ~exist('wantpause','var') || isempty(wantpause)
  wantpause = 0;
end

% setup
nsdsetup;
outputdir = sprintf('%s/ppdata/%s-structurals',nsddir,allnsdids{subjix});
mkdirquiet(sprintf('%s/mode%d',outputdir,mode));

% figure out files to process
switch mode
case 1
  
  % files to align together
  nifti2files = matchfiles(sprintf('%s/*T1w*_gradunwarped.nii.gz',outputdir));
  nifti2fileswarp = matchfiles(sprintf('%s/*T1w*_warp.nii.gz',outputdir));
  
  % reject any?
  ix = setdiff(1:length(nifti2files),rejectt1{subjix});
  nifti2files = nifti2files(ix);
  nifti2fileswarp = nifti2fileswarp(ix);

  % final target to get to
  nifti1file = nifti2files{1};
  
  % final filename
  outfileprefix = 'T1';
  
  % alignment flavor
  amode = 1;

case 2
  
  % files to align together
  nifti2files = matchfiles(sprintf('%s/*T2w*_gradunwarped.nii.gz',outputdir));
  nifti2fileswarp = matchfiles(sprintf('%s/*T2w*_warp.nii.gz',outputdir));
  
  % reject any?
  ix = setdiff(1:length(nifti2files),rejectt2{subjix});
  nifti2files = nifti2files(ix);
  nifti2fileswarp = nifti2fileswarp(ix);

  % final target to get to
  nifti1file = sprintf('%s/T1-originalspace.nii.gz',outputdir);
  
  % final filename
  outfileprefix = 'T2';
  
  % alignment flavor
  amode = 2;

end

%%%%%%%%%% ALIGN EACH TO THE TARGET

% init
regularvol = single(0);
resampledvol = single(0);
Ts = [];

% load nifti1
vol1orig = load_untouch_nii(nifti1file);
vol1size = vol1orig.hdr.dime.pixdim(2:4);
vol1 = double(vol1orig.img);
vol1(isnan(vol1)) = 0;
fprintf('vol1 has dimensions %s at %s mm.\n',mat2str(size(vol1)),mat2str(vol1size));

% homogenize the volume
[vol1h,brainmask,polymodel] = homogenizevolumes(vol1);

% prepare coords as where, relative to original matrix space, to sample from
% in order to achieve the desired matrix size and resolution
indices = resamplingindices(1,256,-2);  % 0.5 mm over 256 mm FOV
[xx,yy,zz] = ndgrid(indices,indices,indices);
xx = (xx - (1+256)/2) / vol1size(1);  % center, convert to voxel indices
yy = (yy - (1+256)/2) / vol1size(2);
zz = (zz - (1+256)/2) / vol1size(3);
xx = xx + (1+size(vol1,1))/2;
yy = yy + (1+size(vol1,2))/2;
zz = zz + (1+size(vol1,3))/2;
origsz = size(xx);
coords = [flatten(xx); flatten(yy); flatten(zz)];
coords(4,:) = 1;  % 4 x N
clear xx yy zz;

% do it
tr = [];
for p=1:length(nifti2files)

  %%%%% LOAD

  % load nifti2
  vol2orig = load_untouch_nii(nifti2files{p});
  vol2size = vol2orig.hdr.dime.pixdim(2:4);
  vol2 = double(vol2orig.img);
  vol2(isnan(vol2)) = 0;
  fprintf('vol2 has dimensions %s at %s mm.\n',mat2str(size(vol2)),mat2str(vol2size));

  % homogenize the volume
  [vol2h,brainmask,polymodel] = homogenizevolumes(vol2);

  %%%%% DO ALIGNMENT

  % start the alignment using homogenized volumes
  % [note: vol2 is "reference"; vol1 is "target"]
  % [this is to ensure we interpolate through the new session]
  alignvolumedata(vol2h,vol2size,vol1h,vol1size,tr);

  % allow some manual alignment?
  if wantpause && isempty(tr)
    keyboard;
  end

  % define ellipse or load it
  if ismember(mode,[1 2])
    file0 = sprintf('%s/mode1_ellipse.mat',outputdir);
%    file0 = sprintf('%s/mode%d_ellipse.mat',outputdir,mode);
  end
  if exist(file0,'file')
    load(file0,'mn','sd');
  else
    [f,mn,sd] = defineellipse3d(vol1h);  % note on the "target" [this is because auto expects that]
    save(file0,'mn','sd');
  end

  % auto-align
  switch amode

  case 1
    % auto-align (correlation, RIGID-BODY, 'linear' interpolation)
    alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[4 4 4]);
    alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[2 2 2]);
    alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);

  case 2
    % auto-align (mutual info, RIGID-BODY, 'linear' interpolation)
    alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[4 4 4],[],[],[],1);
    alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[2 2 2],[],[],[],1);
    alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1],[],[],[],1);

  end

        %   % auto-align (correlation, AFFINE, 'linear' interpolation, many iterations to ensure convergence)
        %   alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[4 4 4]);
        %   alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[2 2 2]);
        %   alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);   % just that?
        %   alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
        %   alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);

%     alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[4 4 4],[],[],[],1);
%     alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[4 4 4],[],[],[],1);
%     alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[2 2 2],[],[],[],1);
%     alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[2 2 2],[],[],[],1);
%     alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1],[],[],[],1);
%     alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1],[],[],[],1);
%     alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1],[],[],[],1);
%     alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1],[],[],[],1);
%     alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1],[],[],[],1);
%     alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1],[],[],[],1);

  %%%%% GET OUTPUTS

  % record transformation
  tr = alignvolumedata_exporttransformation;

  % clean up
  close all;

  % convert the transformation to a matrix (from vol2 to vol1)
  T = transformationtomatrix(tr,1,vol2size);

  % get slices from vol2 to match vol1
  match1 =  extractslices(vol2,vol2size,vol1,vol1size,tr,0);
%  match1h = extractslices(vol2h,vol2size,vol1h,vol1size,tr,0);
  regularvol = regularvol + match1 / length(nifti2files);

  % write out orthMID figure for inspection
  mn = nanmean(match1(:));
  mx = max(size(match1));
  im1 = placematrix(mn*ones(mx,mx),rotatematrix(match1(:,:,round(end/2)),1,2,1),[1 1]);
  im2 = placematrix(mn*ones(mx,mx),rotatematrix(squish(match1(:,round(end/2),:),2),1,2,1),[1 1]);
  im3 = placematrix(mn*ones(mx,mx),rotatematrix(squish(match1(round(end/2),:,:),2),1,2,1),[1 1]);
  im = cat(3,im1,im2,im3);
  imwrite(uint8(255*makeimagestack(im,1)),gray(256),sprintf('%s/mode%d/orthMID%02d.png',outputdir,mode,p));

  % record T
  Ts(:,:,p) = T;

  %%%%% PROCEED TO RESAMPLE AT FINE SCALE

  % load warp coordinates (this tells us where to sample in the original volume to get the unwarped volume)
  a1 = load_untouch_nii(nifti2fileswarp{p});
  warpcoords = double(a1.img);  % X x Y x Z x 3; where to sample from ... seems to be mm units
  for zz=1:3
    warpcoords(:,:,:,zz) = warpcoords(:,:,:,zz) / a1.hdr.dime.pixdim(1+zz) + ...
                            a1.hdr.dime.pixdim(1+zz)/2;  % convert from mm to matrix units
  end
  
  
  % transform coordinates to the matrix space of the source volume
  coords2 = inv(T)*coords;  % 4 x N
  
  % use warpcoords to figure out where to sample in the original volume
  coords3 = [];  % 3 x N
  for zz=1:3
    coords3(zz,:) = ba_interp3_wrapper(warpcoords(:,:,:,zz),coords2(1:3,:),'cubic');
  end
  clear coords2;

  % go ahead and resample the volume
  voltemp = reshape(ba_interp3_wrapper(vol2,coords3(1:3,:),'cubic'),origsz);
  resampledvol = resampledvol + voltemp / length(nifti2files);
  clear coords3;
  
  %%%% WARNING: the above step is an error. The above interpolation should have
  %%%% resampled through the non-gradunwarped volume (the original volume).
  %%%% The consequence of the error is that we have overcompensated for the
  %%%% gradwarp effect. This is ultimately a very small spatial effect,
  %%%% so the results have NOT been revisited and re-corrected.

  % for individual volumes, save a 0.8 version
  voltempB = changevolumeres(flipdim(nanreplace(double(voltemp),0,3),1),[.5 .5 .5],[320 320 320]);
  nsd_savenifti(int16(voltempB),[.8 .8 .8],sprintf('%s/%s_rep%02d_0pt8.nii.gz',outputdir,outfileprefix,p));

end

%%%%%%%%%% FINISH UP

% save parameters
save(sprintf('%s/mode%d_alignment.mat',outputdir,mode),'Ts');

% save original space version
vol1orig.img = regularvol;  % this is single format, I think
file0 = sprintf('%s/%s-originalspace.nii',outputdir,outfileprefix);
save_untouch_nii(vol1orig,file0); gzip(file0); delete(file0);

% save fancy version
volfun = @(x) flipdim(x,1);   % we want to save in LPI!
sz = sizefull(resampledvol,3);
file0 = sprintf('%s/%s_0pt5.nii',outputdir,outfileprefix);
mkdirquiet(stripfile(file0));
save_nii(make_nii(volfun(int16(resampledvol)),[.5 .5 .5],([1 1 1]+sz)/2),file0);  % notice int16 format
gzip(file0); delete(file0);

% convert to other resolutions
a1 = load_untouch_nii([file0 '.gz']);
  % 0.8
[newvol,newvolsize] = changevolumeres(double(a1.img),[.5 .5 .5],[320 320 320]);
assert(isequal(newvolsize,[.8 .8 .8]));
sz = sizefull(newvol,3);
file1 = sprintf('%s/%s_0pt8.nii',outputdir,outfileprefix);
save_nii(make_nii(int16(newvol),[.8 .8 .8],([1 1 1]+sz)/2),file1);  % notice int16 format
gzip(file1); delete(file1);
  % 1.0
[newvol,newvolsize] = changevolumeres(double(a1.img),[.5 .5 .5],[256 256 256]);
assert(isequal(newvolsize,[1 1 1]));
sz = sizefull(newvol,3);
file1 = sprintf('%s/%s_1pt0.nii',outputdir,outfileprefix);
save_nii(make_nii(int16(newvol),[1 1 1],([1 1 1]+sz)/2),file1);  % notice int16 format
gzip(file1); delete(file1);
