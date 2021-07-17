%% %%%%%%%% Determine alignment between HRT2 and T2 and warp labels

for subjix=1:8

  % setup
  nsdsetup;
  outputdir = sprintf('%s/ppdata/%s-structurals',nsddir,allnsdids{subjix});

  % define
  file1 = sprintf('~/nsddata/ppdata/subj%02d/anat/T2_0pt5_masked.nii.gz',subjix);  % official T2 0.5-mm masked volume
  file2 = sprintf('~/nsddata/ppdata/subj%02d/anat/HRT2/HRT2_raw.nii.gz',subjix);   % raw HRT2 volume
  file3 = sprintf('~/nsd/ppdata/MTL/subj%02d_MTL.nii.gz',subjix);                  % manually labeled volume
  pfile = sprintf('%s/HRT2alignment.mat',outputdir);
  
  % copy over the manual ROI file to official location
  copyfile(file3,sprintf('~/nsddata/ppdata/subj%02d/anat/HRT2/MTL_rawlabels.nii.gz',subjix));

  % load vol1
  vol1orig = load_untouch_nii(file1);
  vol1size = vol1orig.hdr.dime.pixdim(2:4);
  vol1 = double(vol1orig.img);
  vol1(isnan(vol1)) = 0;
  fprintf('vol1 has dimensions %s at %s mm.\n',mat2str(size(vol1)),mat2str(vol1size));

  % load vol2
  vol2orig = load_untouch_nii(file2);
  vol2size = vol2orig.hdr.dime.pixdim(2:4);
  vol2 = double(vol2orig.img);
  vol2(isnan(vol2)) = 0;
  fprintf('vol2 has dimensions %s at %s mm.\n',mat2str(size(vol2)),mat2str(vol2size));

  % load ROI
  a3 = load_untouch_nii(file3);

  % homogenize volumes (for purposes of alignment)
  [vol1h,brainmask,polymodel] = homogenizevolumes(vol1);
  [vol2h,brainmask,polymodel] = homogenizevolumes(vol2);

  %%%%% DO ALIGNMENT

  % get the alignment started with an initial seed
  if exist(pfile,'file')
    trinitial = loadmulti(pfile,'trinitial');
    alignvolumedata(vol1h,vol1size,vol2h,vol2size,trinitial);
  else
    alignvolumedata(vol1h,vol1size,vol2h,vol2size);

    % allow some manual alignment (sx should be -1)
    keyboard;
    trinitial = alignvolumedata_exporttransformation;
    save(pfile,'trinitial');
  end
  
  % compute a bounding box
  mask0 = double(a3.img)>0;  % all labeled voxels
  [d1,d2,d3,ii] = computebrickandindices(mask0);
  fct = 0.1;  % extra 10% on each side
  d1mn = max(1,round(d1(1)-range(d1)*fct));
  d1mx = min(size(mask0,1),round(d1(end)+range(d1)*fct));
  d2mn = max(1,round(d2(1)-range(d2)*fct));
  d2mx = min(size(mask0,2),round(d2(end)+range(d2)*fct));
  d3mn = max(1,round(d3(1)-range(d3)*fct));
  d3mx = min(size(mask0,3),round(d3(end)+range(d3)*fct));
  mask0 = zeros(size(mask0));
  mask0(d1mn:d1mx,d2mn:d2mx,d3mn:d3mx) = 1;

  % affine transformation using the bounding box; correlation metric
  mn = mask0;
  sd = [];
  alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[4 4 2]);
  alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[2 2 1]);
  alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[2 2 1]);
  alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
  alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
  alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
  alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
  alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
  alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
  
  %%%%% GET OUTPUTS

  % record transformation
  tr = alignvolumedata_exporttransformation;

  % clean up
  close all;

  % convert the transformation to a matrix (from HRT2 to T2)
  T = transformationtomatrix(tr,0,vol1size);

  % save parameters
  save(pfile,'tr','T','-append');

  % get slices from T2 to match HRT2
  match1 = extractslices(vol1,vol1size,vol2,vol2size,tr,0);

  %%%%% INSPECT
  
  % define some stuff (RIA and then CCW rotation means left->right in the image is R->L)
  skip = 4;
  
  % calc some
  vol2M = vol2 .* mask0;
  match1M = match1 .* mask0;
  
  % write images
  imwrite(uint8(255*makeimagestack(rotatematrix(vol2(:,:,1:skip:end),1,2,1),1)), gray(256),             sprintf('~/nsddata/inspections/HRT2/subj%02d_raw.png',subjix));
  imwrite(uint8(255*makeimagestack(rotatematrix(vol2M(:,:,1:skip:end),1,2,1),1)),gray(256),             sprintf('~/nsddata/inspections/HRT2/subj%02d_masked.png',subjix));
  imwrite(uint8(255*makeimagestack(rotatematrix(double(a3.img(:,:,1:skip:end)),1,2,1),[0 20])),jet(256),sprintf('~/nsddata/inspections/HRT2/subj%02d_rawlabels.png',subjix));
  imwrite(uint8(255*makeimagestack(rotatematrix(mask0(:,:,1:skip:end),1,2,1),[0 1])),gray(256),         sprintf('~/nsddata/inspections/HRT2/subj%02d_mask.png',subjix));
  imwrite(uint8(255*makeimagestack(rotatematrix(match1(:,:,1:skip:end),1,2,1),0.1)),gray(256),          sprintf('~/nsddata/inspections/HRT2/subj%02d_T2matched.png',subjix));
  imwrite(uint8(255*makeimagestack(rotatematrix(match1M(:,:,1:skip:end),1,2,1),0.1)),gray(256),         sprintf('~/nsddata/inspections/HRT2/subj%02d_T2matched_masked.png',subjix));

  % save NIFTIs
  voltemp = vol2orig;
  voltemp.img = int16(mask0);
  save_untouch_nii(voltemp,sprintf('~/nsddata/ppdata/subj%02d/anat/HRT2/HRT2_mask.nii.gz',subjix));  % mask
  voltemp.img = int16(match1);
  save_untouch_nii(voltemp,sprintf('~/nsddata/ppdata/subj%02d/anat/HRT2/T2matched.nii.gz',subjix));  % match from T2
  
  %%%%% WARP LABELS
  
  % slice through ROI labels to bring them to 0.5mm T2 space, using winner-take-all
  labelvol = [];
  for labelix=0:20
    labelvol(:,:,:,labelix+1) = extractslices(vol1,vol1size,double(double(a3.img)==labelix),vol2size,tr,1,'linear');
  end
  [~,iix] = max(labelvol,[],4);
  iix = iix-1;
  
  % save lh
  temp = iix;
  temp(mod(temp,2)==1) = 0;  % set odd to 0
  temp = temp/2;
  nsd_savenifti(int16(temp),[.5 .5 .5],sprintf('~/nsd/ppdata/MTL/subj%02d_MTLinT2_lh.nii.gz',subjix));

  % save rh
  temp = iix;
  temp(mod(temp,2)==0) = 0;  % set even to 0
  temp = ceil(temp/2);
  nsd_savenifti(int16(temp),[.5 .5 .5],sprintf('~/nsd/ppdata/MTL/subj%02d_MTLinT2_rh.nii.gz',subjix));

end

%% %%%%%% NOTES:

% Notes:
% - Anisotropy can produce "jagged" displays
% - The edge effects of the first or last slice can lead to weird visualizations
% - 1.5mm is fairly discrete...

%% %%%%%% POST-PROCESS 

% see external_subcortical.m => IMPORT AND VET
