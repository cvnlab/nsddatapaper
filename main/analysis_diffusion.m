% setup
nsdsetup;

% find all DWI preprocessed results
files0 = matchfiles(sprintf('~/nsd/nsddata_diffusion/ppdata/subj*/run*/dwi/dwi.nii.gz'));
assert(length(files0)==8*2);

% load in the data
allvols = {};
for zz=1:length(files0)
  a1 = load_untouch_nii(files0{zz})  % img: [320x320x320x98 single] for Run 1

  % average across b0 images
  ix = [1 2 18 34 50 66 82];  % b0 images only
  data = mean(a1.img(:,:,:,ix),4);
  
  % record
  allvols{zz} = data;
end

% average across run 1 and run 2 and write out NIFTIs (int16 format)
for zz=1:8
  outputdir = sprintf('~/Dropbox/nsddata/ppdata/subj%02d/anat/',zz);
  vol = (allvols{(zz-1)*2+1} + allvols{(zz-1)*2+2})/2;

  % save
  nsd_savenifti(int16(vol),   [.8 .8 .8],sprintf('%s/DWI_0pt8.nii.gz',outputdir));

  % convert to other resolutions
    % 0.5
  [newvol,newvolsize] = changevolumeres(double(vol),[.8 .8 .8],[512 512 512]);
  nsd_savenifti(int16(newvol),[.5 .5 .5],sprintf('%s/DWI_0pt5.nii.gz',outputdir));
    % 1.0
  [newvol,newvolsize] = changevolumeres(double(vol),[.8 .8 .8],[256 256 256]);
  nsd_savenifti(int16(newvol),[1 1 1],   sprintf('%s/DWI_1pt0.nii.gz',outputdir));

end
