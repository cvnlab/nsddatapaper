% Deal with merging fMRI data from partial sessions.

cd /home/surly-raid1/kendrick-data/nsd/ppdata/-NSD929-nsd02/  % DATE CENSORED

%%%%% preprocessVER1

a1 = load_untouch_nii('preprocessVER1/run04.nii');
a2 = load_untouch_nii('preprocessVER1/run05.nii');

% 50 trials from run04a; ignore first 3 blank trials from run04b
a1.img = cat(4,a1.img(:,:,:,1:50*4),a2.img(:,:,:,3*4+1:end));
a1.img = a1.img(:,:,:,1:301);
a1.hdr.dime.dim(5) = 301;

% deal with scaling
[~,~,polymodel1] = homogenizevolumes(mean(double(a1.img(:,:,:,1:200)),4));
[~,~,polymodel2] = homogenizevolumes(mean(double(a1.img(:,:,:,201:end)),4));
pp = polymodel1./polymodel2;
pp(pp>2) = 2;
a1.img(:,:,:,201:end) = bsxfun(@times,double(a1.img(:,:,:,201:end)),pp);

save_untouch_nii(a1,'preprocessVER1/run04.nii');
for p=6:13
  assert(movefile(sprintf('preprocessVER1/run%02d.nii',p),sprintf('preprocessVER1/run%02d.nii',p-1)));
end
