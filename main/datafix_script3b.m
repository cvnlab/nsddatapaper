% Deal with merging fMRI data from partial sessions.

cd /home/surly-raid1/kendrick-data/nsd/ppdata/NSD929-nsd02/  % DATE CENSORED

%%%%% preprocessVER1SECONDSTAGE

a1 = load_untouch_nii('preprocessVER1SECONDSTAGE/run04.nii');
a2 = load_untouch_nii('preprocessVER1SECONDSTAGE/run05.nii');

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

save_untouch_nii(a1,'preprocessVER1SECONDSTAGE/run04.nii');
for p=6:13
  assert(movefile(sprintf('preprocessVER1SECONDSTAGE/run%02d.nii',p),sprintf('preprocessVER1SECONDSTAGE/run%02d.nii',p-1)));
end

%%%%% preprocessVER1SECONDSTAGELOW

a1 = load_untouch_nii('preprocessVER1SECONDSTAGELOW/run04.nii');
a2 = load_untouch_nii('preprocessVER1SECONDSTAGELOW/run05.nii');

% 50 trials from run04a; ignore first 3 blank trials from run04b
a1.img = cat(4,a1.img(:,:,:,1:(50*4/(4/3))),a2.img(:,:,:,(3*4)/(4/3)+1:end));
a1.img = a1.img(:,:,:,1:226);
a1.hdr.dime.dim(5) = 226;

% deal with scaling
[~,~,polymodel1] = homogenizevolumes(mean(double(a1.img(:,:,:,1:(50*4/(4/3)))),4));
[~,~,polymodel2] = homogenizevolumes(mean(double(a1.img(:,:,:,(3*4)/(4/3)+1:end)),4));
pp = polymodel1./polymodel2;
pp(pp>2) = 2;
a1.img(:,:,:,(3*4)/(4/3)+1:end) = bsxfun(@times,double(a1.img(:,:,:,(3*4)/(4/3)+1:end)),pp);

save_untouch_nii(a1,'preprocessVER1SECONDSTAGELOW/run04.nii');
for p=6:13
  assert(movefile(sprintf('preprocessVER1SECONDSTAGELOW/run%02d.nii',p),sprintf('preprocessVER1SECONDSTAGELOW/run%02d.nii',p-1)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALSO, FIX MOTION EXPORTED FILES

% for the exported motion tsv files, we do (1mm)
take first 200 from run04
start from 13th from run05
and we then have 301.
then adjust all of the run numbers.

% for the exported motion tsv files, we do (1pt8mm)
take first 150 from run04
start from 10th from run05
and we then have 226.
then adjust all of the run numbers.
