% this is an interactive script to determine crop ranges.

% download mean.nii locally

cd ~/Desktop

% load gradunwarp results
a1 = load_untouch_nii('mean.nii');

vol = double(a1.img);

imagesc(max(vol,[],3).^.4); axis image;
imagesc(squeeze(max(vol,[],2)).^.4); axis image;

temp = {[12 110] [20 104] [4 83]};

% test it
vol = vol(temp{1}(1):temp{1}(2),temp{2}(1):temp{2}(2),temp{3}(1):temp{3}(2));
