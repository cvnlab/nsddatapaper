function preprocess_nsd_epialignment(sessix)


% idea is to get as close as you can with simple low-dimensional (and principled/controlled) linear methods...



%%%%% SETUP

% setup
nsdsetup;
datadir = datadirs{sessix};
ppdir = regexprep(datadir,'rawdata','ppdata');
outputdir = sprintf('%s/../%s-sessioncoregistration',ppdir,nsdidfun(sessix));
outputdir2 = sprintf('%s/preprocessVER1_sessioncoregistration',ppdir);

% make output dirs
mkdirquiet(outputdir);
mkdirquiet(outputdir2);

% inputs
nifti2 = [ppdirs{sessix} '/preprocessVER1_gradunwarp/mean.nii'];   % the thing to align

% figure out reference
cfiles = matchfiles([outputdir '/cumvol*.nii']);
cumcnt = length(cfiles)+1;  % figure out what count we are on
if cumcnt==1
  nifti1 = nifti2;       % if there is no existing cumvol, this means the reference is the thing to align
else
  nifti1 = cfiles{end};  % otherwise, use the last cumvol
end

%%%%% LOAD

% load nifti1
vol1orig = load_untouch_nii(gunziptemp(nifti1));
vol1size = vol1orig.hdr.dime.pixdim(2:4);
vol1 = double(vol1orig.img);
vol1(isnan(vol1)) = 0;
fprintf('vol1 has dimensions %s at %s mm.\n',mat2str(size(vol1)),mat2str(vol1size));

% load nifti2
vol2orig = load_untouch_nii(gunziptemp(nifti2));
vol2size = vol2orig.hdr.dime.pixdim(2:4);
vol2 = double(vol2orig.img);
vol2(isnan(vol2)) = 0;
fprintf('vol2 has dimensions %s at %s mm.\n',mat2str(size(vol2)),mat2str(vol2size));

% homogenize the volumes
[vol1h,brainmask,polymodel] = homogenizevolumes(vol1);
[vol2h,brainmask,polymodel] = homogenizevolumes(vol2);

%%%%% DO ALIGNMENT

% start the alignment using homogenized volumes
% [note: vol2 is "reference"; vol1 is "target"]
% [this is to ensure we interpolate through the new session]
alignvolumedata(vol2h,vol2size,vol1h,vol1size);   %,tr

% define ellipse or load it
if exist([outputdir '/ellipse.mat'],'file')
  load([outputdir '/ellipse.mat'],'mn','sd');
else
  [f,mn,sd] = defineellipse3d(vol1h);  % note on the "target" [this is because auto expects that]
  save([outputdir '/ellipse.mat'],'mn','sd');
end

% auto-align (correlation, AFFINE, 'linear' interpolation, many iterations to ensure convergence)
alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[4 4 4]);
alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[2 2 2]);
alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);
alignvolumedata_auto(mn,sd,[1 1 1 1 1 1 0 0 0 0 0 0],[1 1 1]);
alignvolumedata_auto(mn,sd,[0 0 0 0 0 0 1 1 1 1 1 1],[1 1 1]);

%%%%% GET OUTPUTS

% record transformation
tr = alignvolumedata_exporttransformation;

% convert the transformation to a matrix
T = transformationtomatrix(tr,1,vol2size);

% get slices from vol2 to match vol1
match1 =  extractslices(vol2,vol2size,vol1,vol1size,tr,0);
match1h = extractslices(vol2h,vol2size,vol1h,vol1size,tr,0);

%%%%% SAVE

% calculate and save cumulative volume [to be used for subsequent alignment] [not the homogenized volume!]
cumvol = cast((vol1*(cumcnt-1) + match1)/cumcnt,class(vol1orig.img));
vol1orig.img = cumvol;
save_untouch_nii(vol1orig,sprintf([outputdir '/cumvol%02d.nii'],cumcnt));

% write out inspection figures
imwrite(uint8(255*makeimagestack(double(cumvol),1)), gray(256),sprintf([outputdir '/cumulative%02d.png'],cumcnt));
imwrite(uint8(255*makeimagestack(double(match1),1)), gray(256),sprintf([outputdir '/match%02d.png'],cumcnt));
imwrite(uint8(255*makeimagestack(double(match1h),1)),gray(256),sprintf([outputdir '/matchhom%02d.png'],cumcnt));

% write out orthMID figure
mn = nanmean(match1(:));
mx = max(size(match1));
im1 = placematrix(mn*ones(mx,mx),match1(:,:,round(end/2)),[1 1]);
im2 = placematrix(mn*ones(mx,mx),rotatematrix(squish(match1(:,round(end/2),:),2),1,2,1),[1 1]);
im3 = placematrix(mn*ones(mx,mx),rotatematrix(squish(match1(round(end/2),:,:),2),1,2,1),[1 1]);
im = cat(3,im1,im2,im3);
imwrite(uint8(255*makeimagestack(im,1)),gray(256),sprintf([outputdir '/orthMID%02d.png'],cumcnt));

% save alignment results to a .mat file
save(sprintf('%s/alignment.mat',outputdir2),'tr','T');
