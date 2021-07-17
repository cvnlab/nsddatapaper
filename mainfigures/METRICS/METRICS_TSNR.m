% this script generates, for all subjects, a depiction of tsnr results for one run.
% this includes tsnr related stuff (i.e. mean and valid).

% define
massiveim = [];
for subjix=1:8

  sessix = 20;
  runix = 1;

  % setup
  dirs0 = nsdallsessfun2(subjix);
  ppdir = regexprep(dirs0{sessix},'rawdata','ppdata');

  % load EPI run
  epidir = matchfiles([dirs0{sessix} '/dicom/*cmrr_mbep2d_bold*']);  % first EPI run
  epidir = epidir{1};
  [vols,volsizes,inplanematrixsizes,trs] = dicomloaddir(epidir);

  % load tsnr stuff
  a1 = load(sprintf('%s/autoqc_fmri.mat',ppdir),'tsnr','tsnrvol','mnvol');
  tsnr = a1.tsnr(1,:);  % raw volumes, quadratic detrend, mean/std, median over voxels with valid intensity

  % check the value
  fprintf('subj = %d, %.3f\n',subjix,tsnr(runix));

  % pull out the tsnr volume
  vol = a1.tsnrvol(:,:,:,1);

  % recompute valid
  valid = a1.mnvol{1} > 0.1 * prctile(flatten(a1.mnvol{1}),99);

  %% PROCEED TO VISUALIZE

  % notes:
  %   matlab matrix: dim1 is A->P, dim2 is R->L, dim3 is I->S

  rng = [0 prctile(vol(:),[99])]   % [0 81.71082]
  rng = [0 80];
  vol0 = vol; cmap0 = jet(256);
  im1a = cmaplookup(makeimagestack(vol0(:,:,round(end/2)),rng,0),0,1,[],cmap0);
  im2a = cmaplookup(makeimagestack(rotatematrix(squish(vol0(:,round(end/2),:),2),1,2,1),rng,0),0,1,[],cmap0);
  im3a = cmaplookup(makeimagestack(rotatematrix(squish(vol0(round(end/2),:,:),2),1,2,1),rng,0),0,1,[],cmap0);

  figureprep([100 100 200 200]);
  imagesc(randn(10,10),rng);
  colormap(jet(256));
  colorbar;
  figurewrite('tsnrlegend',[],-1,'~/Dropbox/KKTEMP');

  rng = [0 1];
  vol0 = valid; cmap0 = gray(256);
  im1b = cmaplookup(makeimagestack(vol0(:,:,round(end/2)),rng,0),0,1,[],cmap0);
  im2b = cmaplookup(makeimagestack(rotatematrix(squish(vol0(:,round(end/2),:),2),1,2,1),rng,0),0,1,[],cmap0);
  im3b = cmaplookup(makeimagestack(rotatematrix(squish(vol0(round(end/2),:,:),2),1,2,1),rng,0),0,1,[],cmap0);

  epivol = vols{1}(:,:,:,round(end/2));
  rng = [0 prctile(epivol(:),[99])];
  vol0 = epivol; cmap0 = gray(256);
  im1c = cmaplookup(makeimagestack(vol0(:,:,round(end/2)),rng,0),0,1,[],cmap0);
  im2c = cmaplookup(makeimagestack(rotatematrix(squish(vol0(:,round(end/2),:),2),1,2,1),rng,0),0,1,[],cmap0);
  im3c = cmaplookup(makeimagestack(rotatematrix(squish(vol0(round(end/2),:,:),2),1,2,1),rng,0),0,1,[],cmap0);

  massiveim = cat(2,massiveim,uint8(255*cat(1,im1c,im1a,im1b,im2c,im2a,im2b,im3c,im3a,im3b)));

end
imwrite(massiveim,'~/Dropbox/KKTEMP/allsubjects.png');

%   imwrite(uint8(255*cat(1,im1a,im1b,im1c)),sprintf('~/Dropbox/KKTEMP/subj%02d_1.png',subjix));
%   imwrite(uint8(255*cat(1,im2a,im2b,im2c)),sprintf('~/Dropbox/KKTEMP/subj%02d_2.png',subjix));
%   imwrite(uint8(255*cat(1,im3a,im3b,im3c)),sprintf('~/Dropbox/KKTEMP/subj%02d_3.png',subjix));

%   imwrite(uint8(255*makeimagestack(vol0(:,:,round(end/2)),rng,0)),cmap0,                              sprintf('~/Dropbox/KKTEMP/subj%02d_tsnr1.png',subjix))
%   imwrite(uint8(255*makeimagestack(rotatematrix(squish(vol0(:,round(end/2),:),2),1,2,1),rng,0)),cmap0,sprintf('~/Dropbox/KKTEMP/subj%02d_tsnr2.png',subjix))
%   imwrite(uint8(255*makeimagestack(rotatematrix(squish(vol0(round(end/2),:,:),2),1,2,1),rng,0)),cmap0,sprintf('~/Dropbox/KKTEMP/subj%02d_tsnr3.png',subjix))
% 
%   imwrite(uint8(255*makeimagestack(vol0(:,:,round(end/2)),rng,0)),cmap0,                              sprintf('~/Dropbox/KKTEMP/subj%02d_valid1.png',subjix))
%   imwrite(uint8(255*makeimagestack(rotatematrix(squish(vol0(:,round(end/2),:),2),1,2,1),rng,0)),cmap0,sprintf('~/Dropbox/KKTEMP/subj%02d_valid2.png',subjix))
%   imwrite(uint8(255*makeimagestack(rotatematrix(squish(vol0(round(end/2),:,:),2),1,2,1),rng,0)),cmap0,sprintf('~/Dropbox/KKTEMP/subj%02d_valid3.png',subjix))
% 
%   imwrite(uint8(255*makeimagestack(vol0(:,:,round(end/2)),rng,0)),cmap0,                              sprintf('~/Dropbox/KKTEMP/subj%02d_epi1.png',subjix))
%   imwrite(uint8(255*makeimagestack(rotatematrix(squish(vol0(:,round(end/2),:),2),1,2,1),rng,0)),cmap0,sprintf('~/Dropbox/KKTEMP/subj%02d_epi2.png',subjix))
%   imwrite(uint8(255*makeimagestack(rotatematrix(squish(vol0(round(end/2),:,:),2),1,2,1),rng,0)),cmap0,sprintf('~/Dropbox/KKTEMP/subj%02d_epi3.png',subjix))

%%%%%%%%%%%%%%%%%%%%%

% for example subject in main figure, choose subj02. manually crop.
