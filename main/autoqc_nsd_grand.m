function autoqc_nsd_grand(subjix,glmdirname)

% setup
nsdsetup;

% get all scan directories
alldirs0 = [nsdscreeningfun2(subjix) nsdallsessfun2(subjix)];

%%%%% get mean and valid from every session

% init
meanvols  = single([]);
validvols = single([]);

% loop over directories
for p=1:length(alldirs0)

  % calc
  ppdir = regexprep(alldirs0{p},'rawdata','ppdata');
  outputdir = sprintf('%s/../%s-%s-glmBASIC',ppdir,allnsdids{subjix},glmdirname);

  % load mean and valid
  if isequal(glmdirname,'glmdata')
    a1 = load_untouch_nii([ppdir '/preprocessVER1SECONDSTAGE/mean.nii']);
    a2 = load_untouch_nii([ppdir '/preprocessVER1SECONDSTAGE/valid.nii']);
  elseif isequal(glmdirname,'glmdataLOW')
    a1 = load_untouch_nii([ppdir '/preprocessVER1SECONDSTAGELOW/mean.nii']);
    a2 = load_untouch_nii([ppdir '/preprocessVER1SECONDSTAGELOW/valid.nii']);
  end

  % record
  meanvols(:,:,:,p) =  single(a1.img);
  validvols(:,:,:,p) = single(a2.img);

end

%%%%% compute correlation matrix (homogenize the mean volumes and then restrict to valid voxels)

% compute grand
grandmean =  mean(meanvols,4);
grandvalid = sum(validvols,4);

% compute homogenized volumes
homvols = homogenizevolumes(meanvols,[99 1/4 5 3]);

% define good
good = find(all(validvols==1,4));
  %%%% should be restricted to ellipse or bet or something???

% compute confusion
cmatrix = calcconfusionmatrix(subscript(squish(homvols,3),{good ':'}),[],2);

%%%%%

% make output dirs
mkdirquiet(outputdir);

% write niftis
a1.img = int16(grandmean);
save_untouch_nii(a1,[outputdir '/grandmean.nii']);
a2.img = int16(grandvalid);
save_untouch_nii(a2,[outputdir '/grandvalid.nii']);

% write inspections
imwrite(uint8(255*makeimagestack(grandmean,1)), gray(256),sprintf('%s/grandmean.png',outputdir));
imwrite(uint8(255*makeimagestack(grandvalid,-2)),jet(256),sprintf('%s/grandvalid.png',outputdir));

% visualize
figureprep([100 100 500 500]); hold on;
imagesc(cmatrix,[.9 1]);
colormap(jet);
colorbar;
axis image tight;
xlabel('Session number');
ylabel('Session number');
set(gca,'XTick',1:size(cmatrix,2),'XTickLabel',0:size(cmatrix,2)-1);
set(gca,'YTick',1:size(cmatrix,1),'YTickLabel',0:size(cmatrix,1)-1);
set(gca,'YDir','reverse');
title('Pairwise correlation of mean volume');
figurewrite('pairwise',[],[],outputdir);

% save
save([outputdir '/grand.mat'],'cmatrix');
