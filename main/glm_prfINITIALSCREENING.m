function glm_prfINITIALSCREENING(ppdir,runix)

% approach: assume HRF, average the 3 reps up front, fit all data, do supergrid.

% define
tr = 1;

% load data
data = {};
for p=1:length(runix)
  a1 = load_untouch_nii(sprintf('%s/preprocessVER1/run%02d.nii',ppdir,runix(p)));
  data{p} = single(a1.img);
end
clear a1;

% load stimulus
a1 = load('/home/surly-raid4/kendrick-data/nsd/nsddata/stimuli/prf/RETBARsmall.mat');
a2 = load('/home/surly-raid4/kendrick-data/nsd/nsddata/stimuli/prf/RETWEDGERINGMASHsmall.mat');
%stimulus = {a1.stim a2.stim};
stimulus = repmat({a1.stim a2.stim},[1 3]);
clear a1 a2;

% prep data
%data{1} = (data{1}+data{3}+data{5})/3;
%data{2} = (data{2}+data{4}+data{6})/3;
%data = data(1:2);

% truncate data
for p=1:length(data)
  assert(size(data{p},4) >= size(stimulus{p},3));
  data{p} = data{p}(:,:,:,1:size(stimulus{p},3));
end

% start parallel MATLAB to speed up execution
if isempty(gcp)
  parpool;
end

% analyze
results = analyzePRF(stimulus,data,tr,struct('seedmode',-2));

% what is the R2 that seems to separate noise from signal?
R2thresh = findtailthreshold(results.R2(:),0);

% define a mask
bad = ~(results.R2 > R2thresh);

% visualize
figdir = sprintf('%s/analyzePRF_prf_figures',ppdir);
rmdirquiet(figdir); mkdirquiet(figdir);
imwrite(uint8(255*makeimagestack(results.meanvol,1)),                     gray(256),sprintf('%s/meanvol.png',figdir));
imwrite(uint8(255*makeimagestack(results.R2,[0 1])),hot(256), sprintf('%s/R2.png',figdir));
temp = uint8(255*cmaplookup(makeimagestack(results.ang),0,360,1,hsv(256)));
imwrite(temp,                                                           sprintf('%s/ang.png',figdir));
imwrite(copymatrix(temp,logical(repmat(makeimagestack(bad),[1 1 3])),0),sprintf('%s/angMASKED.png',figdir));
temp = uint8(255*cmaplookup(makeimagestack(results.ecc),0,100,0,jet(256)));
imwrite(temp,                                                           sprintf('%s/ecc.png',figdir));
imwrite(copymatrix(temp,logical(repmat(makeimagestack(bad),[1 1 3])),0),sprintf('%s/eccMASKED.png',figdir));

% save
save(sprintf('%s/analyzePRF_prf.mat',ppdir),'-struct','results','-v7.3');
save(sprintf('%s/analyzePRF_prf.mat',ppdir),'R2thresh','-append');
