% this was just used internally (for QC and for subject leaderboard)
% and to derive FIR for library purposes.

function glm_nsdSIMPLE(ppdir)

% this is one-regressor GLM for the screening session phase.

% approach: every stimulus is the same, assume HRF, no boots

% define
stimdur = 3;
tr = 1;

% figure out runs
files0 = matchfiles(sprintf('%s/preprocessVER1/run*.nii',ppdir));
if length(files0)==14
  files0 = files0(2:end-1);  % resting-state has two extra runs
else
  assert(length(files0)==12);
end

% load data
data = {};
for p=1:length(files0)
  a1 = load_untouch_nii(files0{p});
  data{p} = single(a1.img);
end
clear a1;

% load design matrix
stimulus = loadmulti(sprintf('%s/designmatrixSINGLETRIAL_nsd.mat',ppdir),'stimulus');

% check sanity
assert(length(data)==length(stimulus));

% truncate data
for p=1:length(data)
  assert(size(data{p},4) >= size(stimulus{p},1));
  data{p} = data{p}(:,:,:,1:size(stimulus{p},1));
end

% make a single condition
stimulus = cellfun(@(x) sum(x,2),stimulus,'UniformOutput',0);

%%%%%

% run GLMdenoise
figdir = sprintf('%s/GLMdenoise_nsdSIMPLE_figures',ppdir);
rmdirquiet(figdir);
results = GLMdenoisedata(stimulus,data,stimdur,tr, ...
                         'assume',[],struct('numpcstotry',20,'numboots',0), ...
                         figdir);

% massage output
results = rmfield(results,{'models' 'modelse'});
md = results.modelmd;
results.modelmd = md{2};
results.modelmd_hrf = md{1};
clear md;

% save
save(sprintf('%s/GLMdenoise_nsdSIMPLE.mat',ppdir),'-struct','results','-v7.3');
clear results;

%%%%%

% run GLMestimatemodel for FIR
results1 = GLMestimatemodel(stimulus,         data,         stimdur,tr,'fir',30,0);
results2 = GLMestimatemodel(stimulus(1:2:end),data(1:2:end),stimdur,tr,'fir',30,0);
results3 = GLMestimatemodel(stimulus(2:2:end),data(2:2:end),stimdur,tr,'fir',30,0);
firs = cat(4,results1.modelmd,results2.modelmd,results3.modelmd);
R2 = cat(4,results1.R2,results2.R2,results3.R2);
clear results1 results2 results3;
save(sprintf('%s/GLMdenoise_nsdSIMPLEFIR.mat',ppdir),'firs','R2');

% figures
outdir0 = sprintf('%s/GLMdenoise_nsdSIMPLEFIR_figures',ppdir);
mkdirquiet(outdir0);
imwrite(uint8(255*makeimagestack(signedarraypower(R2(:,:,:,1)/100,0.5),[0 1])),hot(256), ...
  sprintf('%s/R2.png',outdir0));
