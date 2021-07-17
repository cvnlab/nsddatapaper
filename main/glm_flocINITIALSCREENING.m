function glm_flocINITIALSCREENING(ppdir,runix)

% this is first-pass GLM for the screening session phase.

% approach: assume HRF, condition split, no boots

% define
stimdur = 4;
tr = 1;
reps_per_run = 6;

% load data
data = {};
for p=1:length(runix)
  a1 = load_untouch_nii(sprintf('%s/preprocessVER1/run%02d.nii',ppdir,runix(p)));
  data{p} = single(a1.img);
end
clear a1;

% load design matrix
stimulus = loadmulti(sprintf('%s/designmatrix_floc.mat',ppdir),'stimulus');

% check sanity
assert(length(data)==length(stimulus));

% truncate data
for p=1:length(data)
  assert(size(data{p},4) >= size(stimulus{p},1));
  data{p} = data{p}(:,:,:,1:size(stimulus{p},1));
end

% perform condition-split
[stimulus,reps_per_run] = condition_split(stimulus,reps_per_run);

% run GLMdenoise
figdir = sprintf('%s/GLMdenoise_floc_figures',ppdir);
rmdirquiet(figdir);
results = GLMdenoisedata(stimulus,data,stimdur,tr, ...
                         'assume',[],struct('numpcstotry',20,'numboots',0), ...
                         figdir);

% massage output
results = rmfield(results,{'models' 'modelse'});
md = results.modelmd;
results.modelmd = md{2};
results.modelmd_hrf = md{1};
results.reps_per_run = reps_per_run;

% save
save(sprintf('%s/GLMdenoise_floc.mat',ppdir),'-struct','results','-v7.3');
