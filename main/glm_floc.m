function glm_floc(subjix,ppdirname,tr,runix,stimversion,glmdirname)

% approach: assume HRF, condition split, no boots

% setup
nsdsetup;
ppdir = regexprep(nsdscreeningfun2(subjix),'rawdata','ppdata');
glmdir = regexprep(ppdir,'ppdata',glmdirname);
mkdirquiet(glmdir);
stimdur = 4;
reps_per_run = 6;

% load data
data = {};
for p=1:length(runix)
  a1 = load_untouch_nii(sprintf('%s/%s/run%02d.nii',ppdir,ppdirname,runix(p)));
  data{p} = single(a1.img);
end
clear a1;

% load design matrix
stimulus = loadmulti(sprintf('%s/designmatrix%s_floc.mat',ppdir,stimversion),'stimulus');

% check sanity
assert(length(data)==length(stimulus));

% truncate data to match stimulus preparation
for p=1:length(data)
  size(data{p})
  size(stimulus{p})
  assert(size(data{p},4) >= size(stimulus{p},1));
  data{p} = data{p}(:,:,:,1:size(stimulus{p},1));
end

% perform condition-split
[stimulus,reps_per_run] = condition_split(stimulus,reps_per_run);

% run GLMdenoise
figdir = sprintf('%s/GLMdenoise_flocFULL_figures',glmdir);
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
save(sprintf('%s/GLMdenoise_flocFULL.mat',glmdir),'-struct','results','-v7.3');
