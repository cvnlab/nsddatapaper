function glm_nsdimagery2(ppdir,ppdirname,designfile,hrfbasisfile,tr,glmdirname,hrfmanifoldfile)


% notes:
% - essentially unchanged wrt glm_nsd...
% - 3-s events for vis, imag, and att.
% - att does have two adjacent 3-s trials. it's going to be noisy.


glmdir = regexprep(ppdir,'ppdata',glmdirname);
mkdirquiet(glmdir);

%%%%% LOAD

% define
stimdur = 3;

% figure out runs
files0 = matchfiles(sprintf('%s/%s/run*.nii',ppdir,ppdirname));
assert(length(files0)==12);

% load data
data = {};
for p=1:length(files0)
  a1 = load_untouch_nii(files0{p});
  data{p} = single(a1.img);
end
clear a1;
data

% load design matrix
stimulus = loadmulti(sprintf('%s/%s',ppdir,designfile),'stimulus');
stimulus

% check sanity
assert(length(data)==length(stimulus));

% truncate data
for p=1:length(data)
  assert(size(data{p},4) >= size(stimulus{p},1));
  data{p} = data{p}(:,:,:,1:size(stimulus{p},1));
end

% do it
[results,resultsdesign] = GLMestimatesingletrial(stimulus,data,stimdur,tr,{sprintf('%s/GLMsingleoutputs',glmdir) sprintf('%s/GLMsinglefigures',glmdir)});
