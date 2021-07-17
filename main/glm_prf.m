function glm_prf(subjix,ppdirname,tr,runix,prepname,stimversion,glmdirname,chunkstart)

% notes:
% - assume HRF
% - average the 3 reps of each stimulus type up front
% - fit ALL voxels within brainmask. very expensive.
% - we use a chunking mechanism to do it in chunks.
% - use supergrid to set the initial seed
% - run types are: 93/94 x 3, 8.4° stimulus diameter
%   93=bars
%   94=wedge/ring
% - the model attempts to fit the pRF exponent!
% - we do not use seedmode 0 1; we use only seedmode 2.
% - we limit iterations to 'maxiter',100

% inputs
if ~exist('chunkstart','var') || isempty(chunkstart)
  chunkstart = 1;
end

% setup
nsdsetup;
ppdir = regexprep(nsdscreeningfun2(subjix),'rawdata','ppdata');
glmdir = regexprep(ppdir,'ppdata',glmdirname);
mkdirquiet(glmdir);

% load data
data = {};
for p=1:length(runix)
  a1 = load_untouch_nii(sprintf('%s/%s/run%02d.nii',ppdir,ppdirname,runix(p)));
  data{p} = single(a1.img);
end
clear a1;

% load brainmask (it is in LPI; we are in ARI)
a1 = load_untouch_nii(gunziptemp(sprintf('%s/ppdata/subj%02d/%s/brainmask.nii.gz',nsddatadir,subjix,prepname)));
brainmask = flipdim(flipdim(permute(double(a1.img),[2 1 3]),1),2);
clear a1;
  % figure;imagesc(makeimagestack(brainmask));
  % figure;imagesc(makeimagestack(data{1}(:,:,:,10)));
  % drawnow;

% load stimulus
a1 = load(sprintf('%s/stimuli/prf/RETBAR%s.mat',nsddatadir,stimversion));
a2 = load(sprintf('%s/stimuli/prf/RETWEDGERINGMASH%s.mat',nsddatadir,stimversion));
stimulus = {a1.stim a2.stim};   %%stimulus = repmat({a1.stim a2.stim},[1 3]);
clear a1 a2;

% prep data by averaging across identical runs
data{1} = (data{1}+data{3}+data{5})/3;
data{2} = (data{2}+data{4}+data{6})/3;
data = data(1:2);  % now we have only two runs

% truncate data to match stimulus preparation
for p=1:length(data)
  size(data{p})
  size(stimulus{p})
  assert(size(data{p},4) >= size(stimulus{p},3));
  data{p} = data{p}(:,:,:,1:size(stimulus{p},3));
end

% prepare chunks
ix = flatten(find(brainmask));  % total list of voxels to analyze
length(ix)
chunks = chunking(ix,100000);   % break into chunks

% analyze
for cc=chunkstart:length(chunks)
  
  % prepare subset of the data as voxels x TR
  data0 = {};
  for p=1:length(data)
    data0{p} = subscript(squish(data{p},3),{chunks{cc} ':'});
  end

  % call analyzePRF!
  results = analyzePRF(stimulus,data0,tr,struct('seedmode',2,'maxiter',100,'display','off'));

  % save
  save(sprintf('%s/analyzePRF_prfFULL_chunk%04d.mat',glmdir,cc),'-struct','results','-v7.3');
  
  % clean up
  clear results;
  
end
