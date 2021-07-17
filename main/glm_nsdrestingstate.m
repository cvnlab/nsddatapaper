function glm_nsdrestingstate(ppdir,ppdirname,designfile,hrfbasisfile,tr,glmdirname,hrfmanifoldfile)

glmdir = regexprep(ppdir,'ppdata',glmdirname);
mkdirquiet(glmdir);

%%%%% LOAD

% define
stimdur = 3;

% figure out runs
files0 = matchfiles(sprintf('%s/%s/run*.nii',ppdir,ppdirname));
assert(length(files0)==14);
files0 = files0([1 14]);  % first and last are resting-state

% load data
data = {};
for p=1:length(files0)
  a1 = load_untouch_nii(files0{p});
  data{p} = single(a1.img);
end
clear a1;

% load design matrix
stimulus = loadmulti(sprintf('%s/%s',ppdir,designfile),'stimulus');

% massage design matrix (use design matrix of runs 1 and 12 of NSD)
stimulus = stimulus([1 12]);
good = sum(catcell(1,stimulus),1)~=0;
stimulus = cellfun(@(x) x(:,good),stimulus,'UniformOutput',0);
assert(length(stimulus)==length(data));

% truncate data
for p=1:length(data)
  assert(size(data{p},4) >= size(stimulus{p},1));
  data{p} = data{p}(:,:,:,1:size(stimulus{p},1));
end

%%%%% FIT GLM [single-trial, use HRF from NSD, no boots]

% load in the hrfs
load(hrfmanifoldfile,'hrfs');  % TR x different-hrfs

% what was the optimal HRF for each voxel? [this steals the HRF estimate found from the NSD data]
a1 = load(sprintf('%s/GLMdenoise_nsdBASICsingletrialfithrf.mat',glmdir),'HRFindex');  % X x Y x Z

% calc
nx = size(data{1},1);
ny = size(data{1},2);
nz = size(data{1},3);
nh = size(hrfs,2);
nr = length(data);
nb = size(stimulus{1},2);

% initialize
R2 =    zeros(nx*ny*nz,1,'single');    % X x Y x Z with R2 values (all runs)
R2run = zeros(nx*ny*nz,nr,'single');   % X x Y x Z x 2 runs with R2 separated by runs
modelmd = zeros(nx*ny*nz,nb,'single'); % X x Y x Z x trialbetas

% figure out chunking scheme
totalnum = prod(sizefull(data{1},3));
numvoxchunk = 200000;
chunks = chunking(1:size(data{1},3),ceil(size(data{1},3)/ceil(totalnum/numvoxchunk)));

% loop over chunks
for z=1:length(chunks), z

  % loop over possible HRFs
  for hh=1:size(hrfs,2), hh

    % figure out which voxels to process.
    % this will be a vector of indices into the small chunk that we are processing.
    % our goal is to fully process this set of voxels!
    goodix = flatten(find(a1.HRFindex(:,:,chunks{z})==hh));
    
    % extract the data we want to process.
    data0 = cellfun(@(x) subscript(squish(x(:,:,chunks{z},:),3),{goodix ':'}),data,'UniformOutput',0);
    
    % calculate the corresponding indices relative to the full volume
    temp = zeros(size(a1.HRFindex));
    temp(:,:,chunks{z}) = 1;
    relix = subscript(find(temp),goodix);
      
    % do the GLM
    results = GLMestimatemodel(stimulus,data0,stimdur,tr,'assume',hrfs(:,hh),0);
    
    % save results
    R2(relix)        = results.R2;
    R2run(relix,:)   = results.R2run;
    modelmd(relix,:) = results.modelmd{2};
    clear results;

  end

end

% reshape for output
R2      = reshape(R2,     [nx ny nz]);
R2run   = reshape(R2run,  [nx ny nz nr]);
modelmd = reshape(modelmd,[nx ny nz nb]);

% save
save(sprintf('%s/restingGLMdenoise_nsdBASICsingletrialfithrf.mat',glmdir), ...
  'R2','R2run','modelmd','-v7.3');
