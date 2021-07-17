function glm_nsd(ppdir,ppdirname,designfile,hrfbasisfile,tr,glmdirname,hrfmanifoldfile)

glmdir = regexprep(ppdir,'ppdata',glmdirname);
mkdirquiet(glmdir);

%%%%% LOAD

% define
stimdur = 3;

% figure out runs
files0 = matchfiles(sprintf('%s/%s/run*.nii',ppdir,ppdirname));
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
stimulus = loadmulti(sprintf('%s/%s',ppdir,designfile),'stimulus');

% check sanity
assert(length(data)==length(stimulus));

% truncate data
for p=1:length(data)
  assert(size(data{p},4) >= size(stimulus{p},1));
  data{p} = data{p}(:,:,:,1:size(stimulus{p},1));
end

%%%%% FIT GLM TYPE 1 [on-off, assume HRF]

% approach: every stimulus is the same, assume HRF, no boots
stimulus0 = cellfun(@(x) sum(x,2),stimulus,'UniformOutput',0);
results = GLMestimatemodel(stimulus0,data,stimdur,tr,'assume',[],0);

% massage output
results = rmfield(results,{'models' 'modelse' 'signal' 'noise' 'SNR' });
md = results.modelmd;
results.modelmd = md{2};
results.modelmd_hrf = md{1};

% save
save(sprintf('%s/GLMdenoise_nsdBASIConoff.mat',glmdir),'-struct','results','-v7.3');

% clean up
clear results;

%%%%% FIT GLM TYPE 1b [on-off, hrfbasis, full/odd/even]

% load basis
load(hrfbasisfile,'v');

% approach: every stimulus is the same, use canonical hrf basis, no boots
  % sum across all stimuli
stimulus0 = cellfun(@(x) full(sum(x,2)),stimulus,'UniformOutput',0);
  % make 3 separate columns with the convolution output
stimulus0 = cellfun(@(x) [conv2(x,v(:,1)) conv2(x,v(:,2)) conv2(x,v(:,3))],stimulus0,'UniformOutput',0);
  % truncate the excess at the end
stimulus0 = cellfun(@(x,y) x(1:size(y,1),:),stimulus0,stimulus,'UniformOutput',0);
  % fit GLM
results1 = GLMestimatemodel(stimulus0,         data,         stimdur,tr,'assume',1,0);
results2 = GLMestimatemodel(stimulus0(1:2:end),data(1:2:end),stimdur,tr,'assume',1,0);
results3 = GLMestimatemodel(stimulus0(2:2:end),data(2:2:end),stimdur,tr,'assume',1,0);

% massage output
R2 = cat(4,results1.R2,results2.R2,results3.R2);
meanvol = cat(4,results1.meanvol,results2.meanvol,results3.meanvol);
modelmd = cat(5,results1.modelmd{2},results2.modelmd{2},results3.modelmd{2});

% save
save(sprintf('%s/GLMdenoise_nsdBASIConoffhrfbasis.mat',glmdir),'R2','meanvol','modelmd');

% clean up
clear results1 results2 results3 R2 meanvol modelmd;

%%%%% FIT GLM TYPE 2 [single-trial, assume HRF] [EXPENSIVE WRT DISK SPACE]

% approach: single-trial analysis, assume HRF, no boots
results = GLMestimatemodel(stimulus,data,stimdur,tr,'assume',[],0);

% massage output
results = rmfield(results,{'models' 'modelse' 'signal' 'noise' 'SNR' });
md = results.modelmd;
results.modelmd = md{2};
results.modelmd_hrf = md{1};

% save
save(sprintf('%s/GLMdenoise_nsdBASICsingletrial.mat',glmdir),'-struct','results','-v7.3');
  %   R2                177x152x143                   15389088  single              
  %   R2run             177x152x143x12               184669056  single              
  %   hrffitvoxels        0x0                                0  double              
  %   inputs              1x1                            90756  struct              
  %   meanvol           177x152x143                   15389088  single              
  %   modelmd           177x152x143x750            11541816000  single              
  %   modelmd_hrf        52x1                              416  double              

% clean up
clear results;

%%%%% FIT GLM TYPE 3 [single-trial, fit HRF. idea is to choose the best R2 for each voxel.] [EXPENSIVE WRT DISK SPACE AND TIME]

% load in the hrfs
load(hrfmanifoldfile,'hrfs');  % TR x different-hrfs

% approach: single-trial analysis, systematically try different HRFs, no boots

% calc
nx = size(data{1},1)
ny = size(data{1},2)
nz = size(data{1},3)
nh = size(hrfs,2)
nr = length(data)
nb = size(stimulus{1},2)

% initialize
FitHRFR2 =    zeros(nx,ny,nz,nh,'single');     % X x Y x Z x HRFs with R2 values (all runs)
FitHRFR2run = zeros(nx,ny,nz,nr,nh,'single');  % X x Y x Z x runs x HRFs with R2 separated by runs
modelmd =     zeros(nx,ny,nz,nb,'single');     % X x Y x Z x trialbetas

% figure out chunking scheme
totalnum = prod(sizefull(data{1},3));
numvoxchunk = 200000;
chunks = chunking(1:size(data{1},3),ceil(size(data{1},3)/ceil(totalnum/numvoxchunk)))

% loop over chunks
for z=1:length(chunks), z

  % do the fitting and accumulate all the betas
  modelmd0 = zeros(nx,ny,length(chunks{z}),nb,nh,'single');  % X x Y x someZ x trialbetas x HRFs
  for p=1:size(hrfs,2)
    tic; results = GLMestimatemodel(stimulus,cellfun(@(x) x(:,:,chunks{z},:),data,'UniformOutput',0),stimdur,tr,'assume',hrfs(:,p),0); toc;
    FitHRFR2(:,:,chunks{z},p) = results.R2;
    FitHRFR2run(:,:,chunks{z},:,p) = results.R2run;
    modelmd0(:,:,:,:,p) = results.modelmd{2};
    clear results;
  end
  
  % keep only the betas we want
  tic;
  [~,ii] = max(FitHRFR2(:,:,chunks{z},:),[],4);
  modelmd(:,:,chunks{z},:) = matrixindex(modelmd0,repmat(ii,[1 1 1 size(modelmd0,4)]),5);
  clear modelmd0;
  toc;

end

% use R2 to select the best HRF for each voxel
[R2,HRFindex] = max(FitHRFR2,[],4);  % HRFindex is X x Y x Z

% also, use R2 from each run to select best HRF
[~,HRFindexrun] = max(FitHRFR2run,[],5);

% using each voxel's best HRF, what are the corresponding R2run values?
R2run = matrixindex(FitHRFR2run,repmat(HRFindex,[1 1 1 size(FitHRFR2run,4)]),5);  % R2run is X x Y x Z x 12 runs

% save
tic;
save(sprintf('%s/GLMdenoise_nsdBASICsingletrialfithrf.mat',glmdir), ...
  'FitHRFR2','FitHRFR2run','HRFindex','HRFindexrun','R2','R2run','modelmd','-v7.3');
toc;

  % >> a4=load('GLMdenoise_nsdBASICsingletrialfithrf.mat');
  % >> a4
  % 
  % a4 = 
  % 
  %   struct with fields:
  % 
  %        FitHRFR2: [103x80x78x20 single]
  %     FitHRFR2run: [5-D single]
  %        HRFindex: [103x80x78 double]
  %     HRFindexrun: [103x80x78x12 double]
  %              R2: [103x80x78 single]
  %           R2run: [103x80x78x12 single]
  %         modelmd: [103x80x78x750 single]

  % >> test=load('GLMdenoise_nsdBASIConoff.mat');
  % >> test
  % 
  % test = 
  % 
  %   struct with fields:
  % 
  %               R2: [103x80x78 single]
  %            R2run: [103x80x78x12 single]
  %     hrffitvoxels: []
  %           inputs: [1x1 struct]
  %          meanvol: [103x80x78 single]
  %          modelmd: [103x80x78 single]
  %      modelmd_hrf: [39x1 double]

%%%%% FIT GLM TYPE 3B [single-trial, use optimal HRF for each voxel, GLMdenoise, do ridge regression, no boots] [EXPENSIVE WRT DISK SPACE AND TIME]

% other relevant players: fracridge,olsmatrix2,GLMestimatemodel

% Notes:
% - two different resolutions independently processed
% - there is no random aspect of this analysis
% - we only cross-validate the first training instance
% - we assume 1/20 granularity, but allocate a few more at the top end.
% - the effect, ultimately, looks like some amount of temporal smoothing.
%   but this is only because of the positive correlation between consecutive trials.
% - we assume no interaction between PC and frac
% - it appears that frac is finesse and PC is like a sledgehammer
% - we try up to 10 PCs
% - key question is whether this has any chance to work:
%   noise regressors into a GLM that is high freedom already
%   (apparently yes)
% - what is the bottom line improvement? log squared beta error is hard to look at.
%   the key is the beta stability.
% - note that we first tune the PCs and then we bake that in and then tune frac
% - it appears we fully resolve (estimate) the nuisance (poly + noise regressors),
%   which seems to suggest their weights are with respect to the OLS beta solutions.

% what was the optimal HRF for each voxel?
a1 = load(sprintf('%s/GLMdenoise_nsdBASICsingletrialfithrf.mat',glmdir),'HRFindex');  % X x Y x Z

% load in the hrfs
load(hrfmanifoldfile,'hrfs');  % TR x different-hrfs

% calc
nx = size(data{1},1);
ny = size(data{1},2);
nz = size(data{1},3);
nh = size(hrfs,2);
nr = length(data);
nb = size(stimulus{1},2);

% figure out chunking scheme
totalnum = prod(sizefull(data{1},3));
numvoxchunk = 200000;
chunks = chunking(1:size(data{1},3),ceil(size(data{1},3)/ceil(totalnum/numvoxchunk)));

% define fractions
fracs = fliplr([.05:.05:.90 .91:.01:1]);  % no bias -> huge bias

% load in the stim order for this session
stimorder = loadmulti(sprintf('%s/%s',ppdir,designfile),'stimorder');  % 1 x 750 with image ids

% define xval scheme (6-fold)
xvals = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12]};

% calculate trial information
validcolumns = {};   % 1 x 12 cell, each element is the vector of trial indices (relative to 750) associated with a given run
stimix = {};         % 1 x 12 cell, each element is the vector of actual image ids associated with a given run
for p=1:12
  if mod(p,2)==1
    validcolumns{p} = (ceil(p/2)-1)*125 + (1:63);
  else
    validcolumns{p} = (ceil(p/2)-1)*125 + 63 + (1:62);
  end
  stimix{p} = stimorder(validcolumns{p});
end

% define
numpcstotry = 10;

%%%%% whew, proceed!

% figure out the noise pool
a8 = load(sprintf('%s/GLMdenoise_nsdBASIConoff.mat',glmdir),'meanvol','R2');  % load in results from a previous call
r2thresh = findtailthreshold(a8.R2(:),0);% figure out a threshold. this is pretty aggressive throwing away signal voxels.
thresh = prctile(a8.meanvol(:),99)*0.1;  % threshold for non-brain voxels. this allows a lot in.
bright = a8.meanvol > thresh;            % logical indicating voxels that are bright (brain voxels)
badxval = a8.R2 < r2thresh;              % logical indicating voxels with poor R2 [NOTE: this is imperfect]
noisepool = bright & badxval;            % logical indicating voxels that satisfy all criteria
%   % remove soon:
%   outputtempdir = sprintf('~/Dropbox/KKTEMP/%s',stripfile(ppdir,1));
%   mkdirquiet(outputtempdir);
%   imwrite(uint8(255*makeimagestack(a8.meanvol,1)),gray(256),[outputtempdir '/meanvol.png']);
%   imwrite(uint8(255*makeimagestack(bright,[0 1])),gray(256),[outputtempdir '/bright.png']);
%   imwrite(uint8(255*makeimagestack(a8.R2,[0 10])),hot(256), [outputtempdir '/r2.png']);
%   imwrite(uint8(255*makeimagestack(badxval,[0 1])),gray(256), [outputtempdir '/badxval.png']);
%   imwrite(uint8(255*makeimagestack(noisepool,[0 1])),gray(256), [outputtempdir '/noisepool.png']);

% define
dimdata = 3;

% determine noise regressors
pcregressors = {};
for p=1:length(data)

  % extract the time-series data for the noise pool
  temp = subscript(squish(data{p},dimdata),{find(noisepool) ':'})';  % time x voxels

  % project out polynomials from the data
  temp = projectionmatrix(constructpolynomialmatrix(size(temp,1),0:3)) * temp;

  % unit-length normalize each time-series (ignoring any time-series that are all 0)
  [temp,len] = unitlengthfast(temp,1);
  temp = temp(:,len~=0);

  % perform SVD and select the top PCs
  [u,s,v] = svds(double(temp*temp'),numpcstotry);
  u = bsxfun(@rdivide,u,std(u,[],1));  % scale so that std is 1
  pcregressors{p} = cast(u,'single');

end
clear temp len u s v;

%%%%% okay, start with GLMdenoise. the goal is to figure out the number of PCs.

% initialize
totalbadness = zeros(nx*ny*nz,1+numpcstotry,'single');    % X * Y * Z x 1+npc  [squared beta error for different numbers of PCs]

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

    % perform GLMdenoise
    clear results;
    for ll=0:numpcstotry

      % define options
      opt = struct('lambda',0,'wantpercentbold',0);
      opt.extraregressors = cell(1,length(data0));
      if ll>0
        for rr=1:length(data0)
          opt.extraregressors{rr} = cat(2,opt.extraregressors{rr},pcregressors{rr}(:,1:ll));
        end
      end
      
      % do the GLM
      [results(ll+1),cache] = GLMestimatemodel(stimulus,data0, ...
                                  stimdur,tr,'assume',hrfs(:,hh),0,opt);

      % save some memory
      results(ll+1).models = [];
      results(ll+1).modelse = [];
    
    end

    % compute the cross-validation performance values
    totalbadness(relix,:) = calcbadness(xvals,validcolumns,stimix,results);  % voxels x regularization levels

  end
  
end

% pick number of PCs
r2finalthresh = 5;  % R2 threshold to use for figuring out number of PCs
xvaltrend = -median(totalbadness(a8.R2(:)>r2finalthresh,:),1);  % NOTE: sign flip so that high is good
  % note that r2 threshold here is imperfect, but useful

% choose number of PCs
chosen = 0;  % this is the fall-back
curve = xvaltrend - xvaltrend(1);  % this is the performance curve that starts at 0 (corresponding to 0 PCs)
mx = max(curve);                   % store the maximum of the curve
best = -Inf;                       % initialize (this will hold the best performance observed thus far)
for p=0:numpcstotry

  % if better than best so far
  if curve(1+p) > best

    % record this number of PCs as the best
    chosen = p;
    best = curve(1+p);
  
    % if we are within opt.pcstop of the max, then we stop.
    if best*1.05 >= mx
      break;
    end
  
  end

end

% record the number of PCs
pcnum = chosen;

% DEPRECATED
%
% % visualize
% figureprep; hold on;
% rvals = [1 3 5 10 20 30];
% cmap0 = jet(length(rvals));
% for pp=1:length(rvals)
%   temp = totalbadness(a8.R2(:)>rvals(pp),:);
%   plot(0:numpcstotry,calczscore(median(temp,1)),'-','Color',cmap0(pp,:));
% end
% straightline(pcnum,'v','k-');
% xlabel('number of pcs');
% ylabel('median badness, z-scored');
% figurewrite('checkbadness',[],[],outputtempdir);
% 
% % visualize  [PERHAPS GO BACK TO LINEAR; USE SCATTERSPARSE?]
% rvals = [1 5 20];
% colors = {'r' 'g' 'b'};
% for p=1:numpcstotry
%   figureprep([100 100 900 900]);
%   for cc=1:length(rvals)
%     temp = totalbadness(a8.R2(:)>rvals(cc),:);
%     scatter(log(temp(:,1)),log(temp(:,1+p)),[colors{cc} '.']);
%   end
%   axissquarify;
%   %ax = axis;
%   %axis([0 ax(2) 0 ax(2)]);
%   xlabel('log error for no pcs');
%   ylabel(sprintf('log error for %d pcs',p));
%   figurewrite(sprintf('scatter%02d',p),[],[],outputtempdir);
% end

% deal with dimensions
totalbadness = reshape(totalbadness,nx,ny,nz,[]);

%%%%% okay, proceed to ridge regression

% initialize
modelmd =     zeros(nx*ny*nz,nb,'single');     % X * Y * Z x trialbetas  [the final beta estimates]
R2 =          zeros(nx,ny,nz,'single');        % X x Y x Z               [the R2 for the specific optimal frac]
R2run =       zeros(nx*ny*nz,nr,'single');     % X * Y * Z x runs        [the R2 separated by runs for the optimal frac]
FRACvalue   = zeros(nx,ny,nz,'single');        % X x Y x Z               [best fraction]
scaleoffset = zeros(nx*ny*nz,2,'single');      % X * Y * Z x 2           [scale and offset]

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

    % process each frac
    clear results;
    for ll=1:length(fracs)

      % define options
      opt = struct('frac',fracs(ll),'wantpercentbold',0);
      opt.extraregressors = cell(1,length(data0));
      if pcnum > 0
        for rr=1:length(data0)
          opt.extraregressors{rr} = cat(2,opt.extraregressors{rr},pcregressors{rr}(:,1:pcnum));
        end
      end

      % fit the entire dataset using the specific frac
      [results(ll),cache] = GLMestimatemodel(stimulus,data0, ...
                                  stimdur,tr,'assume',hrfs(:,hh),0,opt);
        % results = 
        % 
        %   struct with fields:
        % 
        %           models: {[31x1 double]  [3902x750 single]}
        %          modelmd: {[31x1 double]  [3902x750 single]}
        %          modelse: {[31x1 double]  [3902x750 single]}
        %               R2: [3902x1 single]
        %            R2run: [3902x12 single]
        %           signal: [3902x1 single]
        %            noise: [3902x1 single]
        %              SNR: [3902x1 single]
        %     hrffitvoxels: []
        %          meanvol: [3902x1 single]
        %           inputs: [1x1 struct]
      
      % save some memory
      results(ll).models = [];
      results(ll).modelse = [];
    
    end
    
    % compute the cross-validation performance values
    badness = calcbadness(xvals,validcolumns,stimix,results);
    
    % pick best frac
    [~,FRACindex0] = min(badness,[],2);  % FRACindex0 is V x 1 with the index of the best frac
    
    % prepare output
    FRACvalue(relix) = fracs(FRACindex0);
    numbetas = size(results(1).modelmd{2},2);
    for ll=1:length(fracs)
      ii = find(FRACindex0==ll);
      
      % scale and offset to match the unregularized result
      for vv=1:length(ii)
        X = [results(ll).modelmd{2}(ii(vv),:); ones(1,numbetas)]';
        h = olsmatrix(X)*results(1).modelmd{2}(ii(vv),:)';
        if h(1) < 0
          h = [1 0]';
        end
        scaleoffset(relix(ii(vv)),:) = h;
        modelmd(relix(ii(vv)),:) = X*h;
      end
      
      R2(relix(ii))        = results(ll).R2(ii);
      R2run(relix(ii),:)   = results(ll).R2run(ii,:);
    end

  end

end

% deal with dimensions
modelmd =         reshape(modelmd,[nx ny nz nb]);
R2run =           reshape(R2run,  [nx ny nz nr]);
scaleoffset = reshape(scaleoffset,[nx ny nz 2]);

% load in mean volume from a previous call
meanvol = loadmulti(sprintf('%s/GLMdenoise_nsdBASICsingletrial.mat',glmdir),'meanvol');

% deal with percent BOLD change
modelmd = bsxfun(@rdivide,modelmd,abs(meanvol)) * 100;

% save
tic;
save(sprintf('%s/GLMdenoise_nsdBASICsingletrialfithrfGLMdenoiseFracridge.mat',glmdir), ...
  'totalbadness','pcnum','noisepool','pcregressors', ...
  'modelmd','R2','R2run','FRACvalue','scaleoffset','-v7.3');
toc;

%%%%%

function badness = calcbadness(xvals,validcolumns,stimix,results)

% initialize
badness = zeros(size(results(1).modelmd{2},1),length(results));  % voxels x hyperparameters with the sum of the squared error from cross-validation

% calc
alltheruns = catcell(2,xvals);

% do cross-validation
for xx=1:length(xvals)

  % calc
  testix = xvals{xx};              % which runs are testing, e.g. [3 4]
  trainix = setdiff(alltheruns,testix);  % which runs are training, e.g. [1 2 5 6 7 8 9 10 11 12]
  
  % calc
  testcols = catcell(2,validcolumns(testix));    % vector of trial indices in the testing data (out of 750)
  traincols = catcell(2,validcolumns(trainix));  % vector of trial indices in the training data (out of 750)
  testids = catcell(2,stimix(testix));           % vector of image-ids in the testing data
  trainids = catcell(2,stimix(trainix));         % vector of image-ids in the training data
  
  % for each trial, figure out what to do
  topullTEST = [];    % vector of beta indices to use from the testing data (these are indices into the 750)
  topullTRAIN = [];   % corresponding vector of beta indices to use from the training data (these are indices into the 750)
  for ttt=1:length(testids)
    haveix = find(trainids==testids(ttt));
    if ~isempty(haveix)
      topullTEST = [topullTEST testcols(ttt)];
      topullTRAIN = [topullTRAIN traincols(haveix(1))];  % ASSUME THE FIRST ONE ONLY
    end
  end
  
  % topullTEST tells us which trials (isolated within the testing runs) to pull betas for
  % topullTRAIN tells us the corresponding trials (isolated within the training runs) to pull betas for
  
  % calculate cross-validation performance
  for ll=1:length(results)
    badness(:,ll) = badness(:,ll) + sum((results(ll).modelmd{2}(:,topullTRAIN) - results(1).modelmd{2}(:,topullTEST)).^2,2);
  end

end
