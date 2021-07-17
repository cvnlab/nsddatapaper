% This makes some nice prepared versions of the data:
%  ~/ext/figurefiles/nsd/datab3nativesurface_subj[01-08].mat
%  ~/ext/figurefiles/nsd/datab3nativesurface_subj[01-08]_[lh,rh,subcortex,func1pt8mm,cerebellum]_betas.hdf5
%
% In the .mat file, we have:
%   <ord> is 1 x NTRIALS with 10k indices
%   <ordU> is 1 x NIMAGES with unique 10k indices of images shown to this subject
%   <allixs> is 3 x NIMAGES indicating how to pull trials from the betas
%   <numtrialsperim> is 1 x NIMAGES with the actual number of trials contributing to each (1, 2, 3)
%   <alldata> is 1 x [lh,rh,subcortex,func1pt8mm,cerebellum] with VERTICES/VOXELS x NIMAGES with the responses
%   <d1> is A:B (indices into 1st dimension of func1mm)
%   <d2> is C:D (indices into 2nd dimension of func1mm)
%   <d3> is E:F (indices into 3rd dimension of func1mm)
%   <dsz> is the subcortex matrix size, like [20 42 34]
%   <bmii> is the logical X x Y x Z brainmask for the 1.8-mm volume preparation
%   <cd1>,<cd2>,<cd3>,<cii>,<cdsz> is for cerebellum
%   <datasizes> is 5 x 2 (hh x nimages)
%   <adjparams> has [meanmu meansigma b s2 reversionscale reversionoffset]

%% %%%%%%%%%%%%%%%%%%%%%%% PREPARE DATA (load all data (b3, nativesurface, 1mm for subcortex/cerebellum))

% define
bdir = 'betas_fithrf_GLMdenoise_RR';

% load some basic stuff
a1 = load('~/nsd/nsddata/experiments/nsd/nsd_expdesign.mat');
nsess = [40 40 32 30 40 32 40 30];
hemis = {'lh' 'rh'};
datanames = {'lh' 'rh' 'subcortex' 'func1pt8mm' 'cerebellum'};

% define
subcortexpad = [5 5 10 5 5 5];

% loop
for subjix=1:8

  % subcortex stuff
  g1 = load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj%02d/func1mm/roi/thalamus.nii.gz',subjix));
  [d1,d2,d3,ii] = computebrickandindices(double(g1.img)>0);
  d1 = d1(1)-subcortexpad(1):d1(end)+subcortexpad(2);
  d2 = d2(1)-subcortexpad(3):d2(end)+subcortexpad(4);
  d3 = d3(1)-subcortexpad(5):d3(end)+subcortexpad(6);
  dsz = [length(d1) length(d2) length(d3)];
  
  % func1pt8mm stuff
  bm1 = load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj%02d/func1pt8mm/brainmask.nii.gz',subjix));
  bmii = logical(double(bm1.img)>0);  % X x Y x Z binary mask
  
  % cerebellum stuff
  c1 = load_untouch_nii('~/Dropbox/KKTEMP/cerebellum/cerebellum_retmaps_1mm_binary_dilM.nii.gz');  % LPI format
  c1new = nsd_mapdata(subjix,'MNI','func1pt0',flipdim(c1.img,1),'wta');
  [cd1,cd2,cd3,cii] = computebrickandindices(double(c1new)==1);
  cdsz = [length(cd1) length(cd2) length(cd3)];

  % experimental design stuff
  ord = a1.masterordering(1:750*nsess(subjix));  % 1 x NTRIALS with 10k indices
  ordU = unique(ord);                            % 1 x IMAGES with all unique 10k images shown to this subject
  allixs = [];                                   % 3 x UNIQUEIMAGES (this sets the order)
  badtrialix = length(ord)+1;
  for qq=1:length(ordU)
    ix = find(ord==ordU(qq));                     % indices of trial numbers that ordU(qq) was shown on
    ix = vflatten(ix);                            % make a column vector
    ix = [ix; badtrialix*ones(3-length(ix),1)];   % include some dummy entries if necessary
    allixs(:,end+1) = ix;                         % record
  end
  numtrialsperim = sum(allixs~=badtrialix,1);     % 1 x UNIQUEIMAGES with number of trials (1/2/3)
  
  % now, allixs can be used to pull trials from the betas
  
  % load betas
  betadir0 =  sprintf('~/nsd/nsddata_betas/ppdata/subj%02d/nativesurface/%s',subjix,bdir);
  betadir0b = sprintf('~/nsd/nsddata_betas/ppdata/subj%02d/func1mm/%s',      subjix,bdir);
  betadir0c = sprintf('~/nsd/nsddata_betas/ppdata/subj%02d/func1pt8mm/%s',   subjix,bdir);
  alldata =   cell(1,length(datanames));
  adjparams = cell(1,length(datanames));
  for hh=1:length(datanames)
    alldata{hh} =   single([]);
    adjparams{hh} = single([]);
    
    % load
    for nn=nsess(subjix):-1:1, nn
      if hh <= 2
        file0 = sprintf('%s/%s.betas_session%02d.hdf5',betadir0,hemis{hh},nn);
        temp = h5read(file0,'/betas',[1 1],[Inf Inf]);  % 227021 x 750, int16
      elseif hh == 3
        file0 = sprintf('%s/betas_session%02d.hdf5',betadir0b,nn);
        temp = squish(h5read(file0,'/betas',[d1(1) d2(1) d3(1) 1],[range(d1)+1 range(d2)+1 range(d3)+1 Inf]),3);  % X*Y*Z x 750, int16
        size(temp)
      elseif hh == 4
        file0 = sprintf('%s/betas_session%02d.hdf5',betadir0c,nn);
        temp = squish(h5read(file0,'/betas',[1 1 1 1],[Inf Inf Inf Inf]),3);  % X*Y*Z x 750, int16
        temp = temp(find(bmii),:);  % voxels x 750, int16
        size(temp)
      elseif hh == 5
        file0 = sprintf('%s/betas_session%02d.hdf5',betadir0b,nn);
        temp = squish(h5read(file0,'/betas',[cd1(1) cd2(1) cd3(1) 1],[range(cd1)+1 range(cd2)+1 range(cd3)+1 Inf]),3);  % X*Y*Z x 750, int16
        temp = temp(cii,:);  % voxels x 750, int16
        size(temp)
      end
      temp = single(temp)/300;  % PSC units; some may be all 0 (invalid voxels)
      alldata{hh}(:,:,nn) = temp;
      clear temp;
    end

    % proceed to normalization

    % calculate the additive adjustment factor
    mus    = mean(alldata{hh},2);   % V x 1 x session
    sigmas = std(alldata{hh},[],2); % V x 1 x session
    b = nanmean(zerodiv(mus,sigmas,NaN),3);  % V x 1
    if ~all(isfinite(b))
      warning('bad! 1');  % TURNS OUT subj06 has a few completely bad voxels on the outskirts.
    end

    % zscore each session and concatenate (NOTE: invalid voxels will be NaN for a given scan session)
    alldata{hh} = reshape(calczscore(alldata{hh},2),size(alldata{hh},1),[]);  % V x 750*nsess

    % add adjustment factor
    alldata{hh} = bsxfun(@plus,alldata{hh},b);
    
    % add a fake column
    alldata{hh}(:,end+1) = NaN;
    
    % trial-average (after this step, there could be NaNs for some images, e.g. if all reps of an image occurred in a bad session)
    alldata{hh} = squish(nanmean(reshape(alldata{hh}(:,flatten(allixs)),size(alldata{hh},1),3,[]),2),2);  % V x UNIQUEIMAGES
    
    % to make life easy, set NaNs to 0
    alldata{hh}(isnan(alldata{hh})) = 0;
    
    % record adjustment factors
    adjparams{hh}(:,[1 2 3]) = [mean(mus,3) mean(sigmas,3) b];

  end
  
  % "re-zscore"
  for hh=1:length(alldata)
    s2 = std(alldata{hh},[],2);
    if any(s2==0)
      warning('bad! 2');  % TURNS OUT subj06 has a few completely bad voxels on the outskirts.
    end
    alldata{hh} = zerodiv(alldata{hh},s2,0,0);
    adjparams{hh}(:,4) = s2;
  end
  
  % now, 0 is still meaningful, and the std dev of each vertex is 1 (but note below).
  % since this last step is just scaling the data, the noise ceiling idea is still valid.
  % HOWEVER, there are a few bogus cases where the entire dataset is all 0.
  
  % tests:
  %roi1 = cvnloadmgz('~/nsd/nsddata/freesurfer/subj01/label/rh.floc-faces.mgz');
  %temp = alldata{2}(roi1==2,:);
  %figure; imagesc(temp);
  
  % compute datasizes
  datasizes = [];
  for hh=1:length(alldata)
    datasizes(hh,:) = sizefull(alldata{hh},2);
  end
  
  % compute "reversion" scale and offset.
  % in a sense, we want (vals*s2 - b) * mean(sigmas,3) + mean(mus,3).
  for hh=1:length(alldata)
    meanmu =    adjparams{hh}(:,1);
    meansigma = adjparams{hh}(:,2);
    b =         adjparams{hh}(:,3);
    s2 =        adjparams{hh}(:,4);
    b(~isfinite(b)) = 0;  % when b is NaN, we actually didn't change anything
    s2(s2==0) = 1;        % when s2 is 0, we actually didn't change anything
    adjparams{hh}(:,5) = s2 .* meansigma;                % scale factor for reversion (V x 1)
    adjparams{hh}(:,6) = -b .* meansigma + meanmu;       % offset for reversion (V x 1)
  end

  % save preprocessed data
  save(sprintf('~/ext/figurefiles/nsd/datab3nativesurface_subj%02d.mat',subjix), ...
    'ord','ordU','allixs','numtrialsperim','d1','d2','d3','dsz','datasizes','bmii','cd1','cd2','cd3','cii','cdsz','adjparams');
  for hh=1:length(alldata)
    outputfile = sprintf('/home/surly-raid2/kendrick-data/figurefiles/nsd/datab3nativesurface_subj%02d_%s_betas.hdf5',subjix,datanames{hh});
    delete(outputfile);
    h5create(outputfile,'/betas',size(alldata{hh}),'Datatype','single','ChunkSize',[1 size(alldata{hh},2)]);
    h5write(outputfile,'/betas',alldata{hh});
  end

end

% ONE-OFF HACK IN
%
% datanames = {'lh' 'rh' 'subcortex'};
% for subjix=1:8
%   datasizes = [];
%   for hh=1:length(datanames)
%     t1 = h5info(sprintf('~/ext/figurefiles/nsd/datab3nativesurface_subj%02d_%s_betas.hdf5',subjix,datanames{hh}),'/betas');
%     datasizes(hh,:) = t1.Dataspace.Size;
%   end
%   datasizes
%   save(sprintf('~/ext/figurefiles/nsd/datab3nativesurface_subj%02d.mat',subjix),'datasizes','-append');
% end
