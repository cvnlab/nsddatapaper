function nsdfaceprf(subjix,p,q,rr,hhs)

% function nsdfaceprf(subjix,p,q,rr,hhs)
%
% For subject <subjix>, model <p>, version <q>, resampling <rr>, and <hhs> vector,
%   try to find chunks to do, do them, and keep trying.
%
% Call analyzePRF on appropriate prepared data and save results of different chunks.

%% %%%%% SETUP

fprintf('*** starting nsdfaceprf on subjix=%d,p=%d,q=%d,rr=%d,hhs=%s\n',subjix,p,q,rr,mat2str(hhs));

% pre-processed stimulus name
switch p
case {1 2 5 6 7 8 12}  % Face and Word and Body and Foreground
  if ismember(p,[7 8])
    ptouse = 6;
  elseif ismember(p,[12])
    ptouse = 1;
  else
    ptouse = p;
  end
  stimname = sprintf('stimprf_p%d_q%d',ptouse,q);
case {3 4 9}      % contrastgrid
  stimname = 'stimcontrastgrid';
case {10 11}         % contrastgridNEW
  stimname = 'stimcontrastgridNEW';
case {13 14}
  stimname = 'stimsalience';
otherwise
  die;
end
switch p
case {1 2 5 9 12}
  maxpolydeg = 0;
case {3 4 6 7 8 10 11 13 14}
  maxpolydeg = NaN;
otherwise
  die;
end
switch p
case {9}
  modelmode = 3;  % do not optimize exponent
case {1 2 3 4 5 6 7 8 10 11 12 13 14}
  modelmode = 2;
otherwise
  die;
end

% load basic information (subjix is defined)
LOADDATA1;

%% %%%%% START MASTER LOOP

while 1

  % are there any jobs that need doing?
  foundone = 0;
  for hh=hhs
    for cc=1:ceil(datasizes(hh,1)/numinchunk)
      datafile = sprintf('~/ext/figurefiles/nsd/TOCONSOLIDATE/nsdfaceprf_subj%02d_p%d_q%d_hh%d_rr%d_cc%05d.mat',subjix,p,q,hh,rr,cc);
      if ~exist(datafile,'file')
        a = 1;
        save(datafile,'a');  % "lock" it
        foundone = 1;
        break;
      end
    end
    if foundone==1
      break;
    end
  end

  % if nothing to be done, we are done!
  if foundone==0
    return;
  end
  tic;
  unix(sprintf('echo "START: computer is %s, parpool size is %d, datafile is %s" >> ~/nsdfaceprflog',gethostname,getparpoolsize,datafile));

  % load in pre-processed stimulus (if necessary)
  if ~exist('stim','var')
    stim = loadmulti(sprintf('~/ext/figurefiles/nsd/%s_versionA.mat',stimname),'stim'); % 50 x 50 x images, range 0-1

    % flip?
    if ismember(p,[8])
      stim = 1-stim;
    end
  end

  % since hh and cc are defined, let's load data (expensive)
  LOADDATA2;

  % prepare stimulus and data
  if q==1
    stim0 = stim(:,:,expdesign0.subjectim(subjix,ordU));  % get the full set (50 x 50 x NIMAGES (0/1))
    stim0 = stim0(:,:,trainixsALL{q,rr});                      % get the trainixs ones (50 x 50 x THISROUND)
  else
    stim0 = stim(:,:,valix2);                              % which of the shared1000 do we have (as validation images)? (50 x 50 x VALIMAGES)
    stim0 = stim0(:,:,ismember(valix,trainixsALL{q,rr}));       % get the trainixs ones (50 x 50 x THISROUND)
  end
  
  % sign change
  if ismember(p,[4 7 11 12 14])
    sgn = -1;
  else
    sgn = 1;
  end

  %%%% BASELINE ADJUSTMENT BEGIN

% actually, we don't like this.  
%   if ismember(p,[3 4])
%     
%     X = squish(double(stim0),2)';                        % images x 50*50
%     y = sgn*double(somedata(:,trainixsALL{q,rr}))';      % images x voxels
%     [ixtest,~,ixtrain] = picksubset(1:size(X,1),[5 1]);  % ixtest is indices in 20%; ixtrain is indices in the remaining 80%
%     [coef,alphas,offset] = fracridge(X(ixtrain,:),[0:.01:.04 .05:.05:1],y(ixtrain,:),[],[],2);
%     errs = [];
%     for ff=1:size(coef,2)
%       errs(ff,:) = sum((y(ixtest,:) - (X(ixtest,:)*permute(coef(:,ff,:),[1 3 2]) + offset(ff,:))).^2,1);  % 1 x voxels
%     end
%     [~,mnix] = min(errs,[],1);                 % 1 x voxels with the index
%     offsetFINAL = matrixindex(offset,mnix,1);  % 1 x voxels with the actual offset
%     clear X y coef alphas offset errs;
%     
%   else
%    offsetFINAL = 0;
%   end
  
  %%%% BASELINE ADJUSTMENT END

  % run analyzePRF (double format; full mode; fake HRF; baseline constant term)
%   if ismember(p,[3 4])
%     typicalgain = 2;     % for contrast model, we want to uniformly set this to 2
%   else
    typicalgain = NaN;   % for all other cases, we want to do the empirical strategy
%   end
  results = analyzePRF(double(stim0), ...
                       sgn*double(somedata(:,trainixsALL{q,rr})), ...    % - offsetFINAL'
                       1,struct('seedmode',2,'modelmode',modelmode,'hrf',1,'maxpolydeg',maxpolydeg,'maxiter',100,'display','off', ...
                                'algorithm','trust-region-reflective','typicalgain',typicalgain));

  % NOTE: WE HAVE TO BE CAREFUL ABOUT DOUBLE vs SINGLE!!

  % This is what a file looks like.
  % >> a2=load('nsdfaceprf_subj01_p3_q1_hh1_rr1_cc00001.mat')
  % 
  % a2 = 
  % 
  %   struct with fields:
  % 
  %           R2: [10000x1 double]
  %          ang: [10000x1 double]
  %          ecc: [10000x1 double]
  %         expt: [10000x1 double]
  %         gain: [10000x1 double]
  %      meanvol: [10000x1 double]
  %     noisereg: []
  %     numiters: {10000x1 cell}
  %      options: [1x1 struct]
  %       params: [1x5x10000 double]
  %     resnorms: {10000x1 cell}
  %       rfsize: [10000x1 double]

  % some model setup
  nn = size(stim0,1);  % NOTE: this can be either 50 or 51
  [~,xx,yy] = makegaussian2d(nn,2,2,2,2);
  modelfun = @(pp,dd) posrect(pp(4)) * (dd*vflatten(makegaussian2d(nn,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))) .^ posrect(pp(5));
  
  % notes:
  % mx = mean of data
  % my = mean of model fit (gain*model)
  % model fit + baseline = data
  % baseline = data - model fit
  % baseline term is just (mean of data) - (mean of model fit)

  % figure out the gain and baseline
  stim0 = squish(double(stim0),2)';  % images x pixels*pixels
  results.baseline = zeros(size(results.params,3),1);  % initialize
  if ~isnan(maxpolydeg)  % in this case, we have to posthoc compute it (since analyzePRF doesn't give it to us)
    fprintf('starting gain/baseline');
    temp = results.baseline;
    temp2 = trainixsALL{q,rr};
    somedata = somedata;
    parfor vx=1:size(results.params,3)
      statusdots(vx,size(results.params,3));
      thefit0 = modelfun(results.params(1,:,vx),stim0);  % images x 1
      temp(vx) = mean(sgn*double(somedata(vx,temp2))) - mean(thefit0);
    end
    results.baseline = temp;
    fprintf('done.\n');
  end
    
  % now, real model setup (with explicit baseline)
  modelfunB = @(pp,dd,bb) bb + posrect(pp(4)) * (dd*vflatten(makegaussian2d(nn,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))) .^ posrect(pp(5));

  % massage outputs
  results.ecc =    results.ecc*pxtodeg;
  results.rfsize = results.rfsize*pxtodeg;
  results.ang(results.ecc==0) = NaN;

  % compute predictions for image-level prediction accuracy
  if q==1
    temp = double(stim(:,:,expdesign0.subjectim(subjix,ordU)));
    temp = squish(temp(:,:,valix),2)';
  else
    temp = squish(double(stim(:,:,valix2)),2)';
  end
  themodelfit = zeros(size(results.params,3),size(temp,1));
  fprintf('image-level prediction');
  parfor vx=1:size(results.params,3)
    statusdots(vx,size(results.params,3));
    themodelfit(vx,:) = modelfunB(results.params(1,:,vx),temp,results.baseline(vx));
  end
  thedata = sgn*double(somedata(:,valix));  % V x VALIMAGES    (NOTE: no -offsetFINAL needed here)
  predR2 = calccod(themodelfit,thedata,2,0,0);     % V x 1
  predr = calccorrelation(themodelfit,thedata,2);  % V x 1
  fprintf('done.\n');
  
  % record
  results.predR2 = predR2;
  results.predr  = predr;
  
  % SAVE
  save(datafile,'-struct','results','-v7.3');

  % report
  unix(sprintf('echo "FINISH: computer is %s, time elapsed is %.1f min, parpool size is %d, datafile is %s" >> ~/nsdfaceprflog',gethostname,toc/60,getparpoolsize,datafile));

end
