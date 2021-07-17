function nsdfaceprffinish(subjix,p,q,rr,hhs,mtodo,wantsavedata)

% function nsdfaceprffinish(subjix,p,q,rr,hhs,mtodo,wantsavedata)
%
% For subject <subjix>, model <p>, version <q>, resampling <rr>, and <hhs> vector...

% inputs
if ~exist('mtodo','var') || isempty(mtodo)
  mtodo = 3;
end
if ~exist('wantsavedata','var') || isempty(wantsavedata)
  wantsavedata = 1;
end

% define
types = {'Faces' 'Words' 'Contrastgrid'};
labeltypes = {'Automated' 'Semi-automated'};
modelnames = {'category' 'prf' 'fullprf'};
if p==1
  if q==1
    finalname = 'faceauto';
  else
    finalname = 'facesemi';
  end
elseif p==2
  if q==1
    finalname = 'wordauto';
  else
    finalname = 'wordsemi';
  end
elseif p==3
  assert(q==1);
  finalname = 'contrast';
elseif p==4
  assert(q==1);
  finalname = 'contrastneg';
elseif p==5
  assert(q==1);
  finalname = 'bodyauto';
elseif p==6
  assert(q==1);
  finalname = 'foregroundauto';
elseif p==7
  assert(q==1);
  finalname = 'foregroundautoneg';
elseif p==8
  assert(q==1);
  finalname = 'backgroundauto';
elseif p==9
  assert(q==1);
  finalname = 'contrastBASELINE';
elseif p==10
  assert(q==1);
  finalname = 'contrastNEW';
elseif p==11
  assert(q==1);
  finalname = 'contrastNEWneg';
elseif p==12
  assert(q==1);
  finalname = 'faceautoneg';
elseif p==13
  assert(q==1);
  finalname = 'salience';
end

% special prefix
if mtodo==1
  finalname = ['category' finalname];
end

% load
LOADDATA1;

% ???
% category model applies to p=1,2 (and q=1,2)
% prf is the quick model and applies to p=1,2 (and q=1,2)
% fullprfcg is the fully-fit contrast grid pRF model and applies to ...

% calc
finalfile0 = sprintf('~/ext/figurefiles/nsd/%sresults_subj%02d.mat',modelnames{mtodo},subjix);

% init
if exist(finalfile0,'file')
  load(finalfile0,'allresults');
else
  allresults = {};
end

% do it
if mtodo==3

  for hh=hhs

    % find all files and load them in
    allfiles0 = matchfiles(sprintf('~/ext/figurefiles/nsd/TOCONSOLIDATE/nsdfaceprf_subj%02d_p%d_q%d_hh%d_rr%d_cc*.mat',subjix,p,q,hh,rr));
    varstoload = {'R2' 'ang' 'ecc' 'expt' 'gain' 'rfsize' 'params' 'baseline' 'predR2' 'predr'};
    dimtoagg = [1 1 1 1 1 1 3 1 1 1];
    clear results;
    for vv=1:length(varstoload)
      results.(varstoload{vv}) = loadmulti(allfiles0,varstoload{vv},dimtoagg(vv));
    end

    % record
    allresults{p,q,hh,rr} = [results.ang results.ecc results.expt results.rfsize results.R2 results.gain results.baseline ...
                                    circulardiff(pi/4,atan2(results.baseline,results.baseline+results.gain),2*pi)/pi*180 ...
                                    results.predR2 results.predr permute(results.params(1,[1 2 3 5],:),[3 2 1])];

  end
  
  % save a simple full .mat version!
  save(finalfile0,'allresults');

end

% DATA SAVING?
if wantsavedata

  % separate files and enforce limits
  switch mtodo
  case 1
    whtodo = 1:7;
    mns = {[] []   []    []    []    -100     []};
    mxs = {[] []   []    []    []    []       []};
    labs = {'tval' 'ang' 'mag' 'magt' 'allt' 'predR2' 'predr'};
  case 3
    whtodo = [1 2 3 4 5 6 7 8 10];
    mns = {[] []   []    []    -100   []     []  [] []};
    mxs = {[] 100  100   100   []     100    []  [] []};
    labs = {'angle' 'eccentricity' 'exponent' 'size' 'R2' 'gain' 'baseline' 'anglemetric' 'predr'};
  end

  % SAVE INTO NEUROIMAGING FILES [MULTIPLE RRs ARE SAVED TOGETHER]
  for hh=hhs
    if ~isempty(allresults{p,q,hh,1})
      for qq=1:length(whtodo)
        vals = [];
        for rii=1:size(allresults,4)
          if isempty(allresults{p,q,hh,rii})
            break;
          end
          vals(:,rii) = allresults{p,q,hh,rii}(:,whtodo(qq));
        end
        if ~isempty(mns{qq})
          vals(vals<mns{qq}) = mns{qq};
        end
        if ~isempty(mxs{qq})
          vals(vals>mxs{qq}) = mxs{qq};
        end
        fprintf('(qq=%d) min is %.4f, max is %.4f, numnotfinite=%d\n',qq,min(vals(:)),max(vals(:)),sum(~isfinite(vals(:))));
        SAVEDATA(hh,vals,sprintf('%s_%s',finalname,labs{qq}));
      end
    end
  end

end
