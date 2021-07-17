function analyzebehavior_nsd(datadir,nsdexpfile)

% history:
% - 2019/08/09 - revamp (first of series; new soak auto-repeat strategy; record 3 new columns; 250 offset)

nsdsetup;
ppdir = regexprep(datadir,'rawdata','ppdata');
mkdirquiet(ppdir);

%% BEHAVIORAL PERFORMANCE

% find behavioral files
files = matchfiles([datadir '/mat_files_from_scan/*nsd*subj*run*.mat']);

% process each file
clear results0;
badtimes = {};
for p=1:length(files)

  % do nsd analysis
  results0(p) = runnsdbehavioralanalysis(files{p},3000,[299.8 300.2],nsdexpfile);
 
  % add column
  results0(p).trialinfo(:,19) = 0;  % by default, NOT missing data
  
  % deal with missing behavioral data
  subj0 = results0(p).trialinfo(1,1);
  sess0 = results0(p).trialinfo(1,2);
  run0  = results0(p).trialinfo(1,3);
  if sessionbehaviormissing(subj0,sess0,run0)
    fprintf('*** BEHAVIOR IS MISSING! ***\n');
    results0(p).trialinfo(:,[9 10 11 15 16 17 18]) = NaN;
    results0(p).trialinfo(:,19) = 1;
    results0(p).responserate = NaN;      % NOTE THIS HANDLING
    results0(p).finalscore = NaN;
    results0(p).finalscoreWITHIN = NaN;
    results0(p).finalscoreEASY = NaN;
    results0(p).medianrt = NaN;
  end

  % do generic analysis
  a1 = load(files{p});
  deltatime = 1;
  badkey = {'5' 't'};
  deltatimeBAD = 0.1;
  [keytimes,badtimes{p},keybuttons] = ptviewmoviecheck(a1.timeframes,a1.timekeys,deltatime,badkey,deltatimeBAD,1);
  close; close;
  clear a1 keytimes keybuttons;

end

% collect outputs
clear results;
results.trialinfo =     cat(1,results0.trialinfo);
results.responserate =  cat(1,results0.responserate);
results.finalscore =    cat(1,results0.finalscore);
results.finalscoreWITHIN = cat(1,results0.finalscoreWITHIN);
results.finalscoreEASY = cat(1,results0.finalscoreEASY);
results.medianrt =      cat(1,results0.medianrt);
results.userkeys =      results0(1).userkeys;
results.userkeycounts = cat(1,results0.userkeycounts);
results.totaldur =      cat(1,results0.totaldur);
results.badtimes =      badtimes;

% sanity check that the right buttons were being pressed in each run
assert(all(sum(results.userkeycounts(:,1:2),2) > .95 * sum(results.userkeycounts,2)));

% calculate overall stuff
results.sessionresponserate = sum(~isnan(results.trialinfo(:,10))) / sum(results.trialinfo(:,19)==0) * 100;
easy1 = nansum(results.trialinfo(results.trialinfo(:,14)==1,9));
easy2 = sum(results.trialinfo(:,14)==1 & results.trialinfo(:,19)==0);
results.sessionfinalscoreEASY = easy1/easy2 * 100;

% save
save([ppdir '/behavioralresults_nsd.mat'],'-struct','results');

%% DESIGN MATRICES

% load
a0 = load(nsdexpfile,'subjectim');

% define
ptpertrials = [4 3];

% loop
for zz=1:length(ptpertrials)
  ptpertrial = ptpertrials(zz);

  % make design matrices:
  %   stimulusSINGLETRIAL is chronological single-trial
  %   stimulusFULL is fully expanded with repeats
  stimulusSINGLETRIAL = {};   % design matrix for each run (750 columms, one 1 per column)
  stimulusFULL = {};          % design matrix for each run (73k columns, one or more 1 per column)
  stimorder = zeros(1,750);   % 1 x 750 with 73k-ID
  cntSESSION = 1;
  for p=1:length(files)

    % load
    a1 = load(files{p});
  
    % do it
    stimulusSINGLETRIAL{p} = sparse(zeros(75*ptpertrial,750));     % TR x condition
    stimulusFULL{p} =        sparse(zeros(75*ptpertrial,73000));   % TR x condition
    cntRUN = 1;
    for q=1:size(a1.trialpattern,1)
      ix = find(a1.trialpattern(q,:));
      if ~isempty(ix)
        finalid = a0.subjectim(a1.setnum(2),a1.curordering(cntRUN));
        stimulusSINGLETRIAL{p}((q-1)*ptpertrial+1,cntSESSION) = 1;
        stimulusFULL{p}((q-1)*ptpertrial+1,finalid) = 1;
        stimorder(cntSESSION) = finalid;
        cntRUN = cntRUN + 1;
        cntSESSION = cntSESSION + 1;
      end
    end

  end
  if ptpertrial==4
    suffix0 = '';
  else
    suffix0 = sprintf('_%dpertrial',ptpertrial);
  end
  stimulus = stimulusSINGLETRIAL;
  save(sprintf('%s/designmatrixSINGLETRIAL%s_nsd.mat',ppdir,suffix0),'stimulus','stimorder');
  stimulus = stimulusFULL;
  save(sprintf('%s/designmatrixFULL%s_nsd.mat',ppdir,suffix0),'stimulus','stimorder');

end
