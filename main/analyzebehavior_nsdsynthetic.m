function analyzebehavior_nsdsynthetic(datadir,expfile)

nsdsetup;
ppdir = regexprep(datadir,'rawdata','ppdata');
mkdirquiet(ppdir);

%% BEHAVIORAL PERFORMANCE

% find behavioral files
files = matchfiles([datadir '/mat_files_from_scan/*nsdsynthetic*subj*run*.mat'])
assert(length(files)==8);

% HACK FOR FIXING DIFFERENT BUTTONS PRESSED
if ~isempty(regexp(files{1},'NSD814')) || ~isempty(regexp(files{1},'NSD168'))
  choicebuttons = {'3' '4'};
else
  choicebuttons = {'1' '2'};
end

% process each file
clear results0;
badtimes = {};
for p=1:length(files)

  % do nsd analysis
  results0(p) = runnsdsyntheticbehavioralanalysis(files{p},4280,[427.8 428.2],nsdsyntheticexpfile,p,choicebuttons);
 
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
results.trialinfo = cat(1,results0.trialinfo);
results.userkeys =   results0(1).userkeys;
results.userkeycounts = cat(1,results0.userkeycounts);
results.totaldur =   cat(1,results0.totaldur);
results.numcorrect = cat(1,results0.numcorrect);
results.numincorrect = cat(1,results0.numincorrect);
results.pctcorrect = cat(1,results0.pctcorrect);
results.sdtmatrix =  cat(3,results0.sdtmatrix);
results.dprime =     cat(1,results0.dprime);
results.badtimes =   badtimes;

% sanity check that the right buttons were being pressed in each run
assert(all(sum(results.userkeycounts(:,cellfun(@str2double,choicebuttons)),2) > .95 * sum(results.userkeycounts,2)));

% check
min(results.totaldur)
max(results.totaldur)
results.pctcorrect
results.dprime

% final pctcorrect
results.pctcorrectTOTAL = sum(results.numcorrect) / (sum(results.numcorrect) + sum(results.numincorrect)) * 100;
results.pctcorrectTOTAL

% final dprime
sdtmatrix = sum(results.sdtmatrix,3);
hitrate = sdtmatrix(1,1) / sum(sdtmatrix(1,:));
falsealarmrate = sdtmatrix(2,1) / sum(sdtmatrix(2,:));
fprintf('FINAL: hitrate is %.2f, falsealarmrate is %.2f\n',hitrate,falsealarmrate);
A = norminv(1-max(min(hitrate,.99),.01),0,1);
B = norminv(1-max(min(falsealarmrate,.99),.01),0,1);
results.dprimeTOTAL = -A + B;
results.dprimeTOTAL

% save
save([ppdir '/behavioralresults_nsdsynthetic.mat'],'-struct','results');

%% DESIGN MATRICES

% 284-element vector with positive integers 1 to 140
maptodistinct = zeros(1,284);
maptodistinct(1:4) = 1;
maptodistinct(5:8) = 2;
maptodistinct(9:12) = 3;
maptodistinct(13:104) = 4:95;
maptodistinct(105:284) = upsamplematrix(96:140,4,2,[],'nearest');

% load
a0 = load(expfile);

% define
ptpertrials = [4 3];

% loop
for zz=1:length(ptpertrials)
  ptpertrial = ptpertrials(zz);

  % make design matrices:
  %   stimulusSINGLETRIAL is chronological single-trial
  %   stimulusFULL is fully expanded with repeats
  stimulusSINGLETRIAL = {};      % design matrix for each run (744 columms, one 1 per column)
  stimulusFULL = {};             % design matrix for each run (284 columns, one or more 1 per column)
  stimulusFULLALT = {};          % design matrix for each run (140 columns, one or more 1 per column)
  stimorder    = zeros(1,744);   % 1 x 744 with 284-ID
  stimorderALT = zeros(1,744);   % 1 x 744 with 140-ID
  cntSESSION = 1;
  for p=1:length(files)

    % load
    a1 = load(files{p});
  
    % do it
    stimulusSINGLETRIAL{p} = sparse(zeros(107*ptpertrial,744));   % TR x condition
    stimulusFULL{p} =        sparse(zeros(107*ptpertrial,284));   % TR x condition
    stimulusFULLALT{p} =     sparse(zeros(107*ptpertrial,140));   % TR x condition
    cntRUN = 1;
    for q=1:size(a1.trialpattern,1)
      ix = find(a1.trialpattern(q,:));
      if ~isempty(ix)
        finalid = a1.curordering(cntRUN);  % the 284-ID
        stimulusSINGLETRIAL{p}((q-1)*ptpertrial+1,cntSESSION) = 1;
        stimulusFULL{p}((q-1)*ptpertrial+1,finalid) = 1;
        stimulusFULLALT{p}((q-1)*ptpertrial+1,maptodistinct(finalid)) = 1;
        stimorder(cntSESSION) = finalid;
        stimorderALT(cntSESSION) = maptodistinct(finalid);
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
  save(sprintf('%s/designmatrixSINGLETRIAL%s_nsdsynthetic.mat',ppdir,suffix0),'stimulus','stimorder','stimorderALT');
  stimulus = stimulusFULL;
  save(sprintf('%s/designmatrixFULL%s_nsdsynthetic.mat',ppdir,suffix0),'stimulus','stimorder','stimorderALT');
  stimulus = stimulusFULLALT;
  save(sprintf('%s/designmatrixFULLALT%s_nsdsynthetic.mat',ppdir,suffix0),'stimulus','stimorder','stimorderALT');

end
