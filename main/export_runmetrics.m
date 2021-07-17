% This creates run-level metrics that are useful for assessing data quality.
% This comes at the very end.

%% %%%%% SETUP

% setup
nsdsetup;

% init
runmetrics =   NaN*zeros(8,40,12,7);  % 7 metrics are tSNR, FD, R2, RR, PCTCORRECT, EASY, EASYN
runmetricsRS = NaN*zeros(8,40,2,2);   % 2 metrics are tSNR, FD

%% %%%%% LOAD DATA

% load behavioral data
bdata = {};
for subjix=1:8
  file0 = sprintf('~/nsd/nsddata/ppdata/subj%02d/behav/responses.tsv',subjix);
  a1 = importdata(file0);
  bdata{subjix} = a1.data;
end
% bdata =
% 
%   1x8 cell array
% 
%   Columns 1 through 5
% 
%     {30000x19 double}    {30000x19 double}    {24000x19 double}    {22500x19 double}    {30000x19 double}
% 
%   Columns 6 through 8
% 
%     {24000x19 double}    {30000x19 double}    {22500x19 double}

% load motion parameters
motiondata = {};  % cell matrix 8 x 40 x 14, each element is TR x 6
for subjix=1:8
  for sessix=1:40
    file0 = sprintf('~/nsd/nsddata_timeseries/ppdata/subj%02d/func1pt8mm/motion/motion_session%02d_run*.tsv',subjix,sessix);
    files = matchfiles(file0);
    for runii=1:length(files)
      a1 = load(files{runii});
      motiondata{subjix,sessix,runii} = a1;
    end
  end
end
% >> size(motiondata)
% 
% ans =
% 
%      8    40    14

%% %%%%% ANALYZE BEHAVIORAL DATA

% process behavior
for subjix=1:8
  for sessix=1:40
    
    % if there are none in this session, skip it
    if ~any(bdata{subjix}(:,2)==sessix)
      continue;
    end

    for runix=1:12
      
      % in this session, in this run, AND is not missing data
      ix = bdata{subjix}(:,2)==sessix & bdata{subjix}(:,3)==runix & bdata{subjix}(:,19)==0;
      
      % the above AND a button was pressed
      ix2 = ix & isfinite(bdata{subjix}(:,9));
      
      % RECORD RESPONSE RATE
      runmetrics(subjix,sessix,runix,4) = sum(ix2) / sum(ix) * 100;
      
      % the above AND is an easy trial (old with respect to current session)
      ix3 = ix & bdata{subjix}(:,14)==1;
      
      % ... AND subject got it right
      ix4 = ix3 & bdata{subjix}(:,15)==1;
      
      % RECORD "EASY" PERFORMANCE
      runmetrics(subjix,sessix,runix,6) = sum(ix4) / sum(ix3) * 100;
      runmetrics(subjix,sessix,runix,7) = sum(ix3);
      
      % the above AND the subject got the trial correct
      ix5 = ix & bdata{subjix}(:,9)==1;

      % RECORD PERCENT CORRECT
      runmetrics(subjix,sessix,runix,5) = sum(ix5) / sum(ix) * 100;

    end
  end
end

%% %%%%% ANALYZE MOTION DATA

% process motion
for subjix=1:8
  for sessix=1:40
    for runix=1:14

      % get the motion data
      data0 = motiondata{subjix,sessix,runix};

      % if data were not acquired, just skip
      if ~isempty(data0)

        % compute mean FD in the run (note that one run was discontinuous)
        temp = diff(data0,[],1);
        if subjix==8 && sessix==2 && runix==4
          temp(150,:) = [];  % there was a discontinuity, so we need to ignore
        end
        fd = mean(abs(temp) * [1 1 1 50 50 50]');

        % for resting-state sessions, we have to adjust the run indexing
        isrest = (sessix>=21 && sessix<=30) || ...
                 (ismember(subjix,[1 5]) && (sessix>=31 && sessix<=38));
        if isrest
          switch runix
          case 1
            runmetricsRS(subjix,sessix,1,2) = fd;
          case 14
            runmetricsRS(subjix,sessix,2,2) = fd;
          otherwise
            runmetrics(subjix,sessix,runix-1,2) = fd;
          end
        else
          runmetrics(subjix,sessix,runix,2) = fd;
        end

      end

    end
  end
end

%% %%%%% ANALYZE TSNR

for subjix=1:8
  dirs0 = nsdallsessfun2(subjix); length(dirs0)
  for sessix=1:length(dirs0)
    ppdir = regexprep(dirs0{sessix},'rawdata','ppdata');
    
    % load tsnr
    a1 = load(sprintf('%s/autoqc_fmri.mat',ppdir),'tsnr');
    tsnr = a1.tsnr(1,:);  % raw volumes, quadratic detrend, mean/std, median over intensity valid
    if subjix==8 && sessix==2
      tsnr(5) = [];  % for this session, just drop the vestigial 5th run
    end

    % for resting-state sessions, we have to adjust the run indexing
    isrest = (sessix>=21 && sessix<=30) || ...
             (ismember(subjix,[1 5]) && (sessix>=31 && sessix<=38));
    if isrest
      runmetricsRS(subjix,sessix,1:2,1) = tsnr([1 14]);
      runmetrics(subjix,sessix,1:12,1) = tsnr(2:13);
    else
      runmetrics(subjix,sessix,1:12,1) = tsnr;
    end

  end
end

%% %%%%% ANALYZE BOLD

% define
volfun = @(x) flipdim(flipdim(permute(x,[2 1 3 4 5 6 7 8 9 10]),2),1);

% do it
for subjix=1:8
  dirs0 = nsdallsessfun2(subjix);

  % load ROI
  roi1 = load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj%02d/func1pt8mm/roi/nsdgeneral.nii.gz',subjix));

  % for each NSD session, load run-level ON-OFF R2
  for sessix=1:length(dirs0)
    glmdir = regexprep(dirs0{sessix},'rawdata','glmdataLOW');
    
    % load R2
    a1 = load(sprintf('%s/GLMdenoise_nsdBASIConoff.mat',glmdir),'R2run');
    r2run = volfun(a1.R2run);  % X x Y x Z x 12
    
    % record (median R2 across the voxels in the ROI)
    runmetrics(subjix,sessix,:,3) = nanmedian(subscript(squish(r2run,3),{find(roi1.img==1) ':'}),1);

  end
end

%% %%%%% SAVE

% save
save('~/nsd/nsddata/information/runmetrics.mat','runmetrics','runmetricsRS','-v7.3');
