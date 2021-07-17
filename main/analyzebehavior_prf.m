function analyzebehavior_prf(datadir)

ppdir = regexprep(datadir,'rawdata','ppdata');
mkdirquiet(ppdir);

% define some constants
deltatime = 1;    % holding the key down for less than this time will be counted as one button press
badkey = {'5' 't'};
deltatimeBAD = 0.1;

% load behavioral files
stimulusfiles = matchfiles([datadir '/mat_files_from_scan/*subj*run*.mat']);

% analyze the data
[results,resultsrt] = calcfixationdotchangeperformance(stimulusfiles);

% loop over files
totaldur = [];
badtimes = {};
for zz=1:length(stimulusfiles)

  % load
  a1 = load(stimulusfiles{zz});
  
  % calc
  totaldur(zz) = mean(diff(a1.timeframes)) * length(a1.timeframes);

  % do ptviewmoviecheck
  [keytimes,badtimes{zz},keybuttons] = ptviewmoviecheck(a1.timeframes,a1.timekeys,deltatime,badkey,deltatimeBAD,1);
  close; close;

end

% save
save([ppdir '/behavioralresults_prf.mat'],'results','resultsrt','totaldur','badtimes');
