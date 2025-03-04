function results = runnsdsyntheticbehavioralanalysis(filename,expectedframes,expectedduration,nsdexpfile,crunnum,choicebuttons);

% function results = runnsdsyntheticbehavioralanalysis(filename,expectedframes,expectedduration,nsdexpfile,crunnum,choicebuttons);
%
% <filename> is the .mat file saved from runnsdsynthetic.m
% <expectedframes> is the expected number of frames in the experiment
% <expectedduration> is [A B] with the expected range of the total experiment duration in seconds
% <nsdexpfile> is the experimental design .mat file that controls the nsdsynthetic experiment
% <crunnum> is 1-8 with the chronological run number (odds are fixation task, evens are memory task)
% <choicebuttons> (optional) is {C1 C2} such that C1 means the first button and C2 means the second button.
%   Default: {'1' '2'}.
%
% Check timing and analyze the behavioral data.
% We report information to the command window and
% also return <results>.
%
% <results> is a struct with fields:
%   <trialinfo> is N x 5 where N is the number of stimulus trials in the run.
%     column 1 is subject number (1-8)
%     column 2 is chronological run number (1-8)
%     column 3 is stimulus trial number (1-93)
%     column 4 is the 284-ID of the presented image
%     column 5 is the trial start time (i.e. time that image comes on) as 
%                  a MATLAB serial date number. the units are days.
%   <userkeys> is a cell vector of possible user keys
%   <userkeycounts> is a vector of number of times that the keys were pressed
%   <totaldur> is the total empirical duration of the experiment
%
% in the case of the fixation task, <results> will also have:
%   <numcorrect> as the number of correct
%   <numincorrect> as the number of incorrect
%   <pctcorrect> as the percent correct of luminance change events
%
% in the case of the memory task, <results> will also have:
%   <sdtmatrix> as 2 x 2 with [HIT MISS; FA CR]
%   <dprime> with the d-prime performance
%
% the above fields will be [] if not applicable.
%
% See code for details.

%% Setup

% input
if ~exist('choicebuttons','var') || isempty(choicebuttons)
  choicebuttons = {'1' '2'};
end

% initialize results
clear results;

% load the data
a1 = load(filename);     % the behavioral data file
a2 = load(nsdexpfile);   % the experimental design file

%% Simple timeframes and duration stuff

% calc
totaldur = mean(diff(a1.timeframes)) * length(a1.timeframes);

% prettify
fprintf('==============================================================\n');

% report the timing
fprintf('Experiment duration was %.3f.\n',totaldur);

% check the timing
if length(a1.timeframes)==expectedframes && ...
   totaldur > expectedduration(1) && totaldur < expectedduration(2)
  fprintf('Timing was fine!\n');
else
  fprintf('ERROR! Timing was off!!!!');
end

% store
results.totaldur = totaldur;

%% Button pre-processing stuff

% internal constants
deltatime = 20;    % buttons are automatically extended if this many milliseconds elapse while still the same
validkeys = {'1!' '2@' '3#' '4$' '5%' 'r' 'y' 'g' 'b' 't' 'absolutetimefor0' 'trigger'};
userkeys = {'1' '2' '3' '4' 'r' 'y' 'g' 'b'};

% expand the multiple-keypress cases [create timekeysB]
timekeysB = {};
for p=1:size(a1.timekeys,1)
  if iscell(a1.timekeys{p,2})
    for pp=1:length(a1.timekeys{p,2})
      timekeysB{end+1,1} = a1.timekeys{p,1};
      timekeysB{end,2} = a1.timekeys{p,2}{pp};
    end
  else
    timekeysB(end+1,:) = a1.timekeys(p,:);
  end
end

% figure out absolute time for first frame
matlabnowtime = timekeysB{find(ismember(timekeysB(:,2),'absolutetimefor0')),1};

% clean up all the button stuff [create buttontimes and buttonpressed]
oldkey = ''; oldkeytime = -Inf;
buttontimes = [];
buttonpressed = {};
for p=1:size(timekeysB,1)

  % warn if weird key found
  if ~ismember(timekeysB{p,2},validkeys)
    fprintf('*** Unknown key detected (%s); ignoring.\n',timekeysB{p,2});
    continue;
  end

  % is this a bogus key press?
  bad = isequal(timekeysB{p,2},'absolutetimefor0') | ...
        isequal(timekeysB{p,2},'trigger') | ...
        isequal(timekeysB{p,2}(1),a1.triggerkey);
  
  % is this a "held down" case?
  bad2 = (isequal(timekeysB{p,2},oldkey) & timekeysB{p,1}-oldkeytime <= deltatime/1000);

  % if it appears to be a new key, we should do a special case check for simultaneity.
  % bad3 indicates if we should ignore the current key because it comes while some 
  % other key was originally held down first.
  bad3 = 0;
  if ~bad && ~isequal(timekeysB{p,2},oldkey)
  
    % scan ahead...
    q = p+1;
    while 1
      if q>size(timekeysB,1)
        break;  % if we run out, just break
      end
      if timekeysB{q,1} > oldkeytime + deltatime/1000
        break;  % if we are past the window, just break
      end
      if isequal(timekeysB{q,2},oldkey)
        bad3 = 1;  % if we are within the deltatime AND it is the same as the old key, then mark the current one for ignoral
        break;
      end
      q = q + 1;
    end

  end
  
  % if still a held-down button, extend the time
  if bad2
    oldkeytime = timekeysB{p,1};
  end

  % if not bogus, then record the button time
  if ~(bad | bad2 | bad3)
    buttontimes = [buttontimes timekeysB{p,1}];
    buttonpressed = [buttonpressed {timekeysB{p,2}(1)}];
    oldkey = timekeysB{p,2};
    oldkeytime = timekeysB{p,1};
  end

end

% report buttons pressed to the screen (for sanity checking!!)
fprintf('Number of pressed buttons:\n');
userkeycounts = [];
for pp=1:length(userkeys)
  userkeycounts(pp) = sum(ismember(buttonpressed,userkeys{pp}));
  fprintf('%s: %d\n',userkeys{pp},userkeycounts(pp));
end

% store
results.userkeys = userkeys;
results.userkeycounts = userkeycounts;

%% Now analyze the behavior!

% calc
subjnum = a1.setnum(2);  % 1-8
runnum =  a1.setnum(3);  % 1-4
stimids = a2.masterordering((crunnum-1)*93 + (1:93));  % sequence of stimulus ids (which range from 1-284)

% sanity checks
assert(subjnum>=1 && subjnum<=8);
assert(runnum>=1  && runnum<=4);

% init
results.numcorrect = [];
results.numincorrect = [];
results.pctcorrect = [];
results.sdtmatrix = [];
results.dprime = [];

%% COMMON STIMULUS TRIAL STUFF

cnt = 1;  % this keeps track of the stimulus trial number
for p=1:size(a1.trialpattern)
  temp = find(a1.trialpattern(p,:));
  
  % if this is a stimulus trial
  if ~isempty(temp)
  
    % record basic information
    results.trialinfo(cnt,[1 2 3 4]) = [subjnum crunnum cnt stimids(cnt)];
        
    % calc some times
    trialstart = a1.timeframes((p-1)*40 + 1);
    
    % what is the absolute time that the trial started?
    results.trialinfo(cnt,5) = matlabnowtime + trialstart/60/60/24;  % convert from seconds to days and then add in
    
    % increment the stimulus-trial counter
    cnt = cnt + 1;
  
  end
  
end

%% Handle the two different tasks

switch mod(crunnum,2)

% this is the fixation-task case
case 1

  % general approach:
  % responsewindow is assumed (0 to 1400 ms after dot color change event)
  % when the dot color actually changes (in some instances the color stays the same!), that triggers an actual TRIAL.
  % for each actual TRIAL:
  %   take first button pressed in its response window and grade it (correct or incorrect).
  %     if there are any other buttons, just leave them alone.
  %   if no button, it counts as incorrect.
  % if any dot color change TRIAL has its full response window that is not within the very last stimulus frame, just ignore.

  % internal constants
  responsewindow = [0 1400+0];  % get buttons from this range of milliseconds after the fixation dot color change

  % calc
  df = diff(a1.fixationorder(2:end-2));  % vector of color change amounts (NOTE THIS HARD CODING of 2:end-2). positive means got darker. negative means got brighter.
  seq = find(df~=0) + 1;   % vector of indices indicating frames on which the dot color changed

  % process each color change
  numcorrect = 0;
  numincorrect = 0;
  for p=1:length(seq)
    windowstart = a1.timeframes(seq(p)) + responsewindow(1)/1000;
    windowend   = a1.timeframes(seq(p)) + responsewindow(2)/1000;

    % if we don't have the full window, just skip
    if windowend > a1.timeframes(end)
      continue;
    end
  
    % find buttons within response window. we only want buttons that are one of the choice buttons.
    okok = find( buttontimes > windowstart & ...
                 buttontimes <= windowend & ...
                 ismember(buttonpressed,choicebuttons) );

    % if no valid buttons found, this is incorrect
    if isempty(okok)
      numincorrect = numincorrect + 1;

    % if we have at least one valid button
    else
  
      % grade the button
      if (df(seq(p)-1) < 0 & (buttonpressed{okok(1)} == choicebuttons{2})) || ...  % if negative/brighter and correct button
         (df(seq(p)-1) > 0 & (buttonpressed{okok(1)} == choicebuttons{1}))         % or positive/darker and correct button
        numcorrect = numcorrect + 1;
      else
        numincorrect = numincorrect + 1;
      end
    
    end
  end

  % record
  results.numcorrect = numcorrect;
  results.numincorrect = numincorrect;
  results.pctcorrect = numcorrect / (numcorrect + numincorrect) * 100;

  % report scores
  fprintf('==============================================================\n');
  fprintf(['Percent correct is %.2f%%.\n'],results.pctcorrect);
  fprintf('==============================================================\n');

case 0

  % general approach:
  % start grading from 2nd stim trial through last stim trial.
  % (the first stim trial doesn't really have an answer, so it's like NaN...)
  % subject has 0.25 to 4.25 s to respond
  % every stim has a correct answer!
  % a repeat is akin to a "signal" and if you say "old" it is a hit.
  % a novel image is akin to "no signal" and if you say "new" it is correct rejection.
  % if no response, that is counted as pressing the wrong answer.
  
  % internal constants
  responsewindow = [250 4250];  % get buttons from this range of milliseconds after stimulus onset

  % process each trial
  cnt = 1;  % this keeps track of the stimulus trial number
  sdtmatrix = zeros(2,2);  % row 1: hit miss, row 2: FA CR
  for p=1:size(a1.trialpattern)
    temp = find(a1.trialpattern(p,:));
  
    % if this is a stimulus trial
    if ~isempty(temp)

      % if first trial, just ignore
      if cnt==1
        cnt = cnt + 1;
        continue;
      end
      
      % calc some times
      windowstart = a1.timeframes((p-1)*40 + 1) + responsewindow(1)/1000;
      windowend =   a1.timeframes((p-1)*40 + 1) + responsewindow(2)/1000;
    
      % find buttons within response window. we only want buttons that are one of the choice buttons.
      okok = find( buttontimes > windowstart & ...
                   buttontimes <= windowend & ...
                   ismember(buttonpressed,choicebuttons) );
      
      % is this a repeat?
      isrep = stimids(cnt)==stimids(cnt-1);
      
      % if no buttons pressed
      if isempty(okok)
        if isrep
          sdtmatrix(1,2) = sdtmatrix(1,2) + 1;  % stim was a repeat, subject missed it
        else
          sdtmatrix(2,1) = sdtmatrix(2,1) + 1;  % stim was a novel one, subject false alarmed
        end
        
      % if at least one button pressed
      else
        if isrep
          if buttonpressed{okok(1)} == choicebuttons{1}  % stim was a repeat, subject said "new"
            sdtmatrix(1,2) = sdtmatrix(1,2) + 1;  % miss
          else
            sdtmatrix(1,1) = sdtmatrix(1,1) + 1;  % hit
          end
        else
          if buttonpressed{okok(1)} == choicebuttons{1}  % stim was a novel one, subject said "new"
            sdtmatrix(2,2) = sdtmatrix(2,2) + 1;  % CR
          else
            sdtmatrix(2,1) = sdtmatrix(2,1) + 1;  % FA
          end
        end
      end
    
      % increment the stimulus-trial counter
      cnt = cnt + 1;
  
    end
  
  end
  assert(cnt==(93+1));

  % compute d'
  hitrate = sdtmatrix(1,1) / sum(sdtmatrix(1,:));
  falsealarmrate = sdtmatrix(2,1) / sum(sdtmatrix(2,:));
  fprintf('hitrate is %.2f, falsealarmrate is %.2f\n',hitrate,falsealarmrate);
  A = norminv(1-max(min(hitrate,.99),.01),0,1);
  B = norminv(1-max(min(falsealarmrate,.99),.01),0,1);

  % record
  results.sdtmatrix = sdtmatrix;
  results.dprime = -A + B;

  % report scores
  fprintf('==============================================================\n');
  fprintf('full 2x2 matrix is %s.\n',mat2str(sdtmatrix));
  fprintf(['d-prime is %.2f.\n'],results.dprime);
  fprintf('==============================================================\n');

end
