function results = runnsdbehavioralanalysis(filename,expectedframes,expectedduration,nsdexpfile,choicebuttons);

% function results = runnsdbehavioralanalysis(filename,expectedframes,expectedduration,nsdexpfile,choicebuttons);
%
% <filename> is the .mat file saved from runnsd.m
% <expectedframes> is the expected number of frames in the experiment
% <expectedduration> is [A B] with the expected range of the total experiment duration in seconds
% <nsdexpfile> is the experimental design .mat file that controls the nsd experiment
% <choicebuttons> (optional) is {C1 C2} such that C1 means novel and C2 means old.
%   Default: {'1' '2'}.
%
% Check timing and analyze the behavioral data.
% We report information to the command window and
% also return <results>.
%
% <results> is a struct with fields:
%   <trialinfo> is N x 18 where N is the number of stimulus trials in the run.
%     column 1 is subject number (1-8)
%     column 2 is session number (1-40)
%     column 3 is run number (1-12)
%     column 4 is stimulus trial number (1-63 (odd runs) or 1-62 (even runs))
%     column 5 is the 73k ID of the presented image
%     column 6 is the 10k ID of the presented image
%     column 7 is the trial start time (i.e. time that image comes on) as 
%                  a MATLAB serial date number. the units are days.
%     column 8 is 0 (novel) or 1 (old)
%     column 9 is 0 (incorrect) or 1 (correct)
%     column 10 is reaction time in milliseconds (time between trial start time
%                  and button-press time)
%     column 11 is whether this is a trial that involved more than one button press
%                  (0 = no, 1 = yes, NaN = no buttons pressed)
%     column 12 is number of stimulus trials in between current
%                  and most recent presentation
%     column 13 is number of stimulus trials in between current
%                  and second most recent presentation. if there
%                  has been only one previous presentation, this is NaN.
%     column 14 is 0 (novel) or 1 (old) with respect to the current session only
%     column 15 is 0 (incorrect) or 1 (correct) with respect to the current session only
%     column 16 is the total number of 1s ("novel") pressed during this trial
%     column 17 is the total number of 2s ("old") pressed during this trial
%     column 18 is the final button that we officially score (this will be either 1, 2, or NaN)
%     note that columns 12-13 are NaN for the case of novel images.
%     note that columns 9-11 and 15 and 18 are NaN if no button was pressed for this trial
%   <responserate> is percent of stimulus trials with at least a button press
%   <finalscore> is percent of stimulus trials with a correct response
%   <finalscoreWITHIN> is percent of stimulus trials with a correct response,
%     acting as if the experiment spanned only the current session
%   <finalscoreEASY> is percent correct for easy trials (only). easy is defined as
%     trials where it is a memory event and an earlier instance occurred sometime
%     during the current session.
%   <medianrt> is the median reaction time in milliseconds
%   <userkeys> is a cell vector of possible user keys
%   <userkeycounts> is a vector of number of times that the keys were pressed
%   <totaldur> is the total empirical duration of the experiment
%
% Note that we detect user buttons from the time that the stimulus comes on
% until the time that the next trial starts (or would have started in the
% case of blank trials), which is a 4-s interval, except that we offset this
% by 250 ms.
%
% We process only the final button detected (if there is more than one), but
% specifically target the first instance of that one (if there is a series).
% 
% We look specifically for buttons named in <choicebuttons>.
%
% Note that this function is not guaranteed to correctly process partial runs!!
%
% Note that outside of this function we handle specially the case where runs were
% conducted, but due to computer error no buttons were actually saved!!

% history:
% - 2019/08/09 - revamp (first of series; new soak auto-repeat strategy; record 3 new columns; 250 offset; fix simultaenous 1/2 problems)

% input
if ~exist('choicebuttons','var') || isempty(choicebuttons)
  choicebuttons = {'1' '2'};
end

%% Simple frame and duration stuff

% load the data
a1 = load(filename);

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

%% Now analyze the behavior

% internal constants
deltatime = 20;    % buttons are automatically extended if this many milliseconds elapse while still the same
validkeys = {'1!' '2@' '3#' '4$' '5%' 'r' 'y' 'g' 'b' 't' 'absolutetimefor0' 'trigger'};
userkeys = {'1' '2' '3' '4' 'r' 'y' 'g' 'b'};
windowoffset = 250;  % milliseconds of offset. decide to get buttons from [trialstart+offset,trialend+offset]

% load
a2 = load(nsdexpfile);

% expand the multiple-keypress cases [create timekeysB]
timekeysB = {};
for p=1:size(a1.timekeys,1)
  if iscell(a1.timekeys{p,2})
    for pp=1:length(a1.timekeys{p,2})

      % in the simultaneous keypress case, if 3 or 4 is one of them, silently ignore
      if ismember(a1.timekeys{p,2}{pp},{'3#' '4$'})
      else
        timekeysB{end+1,1} = a1.timekeys{p,1};
        timekeysB{end,2} = a1.timekeys{p,2}{pp};
      end

    end
  else
    timekeysB(end+1,:) = a1.timekeys(p,:);
  end
end

% figure out when the user pressed a button [buttontimes]
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

% figure out absolute time for first frame
matlabnowtime = timekeysB{find(ismember(timekeysB(:,2),'absolutetimefor0')),1};

% calc
subjnum = a1.setnum(2);
sessnum = a1.setnum(3);
runnum =  a1.setnum(4);
trialnum = a1.setnum(5);

% sanity checks
assert(subjnum>=1 && subjnum<=8);
assert(sessnum>=1 && sessnum<=40);
assert(runnum>=1  && runnum<=12);
assert(trialnum>=0 && trialnum<=75);
if trialnum > 4
  fprintf('\nWARNING!! THIS IS A PARTIAL RUN!! RESULTS MAY BE INVALID!!\n\n');
end

% initialize results
clear results;

% process each trial
cnt = 1;  % this keeps track of the stimulus trial number
for p=1:size(a1.trialpattern)
  temp = find(a1.trialpattern(p,:));
  
  % if this is a stimulus trial
  if ~isempty(temp)
  
    % index of current stimulus trial (relative to 30k)
    ix = (sessnum-1)*750 + sum(a1.numstimperrun(1:runnum-1)) + cnt;
    
    % index of first stimulus trial in this session
    ixthissession = (sessnum-1)*750 + 1;
    
    % record basic information
    results.trialinfo(cnt,[1 2 3 4]) = [subjnum sessnum runnum cnt];
    
    % is current image one of the previously shown ones? (0 means novel, 1 means old)
    results.trialinfo(cnt,8) = ismember(a2.masterordering(ix),a2.masterordering(1:ix-1));

    % is current image one of the previously shown ones in current session? (0 means novel, 1 means old)
    results.trialinfo(cnt,14) = ismember(a2.masterordering(ix),a2.masterordering(ixthissession:ix-1));
    
    % record 73k ID
    results.trialinfo(cnt,5) = a2.subjectim(subjnum,a2.masterordering(ix));
    
    % record 10k ID
    results.trialinfo(cnt,6) = a2.masterordering(ix);
    
    % if current image is novel
    if results.trialinfo(cnt,8)==0

      results.trialinfo(cnt,[12 13]) = NaN;

    else

      % indices of previous stimulus trials (relative to 30k) [there will be either 1 or 2]
      ix2 = find(a2.masterordering(1:ix-1) == a2.masterordering(ix));
      
      % how many stimulus trials separated from the most recent presentation?
      results.trialinfo(cnt,12) = ix - lastel(ix2) - 1;
      
      % how many stimulus trials separated from the first presentation?
      if length(ix2) > 1
        results.trialinfo(cnt,13) = ix - firstel(ix2) - 1;
      else
        results.trialinfo(cnt,13) = NaN;
      end

    end
    
    % calc some times
    trialstart = a1.timeframes((p-1)*40 + 1);       % time of first frame for this stimulus trial
    trialend =   a1.timeframes((p-1)*40 + 1 + 40);  % time of first frame for the next trial
    
    % find the buttons pressed during this time interval
    ok = find((buttontimes > trialstart+windowoffset/1000) & (buttontimes <= trialend+windowoffset/1000));

    % only want buttons that are one of the choice buttons
    good = logical(zeros(1,length(ok)));
    for qq=1:length(ok)
      good(qq) = ismember(buttonpressed{ok(qq)},choicebuttons);
    end
    ok = ok(good);

    % if no buttons detected
    if isempty(ok)

      results.trialinfo(cnt,[9 10 11 15]) = NaN;
      results.trialinfo(cnt,[16 17]) = 0;
      results.trialinfo(cnt,[18]) = NaN;

    % otherwise, record RT (ms) and whether it was correct (0 or 1) in both versions
    else
      
      % if there is more than one unique button pressed, mark it as "change mind"
      if length(union(buttonpressed(ok),[])) > 1
        results.trialinfo(cnt,11) = 1;
      else
        results.trialinfo(cnt,11) = 0;
      end
      
      % tally all the buttons pressed
      results.trialinfo(cnt,16) = sum(ismember(buttonpressed(ok),choicebuttons{1}));
      results.trialinfo(cnt,17) = sum(ismember(buttonpressed(ok),choicebuttons{2}));

      % handle the case of multiple buttons
      cnt0 = length(ok);
      while 1
        if cnt0-1 < 1
          break;  % if we are at the beginning, get out
        end
        if isequal(buttonpressed{ok(cnt0-1)},buttonpressed{ok(cnt0)})
          cnt0 = cnt0 - 1;  % if cnt0-1 is the same button, use that one!
        else
          break;  % if cnt0-1 is a different button, don't go any more into the past
        end
      end
      ok = ok(cnt0);  % score this button press (out of potentially many button presses)

      % calc
      timeofpress = buttontimes(ok);  % time that key was pressed
      thekey = buttonpressed{ok};     % the key that was pressed, e.g. '1'

      % RT (ms)
      results.trialinfo(cnt,10) = (timeofpress - trialstart) * 1000;
      
      % correctness (0 or 1)
      results.trialinfo(cnt,9) = ((results.trialinfo(cnt,8)==0) && (isequal(thekey,choicebuttons{1}))) || ...  % novel and '1' was pressed
                                 ((results.trialinfo(cnt,8)==1) && (isequal(thekey,choicebuttons{2})));        % old and '2' was pressed

      % correctness (0 or 1) considering this session only
      results.trialinfo(cnt,15) = ((results.trialinfo(cnt,14)==0) && (isequal(thekey,choicebuttons{1}))) || ...  % novel and '1' was pressed
                                  ((results.trialinfo(cnt,14)==1) && (isequal(thekey,choicebuttons{2})));        % old and '2' was pressed
      
      % record the button we actually scored
      if isequal(thekey,choicebuttons{1})
        results.trialinfo(cnt,18) = 1;
      else
        assert(isequal(thekey,choicebuttons{2}));
        results.trialinfo(cnt,18) = 2;
      end
      
    end
    
    % what is the absolute time that the trial started?
    results.trialinfo(cnt,7) = matlabnowtime + trialstart/60/60/24;  % convert from seconds to days and then add in
    
    % increment the stimulus-trial counter
    cnt = cnt + 1;
  
  end
  
end

% calculate some summary metrics
results.responserate = sum(~isnan(results.trialinfo(:,10))) / size(results.trialinfo,1) * 100;  % percent stimulus trials with button press
results.finalscore =       nansum(results.trialinfo(:,9)) / size(results.trialinfo,1) * 100;   % percent correct
results.finalscoreWITHIN = nansum(results.trialinfo(:,15)) / size(results.trialinfo,1) * 100;  % percent correct
easy1 = nansum(results.trialinfo(results.trialinfo(:,14)==1,9));
easy2 = sum(results.trialinfo(:,14)==1);
results.finalscoreEASY = easy1/easy2 * 100;  % percent correct
results.medianrt = nanmedian(results.trialinfo(:,10));

% report buttons pressed (for sanity checking)
fprintf('Number of pressed buttons:\n');
userkeycounts = [];
for pp=1:length(userkeys)
  userkeycounts(pp) = sum(ismember(buttonpressed,userkeys{pp}));
  fprintf('%s: %d\n',userkeys{pp},userkeycounts(pp));
end
results.userkeys = userkeys;
results.userkeycounts = userkeycounts;

% store some info
results.totaldur = totaldur;

% report scores
fprintf('==============================================================\n');
fprintf(['Your response rate is %d%%.\n' ...
         'Your median RT is %d ms.\n' ...
         '\n' ...
         'OVERALL (WRT ALL SCAN SESSIONS):\n' ...
         '  Out of %d novel stimulus trials, you had %d correct, %d incorrect, and %d no-response.\n' ...
         '  Out of %d old stimulus trials, you had %d correct, %d incorrect, and %d no-response.\n' ...
         '  Your overall percent correct is %d%%.\n' ...
         '\n' ...
         'SHORT-TERM (WRT THIS SCAN SESSION ONLY):\n' ...
         '  Out of %d novel stimulus trials, you had %d correct, %d incorrect, and %d no-response.\n' ...
         '  Out of %d old stimulus trials, you had %d correct, %d incorrect, and %d no-response.\n' ...
         '  Your overall percent correct is %d%%.\n' ...
         '\n' ...
         'EASY TRIALS percent correct (%d out of %d) is %.3f%%.\n' ...
         ], ...
        round(results.responserate), ...
        round(results.medianrt), ...
        ...
        sum(results.trialinfo(:,8)==0), ...
        sum(results.trialinfo(results.trialinfo(:,8)==0,9)==1), ...
        sum(results.trialinfo(results.trialinfo(:,8)==0,9)==0), ...
        sum(isnan(results.trialinfo(results.trialinfo(:,8)==0,9))), ...
        sum(results.trialinfo(:,8)==1), ...
        sum(results.trialinfo(results.trialinfo(:,8)==1,9)==1), ...
        sum(results.trialinfo(results.trialinfo(:,8)==1,9)==0), ...
        sum(isnan(results.trialinfo(results.trialinfo(:,8)==1,9))), ...
        round(results.finalscore), ...
        ...
        sum(results.trialinfo(:,14)==0), ...
        sum(results.trialinfo(results.trialinfo(:,14)==0,15)==1), ...
        sum(results.trialinfo(results.trialinfo(:,14)==0,15)==0), ...
        sum(isnan(results.trialinfo(results.trialinfo(:,14)==0,15))), ...
        sum(results.trialinfo(:,14)==1), ...
        sum(results.trialinfo(results.trialinfo(:,14)==1,15)==1), ...
        sum(results.trialinfo(results.trialinfo(:,14)==1,15)==0), ...
        sum(isnan(results.trialinfo(results.trialinfo(:,14)==1,15))), ...
        round(results.finalscoreWITHIN), ...
        ...
        easy1,easy2,results.finalscoreEASY);
fprintf('==============================================================\n');
