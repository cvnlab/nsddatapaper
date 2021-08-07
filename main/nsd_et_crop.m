%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIRST, DEAL WITH VIDEO FEED

% files are in nsd/rawdata (both .edf and .mp4 files)

% ignore the following:
% eyevideo-20191013-NSD134-nsdsyntheticpilot.mp4
% eyevideo-20191103-NSD134-nsdimageryFAILED.mp4

% .edf recorded only for actual NSD (experiment) runs, not resting state!
% video clips available only for NSD runs.

% went through excel log and read it.

% NOPE video feed: remove eyevideo-20190822-NSD134-nsd32.mp4
% NOPE video feed white noise: eyevideo-20191110-NSD139-nsdimagery
% NOPE video feed white noise: eyevideo-20191110-NSD400-nsdimagery
% NOPE video feed corruption: eyevideo-20191113-NSD168-nsdimagery

% eyetracking data (.edf files) are present in 2+2+0+2+2+2+2+0 + 2*8 = 28 sessions
% video feed (.mp4 files) are present in 2+1+2+2+2+1+1+2 + (2*8-3) = 26 sessions

% create a cleaned version "eyetrackingtrim.xlsx" for machine readability.
%  search for 05.41 -> 05:41 [typo fix]

% empirical audio check:
% - 31.8 s from buzz to start of first beep
% - total beeping lasts 5m 8.88s = 308.88 s
% - actual data duration is 188*1.6 = 300.8 s
% - let's use: 31.8 (calibration) + 8 (dummy) + 300.8 (actual)

%%%%%%%%%%%

nsdsetup;

[num,txt,raw] = xlsread('eyetrackingtrim.xlsx');
sessix = find(~cellfun(@isempty,txt(:,1)));
for p=1:length(sessix)
  beginix = sessix(p);
  if p==length(sessix)
    endix = size(txt,1);
  else
    endix = sessix(p+1)-1;
  end
  
  filename0 = txt{beginix,1};
  totnum = endix-beginix+1;  % total number of runs to clip
  
  temp = regexp(filename0,'.+(?<subj>NSD.+?)-nsd(?<sess>.+)','names');
  subjix = find(ismember(allnsdids,temp.subj)); assert(~isempty(subjix));
  if isequal(temp.sess,'imagery')
    sessname = 'nsdimagery';
    offset = 0;
  elseif isequal(temp.sess,'synthetic')
    sessname = 'nsdsynthetic';
    offset = 0;
  else
    sessname = sprintf('session%s',temp.sess);
    s0 = str2num(temp.sess);
    if ismember(subjix,[1 5])
      if s0>=21 && s0<=38
        offset = 1;
      else
        offset = 0;
      end
    else
      if s0>=21 && s0<=30
        offset = 1;
      else
        offset = 0;
      end
    end
  end

  for q=1:totnum
    starttime0 = txt{beginix-1+q,3};  % '00:20:27.915'
    dur0 =       txt{beginix-1+q,5};  % '05:41.423'
    
    % rough durations are:
    %05:41.556
    %07:48.807
    %04:39.899
    %08:39.837
    
    % here are expected durations:
    expectedduration = [188 268 150 300]*1.6 + 39.644;
    expectedduration2 = [188 268 150 300]*1.6;
    
    % idealized times for total actual fMRI data
    goodtimes = {'00:05:00.8' '00:07:08.8' '00:04:00' '00:08:00'};
    
    % which run duration are we in
    dur0sec = sum(cellfun(@str2num,strsplit(dur0,':')) .* [60 1]);
    [mn,eix] = min(abs(dur0sec-expectedduration));
    assert(mn <= 10);  % sanity check
    
    % define
    dir0 = sprintf('/Users/kendrick/Desktop/nsddata_timeseries/ppdata/subj%02d/eyevideo/',subjix);
    mkdirquiet(dir0);

    if 0

      % here we strip until 1 second before the start of the beeping
      startsec = sum(cellfun(@str2num,strsplit(starttime0,':')) .* [60*60 60 1]) + 31.8 - 1;
      extrastr = '';

      % clip it!
      hh0 = floor(startsec/(60*60));
      mm0 = floor((startsec - (hh0*60*60))/60);
      ss0 = (startsec - (hh0*60*60) - mm0*60);
      unix_wrapper(sprintf('ffmpeg -y -i /Volumes/Lana/eyevideo/%s.mp4 -ss %02d:%02d:%.4f -t %s -c copy %s %s/eyevideo_%s_run%02d.mp4', ...
                           filename0,hh0,mm0,ss0,goodtimes{eix},extrastr,dir0,sessname,offset+q));

    end

    % here we start exactly with actual fmri data
    startsec = sum(cellfun(@str2num,strsplit(starttime0,':')) .* [60*60 60 1]) + 31.8 + 8;
    extrastr = '-an';

    % clip it!
    endsec = startsec + expectedduration2(eix);
    hh0 = floor(startsec/(60*60));
    mm0 = floor((startsec - (hh0*60*60))/60);
    ss0 = (startsec - (hh0*60*60) - mm0*60);
    hh1 = floor(endsec/(60*60));
    mm1 = floor((endsec - (hh1*60*60))/60);
    ss1 = (endsec - (hh1*60*60) - mm1*60);
    unix_wrapper(sprintf('ffmpeg -y -i /Volumes/Lana/eyevideo/%s.mp4 -ss %02d:%02d:%.4f -to %02d:%02d:%.4f -c copy %s %s/eyevideo_%s_run%02d.mp4', ...
                         filename0,hh0,mm0,ss0,hh1,mm1,ss1,extrastr,dir0,sessname,offset+q));
  end
end

%%%%%%%%%%%

% first pass:
% now, check that each .mp4 has about 1 second until beep.
% if not, adjust the start time in the xlsx doc and save.

% second pass:
% just repeat it and spot check.

% third pass:
% crop to actual fmri data collection, make it actual length, use precise end time, remove audio

% look over files.

% get to surly:
% rsync -av nsddata_timeseries kendrick@surly.cmrr.umn.edu:"/home/surly-raid3/kendrick-data/nsd/"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEXT DEAL WITH .EDF FILES

nsdsetup;

dirs0 = matchfiles('~/nsd/rawdata/*/eyetracking');

for p=1:length(dirs0)
  if ~isempty(regexp(dirs0{p},'nsdsyntheticpilot'))
    continue;
  end
  if ~isempty(regexp(dirs0{p},'nsdimageryFAILED'))
    continue;
  end
  
  files0 = matchfiles([dirs0{p} '/*.edf'])
  length(files0)
  
  if length(files)==0
    continue;
  end

  temp = regexp(dirs0{p},'.+(?<subj>NSD.+?)-nsd(?<sess>.+?)/','names');
  subjix = find(ismember(allnsdids,temp.subj)); assert(~isempty(subjix));
  if isequal(temp.sess,'imagery')
    sessname = 'nsdimagery';
    offset = 0;
  elseif isequal(temp.sess,'synthetic')
    sessname = 'nsdsynthetic';
    offset = 0;
  else
    sessname = sprintf('session%s',temp.sess);
    s0 = str2num(temp.sess);
    if ismember(subjix,[1 5])
      if s0>=21 && s0<=38
        offset = 1;
      else
        offset = 0;
      end
    else
      if s0>=21 && s0<=30
        offset = 1;
      else
        offset = 0;
      end
    end
  end

  % define
  dir0 = sprintf('/home/surly-raid3/kendrick-data/nsd/nsddata_timeseries/ppdata/subj%02d/eyedata/',subjix);
  mkdirquiet(dir0);
  
  % copy
  for q=1:length(files0)
    copyfile(files0{q},sprintf('%s/eyedata_%s_run%02d.edf',dir0,sessname,offset+q));
  end

end
