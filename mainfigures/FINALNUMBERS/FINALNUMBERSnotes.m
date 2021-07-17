% compute
a1 = load('/Users/kendrick/Dropbox/nsddata/experiments/nsd/nsd_expdesign.mat');
wh = [40 40 32 30 40 32 40 30];
unum = [];
hnum = [];
rec = [];   % 8 subjects x 10000 with number of trials for that image
commonset = [];
totaln = 0;
totaltrials = 0;
allshownims = [];
for p=1:8  % loop over subject
  temp = a1.masterordering(1:750*wh(p));  % ordered indices into 10k
  ntrials(p) = length(temp);
  temp2 = union(temp(:),[]);              % unique list of indices relative to 10k
  commonset = union(commonset,find(ismember(1:1000,temp)));  % aggregated across subjects, what is in the shared 1000?
  totaln = totaln + length(setdiff(temp2,1:1000));  % unique images from each subject, count them!
  totaltrials = totaltrials + length(temp);         % just count trials
  unum(p) = length(temp2);                % how many unique images are shown at least once to the subject
  hnum(p) = sum(ismember(1:1000,temp));   % how many of the shared 1,000 were shown at least once to the subject
  for q=1:10000    % loop over 10k
    rec(p,q) = sum(temp==q);  % for each subject and each image, how many trials were experienced?
  end
  allshownims = union(allshownims,a1.subjectim(p,temp2));
end
totaln = totaln + length(commonset);  % add in the aggregated shared 1000 images

% >> sum(rec>=1,2)'
% 
% ans =
% 
%   Columns 1 through 6
% 
%        10000       10000        9411        9209       10000        9411
% 
%   Columns 7 through 8
% 
%        10000        9209
% 
% 
% >> sum(rec>=2,2)'
% 
% ans =
% 
%   Columns 1 through 6
% 
%        10000       10000        8355        7846       10000        8355
% 
%   Columns 7 through 8
% 
%        10000        7846
% 
% >> sum(rec>=3,2)'
% 
% ans =
% 
%   Columns 1 through 6
% 
%        10000       10000        6234        5445       10000        6234
% 
%   Columns 7 through 8
% 
%        10000        5445

% how many unique images are shown at least once to the subject
unum
%10000       10000        9411        9209       10000        9411       10000        9209

% how many trials each subject participated in?
ntrials
%       30000       30000       24000       22500       30000       24000       30000       22500

% how many of the shared 1,000 were shown at least once to the subject
hnum
%        1000        1000         930         907        1000         930        1000         907

% indices into the shared 1,000 such that all 8 subjects got all 3 trials
ix = find(all(rec(:,1:1000)==3,1));
special515 = a1.sharedix(ix)';     % this is the special 515 (73k IDs)

% completed nsd sessions
% 40, 40, 32, 30, 40, 32, 40, 30

% rest in 21-30
% subj01 and subj05 have rest in 31-38

% from the shared 1000 how many reps?
sum(all(rec(:,1:1000)==3,1))  % 515 were shown all 3 times to all subjects
sum(all(rec(:,1:1000)>=2,1))  % 766 were shown at least 2 times to all subjects
sum(all(rec(:,1:1000)>=1,1))  % 907 were shown at least once to all subjects

% total number of distinct images in the dataset (aggregated)
totaln
% 70566

% total number of trials
totaltrials
% 213000

% write out list of images in the 73k that were NOT shown
fid = fopen('notshown.tsv','w');
fprintf(fid,'%d\n',setdiff(1:73000,allshownims));
fclose(fid);

%%%%%%%%%%%%%%%%

%43, 43, 35, 33, 43, 35, 43, and 33 

% prffloc
8*((300*6 + 312*6)/60/60)
8.16 hours

% nsd
(40+40+32+30+40+32+40+30)*1
284 hours

% resting-state
6*(100/60) + 2*(180/60)
16 hours

% syn
8*((268*1.6)/60/60*8)
7.62311111111111 hours

% imagery
8*((150*1.6)*9 + (300*1.6)*3)/60/60
8 hours

% total (not resting)
8.16+284+7.62311111111111+8

total:
16 Hours of resting state
307.783111111111 Hours of task.

on average each subject got:
2 hours of resting-state fMRI
38.4728888888889 hours of task fMRI => 38.47 hours
