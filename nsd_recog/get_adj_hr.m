function [ahr, binvals] = get_adj_hr(subjix,binType)
% [ahr, binvals] = get_adj_hr(subjix,binType)
%
% takes in subject id number, loads in responses.tsv file (for a hard-coded
% directory that needs to be specified within function) and returns
% adjusted hit rate values (hit rate minus false alarm rate) based on
% specified binType. binType can be:
%
% 'trial': does not actually bin, but returns trial specific hit value (1
% for hit, 0 for miss) minus the originating sessions FA rate. These values
% can then be binned later (effectively pooling across sessions) or used 
% for other analyses. ahr output is vector with length of the number of 
% relevant trials. binvals is time since last presentation for each trial 
% (vector of same length). Default option.
%
% 'day': for each session, bins trials into time elapsed since last
% presentation as a function of nearest day, rounded down. ahr output is 
% session x bin.
%
% 'logsec': for each session, bins trials into time elapsed since last
% presentation as a function of rounded log(seconds). ahr output is session 
% x bin.
%
% jbh 2/18/21

%% specify any defaults, custom paths, etc
if ~exist('binType','var')||isempty(binType)
    binType = 'trial';
end
% edit path below to get filename from subjix
tsvFileStr = fullfile('tsv_files',sprintf('subj%02d',subjix),'responses.tsv');

%% load in data
beh = tdfread(tsvFileStr);

%% assign key vars
% sess info
nSess = length(unique(beh.SESSION)); 
assert(nSess == max(beh.SESSION), 'Missing Sessions!');

% lag info
timerecent=nan(size(beh.TIME));
for tt = 1:length(timerecent)
    if isnan(beh.MEMORYRECENT(tt))
        continue
    else
        timerecent(tt) = beh.TIME(tt)-beh.TIME(tt-1-beh.MEMORYRECENT(tt));
    end
end
assert(all(~isnan(timerecent)==(beh.ISOLD==1)),'Old item/lag mismatch!');

%% grab sess fa
sess_fa = arrayfun(@(x)(get_sess_fa(x,beh)),1:nSess);

%% assign bin values (go home if 'binning' by trial)
switch binType
    case 'trial'
        oldInds = beh.ISOLD == 1;
        binvals = timerecent(oldInds);
        tsfa = sess_fa(beh.SESSION(oldInds))';
        tsh = beh.ISCORRECT(oldInds)==1;
        ahr = double(tsh) - tsfa;
        return
    case 'day'
        trialbinvals = floor(timerecent);
    case '10day'
        trialbinvals = ceil(floor(timerecent)./10).*10;
    case 'logsec'
        trialbinvals = round(log(timerecent.*24*60*60));
end
% preallocate
binvals = unique(nonnan(trialbinvals));
ahr = nan(nSess,length(binvals));

%% loop over sessions and bins, getting rates
for ss = 1:nSess
    for bb = 1:length(binvals)
        hitVals = beh.ISCORRECT(and(beh.SESSION==ss,trialbinvals==binvals(bb)));
        if ~isempty(hitVals)
        hr = nanmean(hitVals);
        ahr(ss,bb) = hr - sess_fa(ss);
        end
    end   
end

