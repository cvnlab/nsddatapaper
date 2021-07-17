function [ahr, binvals, subj] = get_allsubs_adj_hr(binType,cellFlag)
% [ahr, binvals, subj] = get_allsubs_adj_hr(binType,cellFlag)
%
% Produces adjusted hit rate values for all NSD subjects as a function of
% different inter-item lag bins specified by binType. if cellFlag is true
% (disabled by default), output is just in cell form with single cell per
% subject with output from get_adj_hr.m. binType can be:
%
% 'trial': does not actually bin, but returns trial specific hit value (1
% for hit, 0 for miss) minus the originating sessions FA rate(effectively 
% pooling across sessions and subjects) ahr output is vector with length of 
% the number of trials. binvals is same shape, but is time since last 
% presentation for each trial. subj is subject ID for each trial.
% Default option.
%
% 'day': based on bins from time elapsed since last presentation as a 
% function of nearest day rounded down, for each session averaged per subject for 
% all sessions with data. ahr output is subjecy x bin. bivals correspond to
% values for each bin. subj is subject id for each row
%
% '10day': values from 'day' are rounded up to nearest 10th day
%
% 'logsec': same as 'day' but bins are based on rounded log(seconds).
%
% jbh 2/22/21

%% specify any high level vars
included_subs = [1 2 3 4 5 6 7 8];
N = length(included_subs);
if ~exist('binType','var')||isempty(binType)
    binType = 'trial';
end

if ~exist('cellFlag','var')||isempty(cellFlag)
    cellFlag = false;
end

%% load data into cell arrays
% using cells for initial loading, then redistribute later, as some info is
% idiosyncratic across subs (bins)
ahr_cell = cell(N,1); binvals_cell = cell(N,1); subj_cell = cell(N,1);
for ss = 1:N
    subjix = included_subs(ss);
   [ahr_cell{ss}, binvals_cell{ss}] = get_adj_hr(subjix,binType);
   subj_cell{ss} = ones(size(ahr_cell{ss},1),1).*subjix;
end

% if you just want celled data, then we're done here
if cellFlag
   ahr = ahr_cell; binvals = binvals_cell; subj = subj_cell;
   return
end

%% summarize based on bin choice
switch binType
    case 'trial'
        ahr = vertcat(ahr_cell{:});
        binvals = vertcat(binvals_cell{:});
        subj = vertcat(subj_cell{:});
    otherwise
        [ahr, binvals] = bincell_to_submat(ahr_cell,binvals_cell);
        subj = included_subs';
end

