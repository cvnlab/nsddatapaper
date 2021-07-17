function [ahr, binvals] = bincell_to_submat(ahr_cell,binvals_cell)
% [ahr, binvals] = bincell_to_submat(ahr_cell,binvals_cell)
%
% summarizes celled ahr data into subject matrix
%
% jbh 2/26/21

N = length(ahr_cell);

binvals = unique(vertcat(binvals_cell{:}))';
ahr = nan(N,length(binvals));
for ss = 1:N
    [~,binptrs] = ismember(binvals_cell{ss},binvals);
    ahr(ss,binptrs) = nanmean(ahr_cell{ss});
end