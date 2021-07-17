function [h, e] = plot_patch(xVals,data)
% [h, e] = plot_patch(xVals,data)
%
% takes in matrix data (data) where the mean of each column is plotted as a point
% in a graph with SEM error bars. outputs mean data handle, h, and errorbar
% handle, e. if xVals is not specified, defaults to 1:data width
%
% modified from matbar.m by jbh 2/26/21 for nsd code release

% color hard coded
clr = 'k';

if nargin == 1
    data = xVals;
    xVals = 1:size(data,2);
end

Y = nanmean(data);
sem = nanstd(data)./sqrt(sum(~isnan(data)));


% toss column nans?
nannies = isnan(Y);
Y = Y(~nannies);
xVals = xVals(~nannies);
sem = sem(~nannies);


h=plot(xVals,Y,clr,'LineWidth',2);
c=get(h,'Color');
hold on
py = [Y+sem fliplr(Y-sem)];
px = [xVals fliplr(xVals)];
e=patch(px,py,c,'EdgeAlpha',0,'FaceAlpha',.3);
hold off