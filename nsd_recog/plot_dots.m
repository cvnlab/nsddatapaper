function h = plot_dots(xVals,inCell)
% h = plot_dots(xVals,inCell)
% plots each value in inCell as dot jittered artificially in x dim. each
% cell is x coord per xVals
%
% jbh 2/27/21

if nargin == 1
    inCell = xVals;
    xVals = 1:length(inCell);
end

xCell = cell(size(inCell));
xw = min(diff(xVals)).*.7;
for xx = 1:length(inCell)
    xCell{xx} = ((rand(size(inCell{xx}))-.5).*xw)+xVals(xx);
end


X = vertcat(xCell{:});Y = vertcat(inCell{:});

h=scatter(X,Y,[],[0 .1 1],'.');



