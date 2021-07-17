function Y = nonnan(X)
% Y = nonnan(X)
% 
% removes 'nan' values.
% 
% jbh 1/6/15

nannies = isnan(X);
Y=X(~nannies);