% load in master information on the manifold
a1 = load('~/Dropbox/nsddata/templates/hrfs_func1mm.mat');  % hrfs is TR x different-hrfs

% fit parameters
options = optimset('Display','final','FunValCheck','on','MaxFunEvals',Inf,'MaxIter',Inf,'TolFun',1e-6,'TolX',1e-6);
fun = @(pp) pp(7)*subscript(conv(spm_hrf(0.1,[pp(1:6) 50]),ones(30,1)),{linspacefixeddiff(1,10,31)});
params = [];
for p=1:size(a1.hrfs,2)
  if p==18
    seed0 = [4.6 20 0.6 3.2 100 0.1 1.1];  % HACK TO REGULARIZE THE BUMP AWAY
    lb = [0 0 0 0 100 0 0];
    ub = [Inf Inf Inf Inf 100 Inf Inf];
  else
    seed0 = [4.6 20 0.6 3.2 3.1 0.1 1.1];
    lb = [0 0 0 0 0 0 0];
    ub = [Inf Inf Inf Inf Inf Inf Inf];
  end
  [params(p,:),d,d,exitflag,output] = lsqnonlin( ...
    @(pp) a1.hrfs(:,p)-fun(pp), ...
    seed0,lb,ub,options);
end

% quick assessment of goodness of fit
figure; hold on;
for p=1:size(a1.hrfs,2)
  subplot(5,5,p); hold on;
  plot(0:30,a1.hrfs(:,p),'r-');
  plot(0:30,fun(params(p,:)),'b-');
end

% another figure (these are taken from Kay et al. Nature Methods 2020)
paramsearly = [7.21 17.6  0.5  4.34 1.82 -3.09   50];  % note we ignore 7th parameter and insert 50 here
paramslate =  [5.76 21.6  1.11 1.72 3.34 0.193   50];

figure; hold on;

temp = spm_hrf(0.1,paramsearly);
h1 = plot(linspacefixeddiff(0,0.1,length(temp)),normalizemax(temp),'r-');

temp = spm_hrf(0.1,paramslate);
h1 = plot(linspacefixeddiff(0,0.1,length(temp)),normalizemax(temp),'b-');

for p=1:size(params,1), p
  temp = spm_hrf(0.1,[params(p,1:6) 50]);
  h1 = plot(linspacefixeddiff(0,0.1,length(temp)),normalizemax(temp),'k-');
  pause;
  set(h1,'Color','g');
end

% save for official records
save('results.mat');

%% EXPORT FOR PUBLIC CONSUMPTION

% save
params(:,7) = 50;
save('hrfparams.mat','params');
