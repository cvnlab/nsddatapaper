function analyzebehavior_prf_visualize(ppdirs,outputdir)

% load
alldata = {};
for p=1:length(ppdirs)
  alldata{p} = load([ppdirs{p} '/behavioralresults_prf.mat']);
end

% internal constants
expectedtriggers = 188;

% look at tricky trigger stuff
totaln = [];
mdf = [];
realstart = [];
for p=1:length(ppdirs)
  for q=1:length(alldata{p}.badtimes)
  
    % calculate
    totaln(p,q) = length(alldata{p}.badtimes{q})-1;           % total number of real triggers detected
    diff0 = diff(alldata{p}.badtimes{q}(2:end-1));            % empirical difference between middle triggers
    mdf(p,q) = median(diff0);                                 % the typical empirical TR diff
    xxx = round(alldata{p}.badtimes{q}(2:end-1)/mdf(p,q))+1;  % integer volume numbers
    pp = polyfit(xxx,alldata{p}.badtimes{q}(2:end-1),1);      % fit a line
    realstart(p,q) = polyval(pp,1);                           % when did the scanner start?
    
    % important sanity check (dropped triggers are innocuous ones in the middle)
    assert(length(find(diff(alldata{p}.badtimes{q}(1:end-1)) > mdf(p,q)*1.1)) == (expectedtriggers-totaln(p,q)));
    
    %figureprep([100 100 1500 200]); hold on;
    %straightline(alldata{p}.badtimes{q},'v','r-');
    %figurewrite(sprintf('ok_p%03d_q%03d',p,q),[],[],[outputdir '/tests']);

  end
end
figureprep([100 100 1200 400]);
subplot(1,3,1); hold on;
plot(flatten(totaln'),'r.-');
ylabel('totaln');
subplot(1,3,2); hold on;
plot(flatten(mdf')*1000);
ylabel('Median empirical TR diff (ms)');
subplot(1,3,3); hold on;
hist(realstart(:)*1000,20);
xlabel('Estimated start of scanner (ms)');
figurewrite('trigger',[],[],outputdir);

% look at duration and trigger-based duration
figureprep([100 100 1000 300]); hold on;
cnt = 0;
for p=1:length(ppdirs)
  n = length(alldata{p}.totaldur);
  plot(cnt+(1:n),alldata{p}.totaldur,'r-');
  plot(cnt+(1:n),cellfun(@(x) x(end)-x(1),alldata{p}.badtimes),'g-');
  cnt = cnt + n;
end
cnt = 0;
for p=1:length(ppdirs)
  n = length(alldata{p}.totaldur);
  straightline([cnt+.5 cnt+n+.5],'v','k-');
  cnt = cnt + n;
end
set(gca,'XTick',[]);
ylabel('Experiment duration (s)');
title('red is stimulus-based duration; green is trigger-based duration')
figurewrite('duration',[],[],outputdir);

% look at accuracy and RT
figureprep([100 100 1000 700]); hold on;
cnt = 0;
for p=1:length(ppdirs)
  n = size(alldata{p}.results,2);
  plot(cnt+(1:n),alldata{p}.results(4,:),'r-');
  plot(cnt+(1:n),cellfun(@median,alldata{p}.resultsrt)/1000 * 100,'g-');
  cnt = cnt + n;
end
cnt = 0;
for p=1:length(ppdirs)
  n = length(alldata{p}.totaldur);
  straightline([cnt+.5 cnt+n+.5],'v','k-');
  cnt = cnt + n;
end
ax = axis;
axis([ax(1:2) 0 105]);
set(gca,'XTick',[]);
ylabel('Accuracy (%) or RT (cs)');
title('red is accuracy; green is RT');
figurewrite('accuracyrt',[],[],outputdir);

figureprep([100 100 1000 700]); hold on;
accuracy = cellfun(@(x) mean(x.results(4,:)),alldata);  % average across runs, percent correct
bar(accuracy);
ax = axis;
axis([ax(1:2) 60 100]);
ylabel('Accuracy (%)');
figurewrite('final',[],[],outputdir);

% save
save('/home/stone/generic/Dropbox/nsd/analyzebehavior_prf/results.mat','accuracy');

%%%%%%%%%%%%%%%%%%%%%%%%

% ADD?
% asserts on acceptable ranges. green line.??
