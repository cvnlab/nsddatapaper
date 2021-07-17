function analyzebehavior_floc_visualize(ppdirs,outputdir)

% load
alldata = {};
for p=1:length(ppdirs)
  alldata{p} = load([ppdirs{p} '/behavioralresults_floc.mat']);
end

% look at duration
figureprep([100 100 1000 300]); hold on;
cnt = 0;
for p=1:length(ppdirs)
  n = length(alldata{p}.totaldur);
  plot(cnt+(1:n),alldata{p}.totaldur,'r-');
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
figurewrite('duration',[],[],outputdir);

% look at accuracy (hitprop)
figureprep([100 100 1000 700]); hold on;
cnt = 0;
for p=1:length(ppdirs)
  n = length(alldata{p}.hitprop);
  plot(cnt+(1:n),alldata{p}.hitprop * 100,'r-');
  cnt = cnt + n;
end
cnt = 0;
for p=1:length(ppdirs)
  n = length(alldata{p}.hitprop);
  straightline([cnt+.5 cnt+n+.5],'v','k-');
  cnt = cnt + n;
end
ax = axis;
axis([ax(1:2) 0 105]);
set(gca,'XTick',[]);
ylabel('Hit proportion (%)');
figurewrite('accuracy',[],[],outputdir);

% summary
figureprep([100 100 1000 700]); hold on;
accuracy = cellfun(@(x) mean(x.hitprop)*100,alldata);  % average across runs, hit proportion
bar(accuracy);
ax = axis;
axis([ax(1:2) 70 100]);
ylabel('Hit proportion (%)');
figurewrite('final',[],[],outputdir);

% save
save('/home/stone/generic/Dropbox/nsd/analyzebehavior_floc/results.mat','accuracy');

%%%%%%%%%%%%%%%%%%%%%%%%

% ADD?
% asserts on acceptable ranges. green line.??
