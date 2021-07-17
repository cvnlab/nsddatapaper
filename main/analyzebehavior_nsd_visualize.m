function analyzebehavior_nsd_visualize(ppdirs,outputdir)

% load
alldata = {};
for p=1:length(ppdirs)
  file0 = [ppdirs{p} '/behavioralresults_nsd.mat'];
  alldata{p} = load(file0);
end

% internal constants
expectedtriggers = 188;

% look at tricky trigger stuff
totaln = [];
mdf = [];
realstart = [];
ok = [];
for p=1:length(ppdirs)
  for q=1:length(alldata{p}.badtimes)
  
    % calculate
    ok = [ok alldata{p}.badtimes{q}(1)];
    totaln(p,q) = length(alldata{p}.badtimes{q})-1;           % total number of real triggers detected
    diff0 = diff(alldata{p}.badtimes{q}(2:end-1));            % empirical difference between middle triggers
    mdf(p,q) = median(diff0);                                 % the typical empirical TR diff
    xxx = round(alldata{p}.badtimes{q}(2:end-1)/mdf(p,q))+1;  % integer volume numbers
    pp = polyfit(xxx,alldata{p}.badtimes{q}(2:end-1),1);      % fit a line
    realstart(p,q) = polyval(pp,1);                           % when did the scanner start?
    
    % important sanity check (dropped triggers are innocuous ones in the middle)
    temp = diff(alldata{p}.badtimes{q}(1:end-1));
    temp2 = temp(temp > mdf(p,q)*1.1) / mdf(p,q);
    if any(  abs(temp2-round(temp2))  > .1  )
      emailme(sprintf('analyzebehavior_nsd_visualize1: p=%d (%s), q=%d',p,ppdirs{p},q));
    elseif sum(round(temp2)-1) ~= (expectedtriggers-totaln(p,q))
      emailme(sprintf('analyzebehavior_nsd_visualize2: p=%d (%s), q=%d',p,ppdirs{p},q));
    end
    
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

% check
min(ok)*1000
max(ok)*1000
%          -22.0085950004432
%          -3.14292199618649

% look at duration
figureprep([100 100 2000 300]); hold on;
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
title('red is stimulus-based duration')
figurewrite('duration',[],[],outputdir);

% look at behavioral performance
figureprep([100 100 2000 700]); hold on;
cnt = 0;
for p=1:length(ppdirs)
  n = length(alldata{p}.responserate);
  plot(cnt+(1:n),alldata{p}.responserate,'m-');
  plot(cnt+(1:n),alldata{p}.finalscore,'r-');
  plot(cnt+(1:n),alldata{p}.finalscoreWITHIN,'c-');
  plot(cnt+(1:n),alldata{p}.finalscoreEASY,'k-');
  plot(cnt+(1:n),alldata{p}.medianrt/1000*10,'g-');
  cnt = cnt + n;
end
cnt = 0;
for p=1:length(ppdirs)
  n = length(alldata{p}.responserate);
  straightline([cnt+.5 cnt+n+.5],'v','k-');
  cnt = cnt + n;
end
ax = axis;
axis([ax(1:2) 0 105]);
set(gca,'XTick',[]);
ylabel('Response rate or Percent correct (%) or RT (ds)');
title('magenta is response rate; red/cyan/black is percent correct; green is RT');
figurewrite('accuracyrt',[],[],outputdir);

% but don't we want an overall finalscore stuff?

%%%%%%%%%%%%%%%%%%%%%%%%

% ADD?
% asserts on acceptable ranges. green line.??


% NOTES:
%%% BENH COMMENT:
%  -normalize each pair of behaviorals to sum to 1?
