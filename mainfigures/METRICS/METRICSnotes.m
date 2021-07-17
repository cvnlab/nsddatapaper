% this script loads runmetrics.mat and makes a gigantic visualization of them.

% load
a1 = load('~/nsd/nsddata/information/runmetrics.mat');

% define
nsess = [40 40 32 30 40 32 40 30];

% compute mean/median across runs in each session
metrics = NaN*zeros(8,40,6);
for mm=1:6
  if mm==6       % re-compute easy trials
    metrics(:,:,mm) = sum(nanreplace(a1.runmetrics(:,:,:,6))/100 .* a1.runmetrics(:,:,:,7),3) ./ sum(a1.runmetrics(:,:,:,7),3) * 100;
  elseif mm==4   % for response rate, use mean
    metrics(:,:,mm) = nanmean(a1.runmetrics(:,:,:,mm),3);  % there is at least one run where we don't have behavioral data...
  else
    metrics(:,:,mm) = nanmedian(a1.runmetrics(:,:,:,mm),3);
  end
end  

% compute mean/median across sessions
subjvals = zeros(8,6);
subjvalspp = zeros(2,6,8);  % [16,84] x 6 metrics x 8 subjects
for mm=1:6
  for p=1:8
    if mm==4  % response rate uses mean
      subjvals(p,mm) = mean(metrics(p,1:nsess(p),mm),2);
      temp = squish(metrics(p,1:nsess(p),mm),2);  % sess x 1
      subjvalspp(:,mm,p) = prctile(bootstrapdim(temp,1,@(x) mean(x,1),1000),[16 84],1);
    else
      subjvals(p,mm) = median(metrics(p,1:nsess(p),mm),2);
      temp = squish(metrics(p,1:nsess(p),mm),2);  % sess x 1
      subjvalspp(:,mm,p) = prctile(bootstrapdim(temp,1,@(x) median(x,1),1000),[16 84],1);
    end
  end
end

%% VISUALIZE

% plot
figureprep([0 0 1200*.8 .75*1200*.8]);

% calc
cmap0 = jet(8);
cmap1 = jet(8);
cmap1(3:5,:) = (cmap1(3:5,:)*5 + [0 0 0])/6;
markers = {'o' 'x' '+' '*' 's' 'd' 'v' '^'};
markers = repmat({''},[1 8]);
labels0 = {'tSNR' 'FD' 'ON-OFF R^2 (%)' 'Response rate (%)' 'Percent correct (%)' 'Easy trials (%)'};

% OLD VERSION BY RUN
% labels0 = {'tSNR' 'FD' 'R2' 'RR' 'PCTCORRECT' 'EASY'};
% for mm=1:6
%   subplot(3,2,mm); hold on;
%   for p=1:8
%     plot(flatten(squish(a1.runmetrics(p,:,:,mm),2)'),'Color',cmap0(p,:),'LineWidth',1);
%   end
%   temp = squish(nanmedian(a1.runmetrics(:,:,:,mm),1),2);  % 40 x 12
%   plot(flatten(temp'),'k-','LineWidth',2);
%   set(gca,'Color',[.4 .4 .4]);
%   ax = axis;
%   axis([0.5 12*40+.5 ax(3:4)]);
%   straightline((0:12:12*40)+.5,'v','w:');
%   set(gca,'XTick',linspacefixeddiff(13/2,12,40));
%   set(gca,'XTickLabel',1:40);
%   xlabel('Session');
%   ylabel(labels0{mm});
% end

ord = [1 4 2 6 3 5];  % order the panels in this way
for mmii=1:6
  mm = ord(mmii);  % which metric

  subplot(3,2*4,(mmii-1)*4+[1 2 3]); hold on;
%   for ff=1:20
%     for p=1:8
%       whix = picksubset(1:40,[20 ff],p);
%       set(scatter(whix,metrics(p,whix,mm),['r' markers{p}]),'LineWidth',2,'MarkerEdgeColor',cmap1(p,:));
%     end
%   end
  for p=1:8
    plot(metrics(p,:,mm),['r-'],'Color',cmap1(p,:),'LineWidth',2);
  end
  temp = nanmedian(metrics(:,:,mm),1);  % median across subjects
  plot(temp,'k-','LineWidth',4);
%  set(gca,'Color',[.5 .5 .5]);
  ax = axis;
  axis([0.5 40+.5 ax(3:4)]);
  if mm==2
    axis([0.5 40+.5 ax(3) 0.24]);
    set(gca,'YTick',0:0.1:0.3);
    ax = axis;
  end
%  straightline(5:10:40,'v','w:');
%  set(gca,'XTick',linspacefixeddiff(13/2,12,40));
  set(gca,'XTick',5:5:40);
  xlabel('Session');
  ylabel(labels0{mm});
  
  subplot(3,2*4,(mmii-1)*4+4); hold on;
  for p=1:8
    h = bar(p,subjvals(p,mm),1);
    set(h,'FaceColor',cmap0(p,:));
    errorbar2(p,subjvals(p,mm),subjvalspp(:,mm,p),'v','k-','LineWidth',2);
  end
  axis([0 9 ax(3:4)]);
  if mm==2
    axis([0 9 ax(3) 0.24]);
    set(gca,'YTick',0:0.1:0.3);
  end
  set(gca,'XTick',1:8);
  xlabel('Subject');
%  set(gca,'YTick',[]);

end

figurewrite('metrics',[],-1,'~/Dropbox/KKTEMP/');

%% COUNTS

% response rate [maybe subjects fell asleep in a run?] [how many runs were below the threshold?]
temp = filterout(flatten(a1.runmetrics(:,:,:,4)),NaN,0);
thresh = 95;
fid = fopen('~/Dropbox/KKTEMP/results.txt','w');
fprintf(fid,'the thresh is %.3f.\n',thresh);
fprintf(fid,'the number of bad is %d/%d (%.1f%%).\n',sum(temp < thresh),length(temp),sum(temp < thresh)/length(temp)*100);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

notes:
- should we compute FD at 90% percentile? [well...]
- should we show resting state? [well...]
