% load
file0 = sprintf('~/nsd/nsddata/information/b3pcnum_%s.tsv','func1pt8mm');
a1 = load(file0);

% compute
tempmd = [];
temppp = [];
for subjix=1:8
  tempmd(subjix) = nanmean(a1(:,subjix));
  temppp(:,subjix) = prctile(bootstrapdim(a1(:,subjix),1,@(x) nanmean(x,1),1000),[16 84]);
end

% visualize
figureprep([100 100 185 110],1); hold on;
cmap0 = jet(8);
cmap0(3:5,:) = (cmap0(3:5,:)*5 + [0 0 0])/6;
alldata = [];
alllabels = {};
for p=1:8
  temp = filterout(a1(:,p),NaN,0)';
  alldata = cat(1,alldata,temp);
  temp2 = repmat({num2str(p)},[length(temp) 1]);
  alllabels = [alllabels; temp2];
end
set(gca,'ColorOrder',upsamplematrix(cmap0,[5 1],[],[],'nearest'));
opt = {'Width',.4,'Bandwidth',1,'ViolinAlpha',0.1,'ShowData',true};  %'BoxWidth',.1,'BoxColor',[0 0 1],'MedianColor',[0 1 0],'ShowNotches',false,'ShowMean',false   
%'ViolinColor',cmap0,   'EdgeColor',(cmap0(p,:)+[0 0 0])/2,
violins = violinplot(alldata,alllabels,opt{:}); 
for p=1:8
  delete(violins(p).WhiskerPlot);
  delete(violins(p).BoxPlot);  % this is the IQR box. delete it!
  delete(violins(p).MedianPlot);  % delete the median dot.
  set(violins(p).ScatterPlot,'SizeData',4^2);
  set(violins(p).ViolinPlot,'EdgeColor',(cmap0(p,:)+[.5 .5 .5])/2);
end
for p=1:8
%  h = bar(p,tempmd(p));
%  set(h,'FaceColor',cmap0(p,:));
%  errorbar2(p,tempmd(p),temppp(:,p),'v','k-','LineWidth',2);
  h = errorbar2(p,tempmd(p),temppp(:,p),'v','k-','LineWidth',2);
  set(h,'Color',cmap0(p,:));
%  h2 = scatter(p,tempmd(p),'rx');%,'filled');
%  set(h2,'CData',cmap0(p,:));
  plot([p-.4 p+.4],repmat(tempmd(p),[1 2]),'k-','Color',cmap0(p,:),'LineWidth',2);
end
ax = axis;
axis([0 9 ax(3:4)]);
%axis equal;
set(gca,'XTick',1:8);
set(gca,'YTick',0:2:10);
xlabel('Subject');
ylabel(sprintf('Number of\nGLMdenoise\nregressors'));
figurewrite('pcnum',[],-2,'~/Dropbox/KKTEMP/');

% example fracvalue
% cp ~/nsddata_betas/ppdata/subj05/func1mm/betas_fithrf_GLMdenoise_RR/FRACvalue_session10.nii.gz  ~/Dropbox/KKTEMP/
% 42 81 55
% 0.7 to 1. copper.
% 0.85 to 1. copperupperhalf

%%%%% JUNK:

% NSD258, nsd10

zz = 103;
glmdir = regexprep(datadirs{zz},'rawdata','glmdata');
a2 = load(sprintf('%s/GLMdenoise_nsdBASICsingletrialfithrfGLMdenoiseFracridge.mat',glmdir), ...
  'pcnum','totalbadness','pcregressors');
a8 = load(sprintf('%s/GLMdenoise_nsdBASIConoff.mat',glmdir),'R2');
totalbadness = squish(a2.totalbadness,3);
xvaltrend = -median(totalbadness(a8.R2(:)>5,:),1);  % NOTE: sign flip so that high is good

figure;plot(a2.pcregressors{1}(:,1:4),'-')
figure;plot(xvaltrend);

% load and visualize the pc curves / pcnum
