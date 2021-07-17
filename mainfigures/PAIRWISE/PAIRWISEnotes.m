% takes files from autoqc_nsd_grand.m...

files0 = matchfiles('~/nsd/ppdata/NSD???-glmdata-glmBASIC/grand.mat');
cmatrix = {};
for p=1:length(files0)
  a0 = load(files0{p});
  cmatrix{p} = a0.cmatrix(2:end,2:end);  % we want only the core NSD
end

% visualize
figureprep([100 100 1200 600]); hold on;
for p=1:8
  subplot(2,4,p);
  imagesc(cmatrix{p},[.9 1]);
  colormap(jet);
%  colorbar;
  axis image tight;
%  xlabel('Session number');
%  ylabel('Session number');
%  set(gca,'XTick',10:10:size(cmatrix{p},2));
%  set(gca,'YTick',10:10:size(cmatrix{p},1));
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  set(gca,'YDir','reverse');
  straightline((5:5:size(cmatrix{p},1)-.1)+.5,'h','w--');
  straightline((5:5:size(cmatrix{p},2)-.1)+.5,'v','w--');
%  title('Pairwise correlation of mean volume');
  title(sprintf('Subject %d',p));
end
figurewrite('pairwise',[],-1,'~/Dropbox/KKTEMP/');

figureprep([100 100 300 300]);
imagesc(randn(10,10),[.9 1]);
colormap(jet);
colorbar;
figurewrite('legend',[],-1,'~/Dropbox/KKTEMP');
