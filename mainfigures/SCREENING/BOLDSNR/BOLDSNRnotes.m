% this script looks at prffloc screening...
%   ranks the 14 subjects according to BOLD...
%   and boil it down to one number per subject.

% visual inspection of variance explained [make sure they look sane]
cd /Users/kendrick/sessions/nsd/
open */analyzePRF_prf_figures/R2.png
open */GLM*figures/FinalModel.png

% make figure showing histograms
pcts = [];   % dataset x [prf/floc] x thresholds
threshs = 1:1:100;
cmap0 = jet(14);
figureprep([100 100 1000 500]); hold on;
cnt = 1;
for zz=6+(1:14)

  % load
  a1 = load(sprintf('%s/analyzePRF_prf.mat',ppdirs{zz}),'R2');
  a1.R2 = signedarraypower(a1.R2,2)*100;
  a2 = load(sprintf('%s/GLMdenoise_floc.mat',ppdirs{zz}),'R2');  %,'modelmd'

  % plots
  subplot(1,2,1); hold on;
  [nn,xx] = hist(a1.R2(:),0:1:100);
  plot(xx,nn .^ 0.5,'r-','Color',cmap0(cnt,:));
  xlabel('Variance explained (%)');
  ylabel('Sqrt of count');
  title('prf');
  subplot(1,2,2); hold on;
  [nn,xx] = hist(a2.R2(:),0:1:100);
  plot(xx,nn .^ 0.5,'r-','Color',cmap0(cnt,:));
  xlabel('Variance explained (%)');
  ylabel('Sqrt of count');
  title('floc');

  % calculate percentages exceeding thresholds
  for tt=1:length(threshs)
    pcts(cnt,1,tt) = sum(a1.R2(:) > threshs(tt)) / numel(a1.R2) * 100;
    pcts(cnt,2,tt) = sum(a2.R2(:) > threshs(tt)) / numel(a2.R2) * 100;
  end

  % increment
  cnt = cnt + 1;

end
figurewrite('hist',[],-1,'~/Dropbox/KKTEMP');

% make figure showing relative dataset performance
figureprep([100 100 1000 500]); hold on;
subplot(1,2,1); hold on;
plot(unitlength(squish(pcts(:,1,:),2)',2));
xlabel('Threshold (%)');
ylabel('L2-normalized percentages exceeding threshold');
title('prf');
subplot(1,2,2); hold on;
plot(unitlength(squish(pcts(:,2,:),2)',2));
xlabel('Threshold (%)');
ylabel('L2-normalized percentages exceeding threshold');
title('floc');
figurewrite('thresh',[],-1,'~/Dropbox/KKTEMP');

% bar chart
tt = 20;
figureprep([100 100 1000 500]); hold on;
subplot(1,2,1); hold on;
bar(pcts(:,1,tt));
xlabel('Participant');
ylabel('Percentage exceeding threshold');
title('prf');
subplot(1,2,2); hold on;
bar(pcts(:,2,tt));
xlabel('Participant');
ylabel('Percentage exceeding threshold');
title('floc');
figurewrite('rank',[],-1,'~/Dropbox/KKTEMP');

% save
save('~/Dropbox/nsd/figures/BOLDSNR/results.mat','pcts');

%%%%%%%%%%%%%%

notes:
- no compensation for brain size
- on psc: ironically, the more you push over the edge, the more biased it is for PSC to be lower.
  - an alternative strategy would be to simply measure in a fixed anatomical ROI. but.
  - like go to MNI and just define a general
- our simple voxel count is the best bet. bigger brains have an easier time of course...so caveat.
  - this is visually evoked BOLD response (relative to noise), statistical measure.
- but what about beta stability (since denoising may be part of the equation)..
  we have pattern stability. [ehh.. deprecate]
• it's shocking how much subject to subject variability there is.
