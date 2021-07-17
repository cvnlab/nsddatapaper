% this script does the final assessment of all of the screening.

%%%%%%%% load all the results!

nsdsetup;

ix = 6+(1:14);

a1 = load('/research/papers/2022 nsd/figures/SCREENING/BOLDSNR/results.mat','pcts');
tt = 20;
boldprf = normalizerange(a1.pcts(:,1,tt)',1,5);
boldfloc = normalizerange(a1.pcts(:,2,tt)',1,5);

a1 = load('~/Dropbox/nsd/analyzebehavior_prf/results.mat','accuracy');
behaviorprf = normalizerange(a1.accuracy,1,5);

a1 = load('~/Dropbox/nsd/analyzebehavior_floc/results.mat','accuracy');
behaviorfloc = normalizerange(a1.accuracy,1,5);

a1 = load('~/Dropbox/nsd/autoqc_ppfmri/results.mat','motionregular','motiondetrend');
motionregular = normalizerange(-a1.motionregular(ix),1,5);
motiondetrend = normalizerange(-a1.motiondetrend(ix),1,5);

%%%%%%%% compute simple summary

overall = (boldprf+boldfloc+behaviorprf+behaviorfloc+motionregular+motiondetrend)/6;

%%%%%%%% plot

% define
toplot = {boldprf boldfloc behaviorprf behaviorfloc motionregular motiondetrend overall};
labels = {'BOLD prf' 'BOLD floc' 'Behavior prf' 'Behavior floc' 'Motion regular' 'Motion detrend' 'Overall'};
names = fsids(6+(1:14));
selectedix = [1 2 4 6 10 11 12 14];
selectedix2 = setdiff(1:14,selectedix);

% do it
figureprep([100 100 1800 250]); hold on;
for p=1:length(toplot)
  subplot(1,length(toplot),p); hold on;
  h1 = bar(selectedix,toplot{p}(selectedix),1);
  set(h1,'FaceColor','y');
  h2 = bar(selectedix2,toplot{p}(selectedix2),1);
  set(h2,'FaceColor',[.5 .5 .5]);
  set(gca,'XTick',1:length(names));%,'XTickLabel',names);
  set(gca,'YTick',1:5);
  axis([0 length(names)+1 0 6]);
  title(labels{p});
end
figurewrite('scores',[],-1,'~/Dropbox/nsd/figures/FINALPARTICIPANTSUMMARY');

%%%%%% subject extrapolation figure

b1 = load('/research/papers/2022 nsd/figures/MAPOFSIGNAL/results.mat');
ix = [1 2 4 6 10 11 12 14];  % 139,149,814,929,400,258,134,168
ix2 = [7 1 2 8 6 5 3 4];
figureprep([0 0 130 130]); hold on;
xx = overall(ix(ix2))';
yy = b1.summary(:,3);
scatter(xx,yy,'ro','filled');
title(sprintf('r=%.2f',calccorrelation(xx,yy)));
pp = polyfit(xx,yy,1);
plot([3.2 5],polyval(pp,[3.2 5]),'k-');
plot([0 3.2],polyval(pp,[0 3.2]),'k:');
xx2 = overall(setdiff(1:14,ix));
yy2 = polyval(pp,xx2);
scatter(xx2,yy2,'bo');
xlabel('Subject ranking');
ylabel('Noise ceiling (b3)');
axis([1 5 -10 60]);
axis square;
set(gca,'XTick',1:5);
figurewrite('extrapolation',[],-2,'~/Dropbox/KKTEMP');
mean(yy)
mean([yy(:); yy2(:)])
%           36.1991070241837
%           25.3725802836003

%%%%%%

notes:
- motivate why tsnr is not a part of this?
