% this is a record of what we generated for the participants.
% this is archived but we may do something different...

scores = [1 1 1 1 1 1
1 1 2 1 1 1
2 2 2 2 2 2
1 1 2 1 1 1
3 3 1 1 2 2
1 1 1 1 2 2
3 3 1 1 3 3
3 3 3 3 2 2
3 3 3 3 3 2
1 1 1 1 1 1
1 1 2 1 1 1
1 1 1 1 1 1
3 3 3 2 2 1
1 2 2 1 1 1];

overall = [
1
1.181818182
1.954545455
1.181818182
1.909090909
1.181818182
2.090909091
2.681818182
2.772727273
1
1.181818182
1
2.363636364
1.363636364];

names = { ...
'NSD139'
'NSD149'
'NSD254'
'NSD814'
'NSD244'
'NSD929'
'NSD433'
'NSD841'
'NSD832'
'NSD400'
'NSD258'
'NSD134'
'NSD432'
'NSD168'
};

figureprep([100 100 900 400]); hold on;

subplot(1,4,1); hold on;
barh(mean(scores(:,1:2),2));
set(gca,'YTick',1:length(names),'YTickLabel',names);
set(gca,'YDir','reverse');
title('Brain response (1 is best)');

subplot(1,4,2); hold on;
barh(mean(scores(:,3:4),2));
set(gca,'YTick',1:length(names),'YTickLabel',names);
set(gca,'YDir','reverse');
title('Behavioral performance (1 is best)');

subplot(1,4,3); hold on;
barh(mean(scores(:,5:6),2));
set(gca,'YTick',1:length(names),'YTickLabel',names);
set(gca,'YDir','reverse');
title('Head motion performance (1 is best)');

subplot(1,4,4); hold on;
barh(overall);
set(gca,'YTick',1:length(names),'YTickLabel',names);
set(gca,'YDir','reverse');
title('Overall performance (1 is best)');

figurewrite('participantsummary',[],-1,'~/Dropbox/nsd/figures/PARTICIPANTSUMMARY');
