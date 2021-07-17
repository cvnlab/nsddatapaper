%% %%% SETUP

% load
rdata = {};
for p=1:8
  file0 = sprintf('~/nsd/nsddata/ppdata/subj%02d/behav/responses.tsv',p);
  rdata{p} = importdata(file0);
end

% define
prffloc = [
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
];

nsdsyn = [
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
];

nsdimag = [
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
];

nsd01 = [
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
DATECENSORED
];

nsess = [40 40 32 30 40 32 40 30];

%% %%% FIGURE

figureprep([100 100 1300 150],1); hold on;
cmap0 = jet(8);
cmap0(3:5,:) = (cmap0(3:5,:)*5 + [0 0 0])/6;
dotsize = 25;
cfun = @(t1) datenum(str2double(t1(1:4)),str2double(t1(5:6)),str2double(t1(7:8)));
for subjix=1:8

  % deal with nsd sessions
  for sess=1:40
    ix = find(rdata{subjix}.data(:,2)==sess);
    if isempty(ix)
      continue;
    end
    temp = unique(floor(rdata{subjix}.data(ix,7)));  % midnight of the day
    if (sess>=21 && sess<=30) || (ismember(subjix,[1 5]) && (sess >= 31 && sess <= 38))  % resting sessions
      h = scatter(flatten(temp),subjix*ones(1,length(temp)),dotsize,'rx');%,'filled');
      set(h,'MarkerEdgeColor',cmap0(subjix,:),'LineWidth',2);
    else
      h = scatter(flatten(temp),subjix*ones(1,length(temp)),dotsize,'ro','filled');
      set(h,'CData',cmap0(subjix,:));
    end
  end

  % special sessions are relative to nsd01 session!
  set(scatter(cfun(num2str(prffloc(subjix))) - cfun(num2str(nsd01(subjix))),subjix,dotsize,'r+'),'MarkerEdgeColor',cmap0(subjix,:),'LineWidth',2);
  set(scatter(cfun(num2str(nsdsyn(subjix)))  - cfun(num2str(nsd01(subjix))),subjix,dotsize,'r+'),'MarkerEdgeColor',cmap0(subjix,:),'LineWidth',2);
  set(scatter(cfun(num2str(nsdimag(subjix))) - cfun(num2str(nsd01(subjix))),subjix,dotsize,'r+'),'MarkerEdgeColor',cmap0(subjix,:),'LineWidth',2);

  % total fMRI sessions
  text(330,subjix,sprintf('1+%d+2 = %d',nsess(subjix),nsess(subjix)+3));

end
axis([-100 400 -1 10])
set(gca,'Color',[1 1 1]);
xlabel('Days');
ylabel('Subject number');
set(gca,'YDir','reverse');
set(gca,'YTick',1:8);
figurewrite('timeline',[],-1,'~/Desktop');

%%%%%%%%%%%%%%%%%

caption:
- small disclaimer about split sessions?
