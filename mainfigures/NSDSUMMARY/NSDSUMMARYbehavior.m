% load behavioral data
bdata = {};
for subjix=1:8
  file0 = sprintf('~/nsd/nsddata/ppdata/subj%02d/behav/responses.tsv',subjix);
  a1 = importdata(file0);
  bdata{subjix} = a1.data;
end
% bdata =
% 
%   1x8 cell array
% 
%   Columns 1 through 5
% 
%     {30000x19 double}    {30000x19 double}    {24000x19 double}    {22500x19 double}    {30000x19 double}
% 
%   Columns 6 through 8
% 
%     {24000x19 double}    {30000x19 double}    {22500x19 double}

% process behavior
allvals = NaN*zeros(8,40,3);
for subjix=1:8
  for sessix=1:40
    
    % if there are none in this session, skip it
    if ~any(bdata{subjix}(:,2)==sessix)
      continue;
    end

    % in this session AND is not missing data
    ix = bdata{subjix}(:,2)==sessix & bdata{subjix}(:,19)==0;
    
%     % ix AND a button was pressed
%     ix2 = ix & isfinite(bdata{subjix}(:,9));
    
    % ix AND is an easy trial (old with respect to current session)
    ix3 =  ix & bdata{subjix}(:,14)==1;
    
    % ix AND is a hard trial (old, but not an easy trial)
    ix3b = ix & bdata{subjix}(:,8)==1 & ~ix3;
    
    % ix AND is novel trial
    ix3c = ix & bdata{subjix}(:,8)==0;
    
    % calculate the three quantities of interest
    allvals(subjix,sessix,1) = sum(bdata{subjix}(ix3,18)==2) / sum(ix3) * 100;    % easy and got it right
    allvals(subjix,sessix,2) = sum(bdata{subjix}(ix3b,18)==2) / sum(ix3b) * 100;  % hard and got it right
    allvals(subjix,sessix,3) = sum(bdata{subjix}(ix3c,18)==2) / sum(ix3c) * 100;  % false alarm
    
  end
end

% visualize
figureprep([100 100 700 110]);
cmap0 = get(gca,'ColorOrder');
for subjix=1:8
  subplot(1,8,subjix); hold on;
  plot(allvals(subjix,:,1),'-','Color',cmap0(1,:));
  plot(allvals(subjix,:,2),'-','Color',cmap0(2,:));
  plot(allvals(subjix,:,3),'-','Color',cmap0(3,:));
  axis([0 40 0 100]);
  xlabel('Scan session');
  title(sprintf('Subject %d',subjix));
  straightline(10:10:40,'v','k:');
  if subjix==1
    ylabel('Response = old (%)');
  else
    set(gca,'YTick',[]);
  end
end
figurewrite('behaviorresults',[],-2,'~/Dropbox/KKTEMP');
