%%%%%%%%%% this section generates sample motion traces for all 8 subjects

% define
sessix = 20;
runix = 1;

figureprep([100 100 600 600],1);
for subjix=1:8
  file0 = sprintf('~/nsd/nsddata_timeseries/ppdata/subj%02d/func1pt8mm/motion/motion_session%02d_run%02d.tsv',subjix,sessix,runix);
  a1 = load(file0);  % 226 x 6
  
  subplot(4,2,subjix); hold on;

  yyaxis left;
  plot(a1(:,1),'r-','LineWidth',1,'Color',[1 0 0]);
  plot(a1(:,2),'r-','LineWidth',1,'Color',[.8 0 0])
  plot(a1(:,3),'r-','LineWidth',1,'Color',[.6 0 0]);
  plot(50*a1(:,4),'r-','LineWidth',1,'Color',[0 0 1]);
  plot(50*a1(:,5),'r-','LineWidth',1,'Color',[0 0 .8])
  plot(50*a1(:,6),'r-','LineWidth',1,'Color',[0 0 .6]);
  axis([1 size(a1,1) -2 2]);
  set(gca,'YTick',-4:4);
  if mod(subjix,2)==1
    ylabel('Distance (mm)');
  end

  yyaxis right;
  fd = abs(diff(a1,[],1)) * [1 1 1 50 50 50]';
  plot(1+(1:length(fd)),1+fd,'k-','LineWidth',2);
  straightline(1,'h','k-');
  mdfd = mean(fd);
%  straightline(1+mdfd,'h','c-');
  axis([1 size(a1,1) -2 2]);
  set(gca,'YTick',[1 2],'YTickLabel',{'0' '1'});
  if mod(subjix,2)==0
    ylabel('FD (mm)');
  end
  
  title(sprintf('Subject %d (mean FD = %.2f)',subjix,mdfd));
  if subjix>=7
    xlabel('Volume');
  else
    set(gca,'XTick',[]);
  end

end
figurewrite('motion',[],-1,'~/Dropbox/KKTEMP');

%%%%%%%%%% this section generates some FD summary plots

%%%%% LOAD

% load motion parameters
motiondata = {};  % cell matrix 8 x 40 x 14, each element is TR x 6
for subjix=1:8
  for sessix=1:40
    file0 = sprintf('~/nsd/nsddata_timeseries/ppdata/subj%02d/func1pt8mm/motion/motion_session%02d_run*.tsv',subjix,sessix);
    files = matchfiles(file0);
    for runii=1:length(files)
      a1 = load(files{runii});
      motiondata{subjix,sessix,runii} = a1;
    end
  end
end

% define
allfd   = cell(2,8);  % 2 x 8 cell. each contains all FD observed in all NSD runs and sessions.
                      % second row contains FD calculated ignoring the AP translation parameter.
allfdRS = cell(2,8);  % 2 x 8 cell. each contains all FD observed in all RS runs and sessions
                      % second row contains FD calculated ignoring the AP translation parameter.
weights = {[1 1 1 50 50 50] [0 1 1 50 50 50]};

% do it
for ww=1:length(weights)
  for subjix=1:8
    for sessix=1:40
      for runix=1:14

        % get the motion data
        data0 = motiondata{subjix,sessix,runix};

        % if data were not acquired, just skip
        if ~isempty(data0)

          % compute mean FD in the run (note that one run was discontinuous)
          temp = diff(data0,[],1);
          if subjix==8 && sessix==2 && runix==4
            temp(150,:) = [];  % there was a discontinuity, so we need to ignore
          end
          temp = abs(temp) * weights{ww}';  % volumes x 1

          % for resting-state sessions, we have to adjust the run indexing
          isrest = (sessix>=21 && sessix<=30) || ...
                   (ismember(subjix,[1 5]) && (sessix>=31 && sessix<=38));
          if isrest
            switch runix
            case 1
              allfdRS{ww,subjix} = [allfdRS{ww,subjix} flatten(temp)];
            case 14
              allfdRS{ww,subjix} = [allfdRS{ww,subjix} flatten(temp)];
            otherwise
              allfd{ww,subjix}   = [allfd{ww,subjix}   flatten(temp)];
            end
          else
            allfd{ww,subjix} =     [allfd{ww,subjix}   flatten(temp)];
          end

        end

      end
    end
  end
end

%%%%% PLOT

% define
edges0 = 0:.01:10;
edges0c = (edges0(1:end-1)+edges0(2:end))/2;
thresh = 0.15;  % choose reasonable data-driven threshold
colors0 = {[0 123 200] [223 93 0]};

% plot it
figureprep([100 100 900 185]);
for typ=1:2
  for p=1:8

    subplot(2,8,p); hold on;
    n = histcounts(allfd{typ,p},edges0);
    h = bar(edges0c,n,1);
    set(h,'FaceColor',colors0{typ}/255);
    xlim([0 .5]);
    title(sprintf('%.0f%%',mean(allfd{2,p} > thresh)*100));
    if p==1
      ylabel('Frequency');
    end
    straightline(thresh,'v','r-');
    set(gca,'XTick',0:.2:1);

    subplot(2,8,8+p); hold on;
    n = histcounts(allfdRS{typ,p},edges0);
    h = bar(edges0c,n,1);
    set(h,'FaceColor',colors0{typ}/255);
    xlim([0 .5]);
    title(sprintf('%.0f%%',mean(allfdRS{2,p} > thresh)*100));
    xlabel('FD (mm)');
    if p==1
      ylabel('Frequency');
    end
    straightline(thresh,'v','r-');
    set(gca,'XTick',0:.2:1);

  end
end
figurewrite('fdsummary',[],-2,'~/Dropbox/KKTEMP');
