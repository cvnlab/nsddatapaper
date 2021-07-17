%% overview
% this script will generate the figure used in the NSD data paper of
% adjusted hit rate as a function of lag between repetitions. 
%
% jbh 2/22/21


%% load in data

% for 10 day bins
[day_ahr, day_binvals, day_subj] = get_allsubs_adj_hr('10day',true);

% for log seconds
[logsec_ahr, logsec_binvals, logsec_subj] = get_allsubs_adj_hr('logsec',true);

N = length(logsec_ahr);

%% generate grp data (removing any bin without all subjects' data)
[grp_day_ahr, day_x] = bincell_to_submat(day_ahr,day_binvals);
bumcols = sum(isnan(grp_day_ahr))>0; grp_day_ahr(:,bumcols) = nan;
[grp_ls_ahr, ls_x] = bincell_to_submat(logsec_ahr,logsec_binvals);
bumcols = sum(isnan(grp_ls_ahr))>0; grp_ls_ahr(:,bumcols) = nan;


%% generate pooled data
pooled_day_ahr = cell(1,length(day_x));
pooled_ls_ahr = cell(1,length(ls_x));
for ss = 1:N
    for bb = 1:length(day_binvals{ss})
        binind = day_x==day_binvals{ss}(bb);
        bindata = nonnan(day_ahr{ss}(:,bb));
        pooled_day_ahr{binind} = vertcat(pooled_day_ahr{binind},bindata);
    end
    for bb = 1:length(logsec_binvals{ss})
        binind = ls_x==logsec_binvals{ss}(bb);
        bindata = nonnan(logsec_ahr{ss}(:,bb));
        pooled_ls_ahr{binind} = vertcat(pooled_ls_ahr{binind},bindata);    
    end
end

%% set up any global vars
chance = 0;
y_limits = [-.6 1];
y_t = min(y_limits):.2:max(y_limits);

%% plot day data
figure;
subplot(1,2,1);
plot_dots(day_x,pooled_day_ahr); % pooled data
hold on
plot_patch(day_x,grp_day_ahr); % grp data
hold on % plot_patch toggles it off
plot(xlim,[chance chance],'k:','LineWidth',1); % chance
hold off

% tweak/gussy
ylim(y_limits);
ylabel('Adjusted hit rate');
xlabel('Time (days)');
xlim([-5 320]);
yticks(y_t);

%% plot logsec data
% figure;
subplot(1,2,2);
plot_dots(ls_x,pooled_ls_ahr); % pooled data
hold on
plot_patch(ls_x,grp_ls_ahr); % grp data
hold on % plot_patch toggles it off
plot(xlim,[chance chance],'k:','LineWidth',1); % chance
hold off

% tweak/gussy
ylim(y_limits);
ylabel('Adjusted hit rate');
[xtdvals, xtdlabels] = get_logsec_xticks;
xlim([min(ls_x)-.6 max(ls_x)+.6]);
xticks(xtdvals);
xticklabels(xtdlabels);
xtickangle(45);
% xlabel('Time');
yticks(y_t);
