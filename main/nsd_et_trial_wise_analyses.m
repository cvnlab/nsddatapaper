% ----------------------------------------------------------------------- %
% Trial-wise time-resolved analyses of the NSD eye-tracking data.
% M.N. July, 2021
% ----------------------------------------------------------------------- %

% settings
clearvars;
path2logs   = '/System/Volumes/Data/misc/data61/naum2/NSD/data_logs';
paths2ET    = '/System/Volumes/Data/misc/data61/naum2/NSD/analysis_output_final';
outpath     = '/Users/naum2/Desktop/NSD_figures';
frame_rate  = 100;  % frame rate of eye-tracking data after preprocessing (Hz)
trial_dur   = 3;    % trial duration in seconds
cutoff      = 33.333; % in percent, exclude runs with fewer valid samples
subj_colors = jet(8); % colors for individual subjects
subj_colors(3:5,:) = (subj_colors(3:5,:)*5 + [0 0 0])/6;

% subjects
test_subs   = [1:2,4:7];
subs        = dir(fullfile(path2logs, 'sub*')); subs = fullfile(subs(1).folder, {subs(1:end).name}');
subs        = subs(test_subs); % exclude the two synthsubs

% log files
files_logs  = arrayfun(@(x) dir(fullfile(subs{x}, 'func1pt8mm/design', 'design*.tsv')), 1:numel(subs), 'uni', 0);
files_logs  = arrayfun(@(y) arrayfun(@(x) fullfile(files_logs{y}(x).folder, files_logs{y}(x).name), 1:numel(files_logs{y}), 'uni', 0)', 1:numel(subs), 'uni', 0);
trialOns    = arrayfun(@(y) arrayfun(@(x) textread(files_logs{y}{x}), 1:numel(files_logs{y}), 'uni', 0)', 1:numel(subs), 'uni', 0);
trialOns    = arrayfun(@(y) arrayfun(@(x) (find(trialOns{y}{x}>0)-1)./0.999878*(4/3), 1:numel(files_logs{y}), 'uni', 0)', 1:numel(subs), 'uni', 0);

% ET files
files_ET    = dir(fullfile(paths2ET,'preproc*')); files_ET = fullfile(files_ET(1).folder, {files_ET(:).name}');
data_ET     = arrayfun(@(y) load(files_ET{y}), test_subs, 'uni', 0);
data_ET     = arrayfun(@(y) data_ET{y}.data(~cellfun(@isempty, data_ET{y}.data)),1:numel(subs), 'uni', 0);
names_ET    = arrayfun(@(y) arrayfun(@(x) data_ET{y}{x}.filename, 1:numel(data_ET{y}), 'uni', 0)', 1:numel(subs), 'uni', 0);

% exlude runs with fewer than 33.333% valid samples
valid_runs  = arrayfun(@(y) cell2mat(arrayfun(@(x) data_ET{y}{x}.valid_ratio, 1:numel(data_ET{y}), 'uni', 0))', 1:numel(subs), 'uni', 0);
valid_runs  = arrayfun(@(y) valid_runs{y}>cutoff, 1:numel(subs), 'uni', 0);
data_ET     = arrayfun(@(y) data_ET{y}(valid_runs{y}), 1:numel(subs), 'uni', 0);
names_ET    = arrayfun(@(y) names_ET{y}(valid_runs{y}), 1:numel(subs), 'uni', 0);

% find sessions and run names
names_logs = arrayfun(@(y) arrayfun(@(x) strsplit(files_logs{y}{x}, {'/','_','.'}), 1:numel(files_logs{y}), 'uni', 0)', 1:numel(subs), 'uni', 0);
names_logs = arrayfun(@(y) arrayfun(@(x) strcat(names_logs{y}{x}{end-2:end-1}), 1:numel(names_logs{y}), 'uni', 0)', 1:numel(subs), 'uni', 0);
names_ET   = arrayfun(@(y) arrayfun(@(x) strsplit(names_ET{y}{x}, {'_'}), 1:numel(names_ET{y}), 'uni', 0)', 1:numel(subs), 'uni', 0);
names_ET   = arrayfun(@(y) arrayfun(@(x) strcat(names_ET{y}{x}{2:3}), 1:numel(names_ET{y}), 'uni', 0)', 1:numel(subs), 'uni', 0);

% loop over subs and runs
X = cell(1,numel(subs)); Y = X; pupil = X; EE = X;
for cSub = 1:numel(subs)
    for cRun = 1:numel(names_ET{cSub})
        
        % match run names
        run_id = strcmp(names_logs{cSub}, names_ET{cSub}{cRun});
        
        if sum(run_id)
            cTrialOns = trialOns{cSub}{run_id};
            
            % select data of correct run
            cETdata      = data_ET{cSub}{cRun}.samples_clean;
            cETdata(:,5) = data_ET{cSub}{cRun}.euclDist;
            cETdata(:,6) = [1:size(cETdata,1)]./frame_rate;
            
            % split ET data into trials
            frame_id     = arrayfun(@(x) find(cETdata(:,6)>cTrialOns(x) & cETdata(:,6)<(cTrialOns(x)+trial_dur)), 1:numel(cTrialOns), 'uni', 0)';
            trialData    = arrayfun(@(x) cETdata(frame_id{x},:), 1:numel(cTrialOns), 'uni', 0)';
            
            % collect pupil size, X and Y over runs of each sub
            X{cSub}      = [X{cSub}; cell2mat(arrayfun(@(x) trialData{x}(:,2), 1:numel(trialData), 'uni', 0))'];
            Y{cSub}      = [Y{cSub}; cell2mat(arrayfun(@(x) trialData{x}(:,3), 1:numel(trialData), 'uni', 0))'];
            pupil{cSub}  = [pupil{cSub}; cell2mat(arrayfun(@(x) trialData{x}(:,4), 1:numel(trialData), 'uni', 0))'];
            EE{cSub}     = [EE{cSub}; cell2mat(arrayfun(@(x) trialData{x}(:,5), 1:numel(trialData), 'uni', 0))'];
        end
    end
end

% compute heatmaps and Euclidean error for different temporal bins
spat_bins = [-4,4,41]; % min, max, nbins
temp_bins = [0:0.3:3]; % temporal bins

% get correct indizes
temp_bins = round(temp_bins.*frame_rate)+1;
xx = arrayfun(@(y) arrayfun(@(x) X{y}(:,temp_bins(x):temp_bins(x+1)-1)', 1:numel(temp_bins)-1, 'uni', 0), 1:numel(X), 'uni', 0);
xx = arrayfun(@(y) arrayfun(@(x) xx{y}{x}(:), 1:numel(temp_bins)-1, 'uni', 0), 1:numel(xx), 'uni', 0);
yy = arrayfun(@(y) arrayfun(@(x) Y{y}(:,temp_bins(x):temp_bins(x+1)-1)', 1:numel(temp_bins)-1, 'uni', 0), 1:numel(Y), 'uni', 0);
yy = arrayfun(@(y) arrayfun(@(x) yy{y}{x}(:), 1:numel(temp_bins)-1, 'uni', 0), 1:numel(yy), 'uni', 0);

% compute 2D histogram
D = arrayfun(@(y) arrayfun(@(x) hist2d(xx{y}{x},yy{y}{x}, linspace(spat_bins(1),spat_bins(2),spat_bins(3)),linspace(spat_bins(1),spat_bins(2),spat_bins(3))), 1:numel(temp_bins)-1, 'uni', 0), 1:numel(X), 'uni', 0);

% plot single-participant heat maps
figure('Position', [100 100 1500 900]);
n = 0;
for cSub = 1:numel(EE)
    for cBin = 1:numel(temp_bins)-1
        n = n+1;
        subplot(6,numel(temp_bins)-1,n);
        imagesc(log(D{cSub}{cBin}));
        colormap(hot)
        set(gca, 'XTick', [1, ceil(spat_bins(3)/2) spat_bins(3)-1]);
        set(gca, 'XTickLabels', {sprintf('%d°', spat_bins(1)), '0°',sprintf('%d°', spat_bins(2))});
        set(gca, 'YTick', [1, ceil(spat_bins(3)/2) spat_bins(3)-1]);
        set(gca, 'YTickLabels', {sprintf('%d°', spat_bins(2)), '0°',sprintf('%d°', spat_bins(1))});
        title(sprintf('Subj%d - %.1f-%.1f s', cSub, (temp_bins(cBin)-1)./frame_rate, (temp_bins(cBin+1)-1)./frame_rate));
        set(gcf,'color','w');
    end
end
% nsd_et_savefig(gcf, 'trial_wise_heatmap', outpath);

% plot group-level heat maps
figure('Position', [100 100 3500 300]);
for cBin = 1:numel(temp_bins)-1
    subplot(1,numel(temp_bins)-1,cBin);
    tmp = arrayfun(@(y) D{y}{cBin}, 1:numel(D), 'uni', 0);
    imagesc((squeeze(median(cat(3,tmp{:}),3))));
    colormap(hot)
    set(gca, 'XTick', [1, ceil(spat_bins(3)/2) spat_bins(3)-1]);
    set(gca, 'XTickLabels', {sprintf('%d°', spat_bins(1)), '0°',sprintf('%d°', spat_bins(2))});
    set(gca, 'YTick', [1, ceil(spat_bins(3)/2) spat_bins(3)-1]);
    set(gca, 'YTickLabels', {sprintf('%d°', spat_bins(2)), '0°',sprintf('%d°', spat_bins(1))});
    title(sprintf('%.1f-%.1f s', (temp_bins(cBin)-1)./frame_rate, (temp_bins(cBin+1)-1)./frame_rate));
    set(gcf,'color','w');
end

% ---------------------------------------------------------------------- %
%                    Eucl. deviation from median positions
% ---------------------------------------------------------------------- %

% plot settings
figure('Position', [400 400 1500, 1000]);
dot_alpha   = 0.8;
dot_sz      = 75;
box_alpha   = 0.5;
box_width   = .9;
line_width  = 2;

% plot Euclidean error over trials
subplot(2,3,1:3); hold all;
arrayfun(@(y) plot(nanmedian(EE{y}), 'LineWidth', 3, 'Color', subj_colors(test_subs(y),:)), 1:numel(EE), 'uni', 0);
title('Deviation from fixation cross');
set(gca, 'XTick', [0:100:300]);
set(gca, 'XTickLabel', [0:trial_dur]);
set(gca,'linewidth',line_width);
set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
box off; set(gcf,'color','w');
xlabel('seconds since image onset', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Euclidean Error (median across trials)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 15, 'FontWeight', 'bold');
hold off;

% ---------------------------------------------------------------------- %
%                               Pupil size
% ---------------------------------------------------------------------- %
% plot pupil size over time
subplot(2,3,4:6); hold all;
arrayfun(@(y) plot(nanmedian(pupil{y}-nanmean(pupil{y},2)), 'LineWidth', 3, 'Color', subj_colors(test_subs(y),:)), 1:numel(EE), 'uni', 0);
title('Pupil size relative to image onset');
set(gca, 'XTick', [0:100:300]);
set(gca, 'XTickLabel', [0:trial_dur]);
set(gca,'linewidth',line_width);
set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
box off; set(gcf,'color','w');
xlabel('seconds since image onset', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('pupil size (mean-centered)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 15, 'FontWeight', 'bold');
hold off;
% nsd_et_savefig(gcf, 'trial_wise', outpath);

% ---------------------------------------------------------------------- %
%           Proportion of time spent fixating the central -+ 1deg
% ---------------------------------------------------------------------- %

tEE = arrayfun(@(y) EE{y}(~isnan(EE{y})), 1:numel(EE), 'uni', 0);
tEE = arrayfun(@(y) sum(tEE{y}(:)<=1)/numel(tEE{y}), 1:numel(tEE), 'uni', 1);
figure('Position', [400 400 450, 700]);
h1 = boxplot(tEE, 'outliersize', .1, 'width', box_width); hold all;
box off; set(gca,'linewidth',line_width);
scatter(1.5 + rand(numel(tEE),1)/1.1, tEE, ...
    dot_sz, subj_colors(1,:), 'filled', 'MarkerFaceAlpha', dot_alpha, 'MarkerEdgeColor', 'k', 'LineWidth',.01);
ylim([0, max(tEE)+0.1]); xlim([0 2.5]);
h2 = findobj(gca,'Tag','Box');set(gcf,'color','w');
arrayfun(@(x) patch(get(h2(x),'XData'),get(h2(x),'YData'), subj_colors(1,:), 'FaceAlpha',box_alpha), 1:length(h2), 'uni', 0);
h3 = boxplot(tEE,ones(1,numel(tEE)), 'Color', 'k', 'symbol','', 'width', box_width); set(h3, 'LineWidth', line_width); delete(h1);
h4 = findobj(gca, 'Tag', 'Lower Whisker'); set(h4,'LineStyle','-');
h4 = findobj(gca, 'Tag', 'Upper Whisker'); set(h4,'LineStyle','-');
ylabel('Fraction of time spent deviating less than 1° from fixation', 'FontSize', 12, 'FontWeight', 'bold');
set(gcf,'color','w');  set(gca, 'XTick', [])
set(gca,'linewidth',line_width); box off; 
set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
xlim([0 3]);
% nsd_et_savefig(gcf, sprintf('euclerror_allSubs.svg', cSub), outpath);

