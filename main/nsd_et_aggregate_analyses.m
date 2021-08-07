% ----------------------------------------------------------------------- %
% Aggregate and group-level analyses of the NSD eye-tracking data.
% M.N. July, 2021
% ----------------------------------------------------------------------- %

% settings
clearvars;
paths2ET    = '/System/Volumes/Data/misc/data61/naum2/NSD/analysis_output_final';
outpath     = '/Users/naum2/Desktop/NSD_final_figs';
test_subs   = 1:8;    % which participants do you want to visualize
FR.new      = 100;    % frame rate of eye-tracking data after preprocessing (Hz)
cutoff      = 33.333; % in percent, exclude runs with fewer valid samples
subj_colors = jet(8); % colors for individual subjects
subj_colors(3:5,:) = (subj_colors(3:5,:)*5 + [0 0 0])/6;
line_width  = 2;

% ET files
files_ET = dir(fullfile(paths2ET,'preproc*'));
files_ET    = dir(fullfile(paths2ET,'preproc*')); files_ET = fullfile(files_ET(1).folder, {files_ET(:).name}');
data_ET     = arrayfun(@(y) load(files_ET{y}), test_subs, 'uni', 0);
data_ET     = arrayfun(@(y) data_ET{y}.data(~cellfun(@isempty, data_ET{y}.data)),1:numel(data_ET), 'uni', 0);
valid_runs  = arrayfun(@(y) cell2mat(arrayfun(@(x) data_ET{y}{x}.valid_ratio>cutoff, 1:numel(data_ET{y}), 'uni', 0))', 1:numel(data_ET), 'uni', 0);
data_ET     = arrayfun(@(y) data_ET{y}(valid_runs{y}), 1:numel(data_ET), 'uni', 0);

% ---------------------------------------------------------------------- %
%                       Effect of preprocessing
% ---------------------------------------------------------------------- %
% plot examplary subject and scanning run before & after preprocessing
cSub     = 1;
cRun     = 5;
zoom_win = [1:900]+5600;
FR.plot  = 50;

% downsample again for plotting
data_ds    = downsample(data_ET{cSub}{cRun}.samples,2000/FR.new);
data_ds_ds = downsample(data_ds,FR.new/FR.plot);
data_clean = downsample(data_ET{cSub}{cRun}.samples_clean,FR.new/FR.plot);

% raw data all
figure('Position',[50,50,1350,350]);
subplot(2,8,1:6); hold all;
plot(data_ds_ds(:,3)-nanmedian(data_ds_ds(:,3)), 'Color', [0.65 0.65 0.65], 'LineWidth', line_width);
h = title('Raw eye-tracking data');
set(gca, 'XTick', []); set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
xlim([1, numel(data_ds_ds(:,2))]); ylim([-2 2]);
box off; ylabel('Y°', 'FontSize', 12, 'FontWeight', 'bold');
set(h, 'FontName', 'Gill Sans MT', 'FontSize', 15, 'FontWeight', 'bold');
set(gca,'linewidth',line_width);
hold off;

% raw data zoom in
subplot(2,8,7:8); hold all;
plot(data_ds(zoom_win,3)-nanmedian(data_ds(:,3)), 'Color', [0.65 0.65 0.65], 'LineWidth', line_width);
h = title('Zoom in on blinks'); set(h, 'FontName', 'Gill Sans MT', 'FontSize', 15, 'FontWeight', 'bold');
set(gca, 'XTick', []); set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
xlim([1, numel(zoom_win)]); ylim([-2 2]);
set(gca,'linewidth',line_width)
box off; hold off;

% clean data all
subplot(2,8,9:14); hold all;
xticks = [1, ceil(numel(data_clean(:,3))/2), numel(data_clean(:,3))];
plot(data_clean(:,3), 'Color', subj_colors(cSub,:), 'LineWidth', line_width);
h = title('Cleaned eye-tracking data');
set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
xlim([1, xticks(end)]); ylim([-2 2]);
set(gca, 'XTick', xticks); set(gca, 'XTickLabels', {round(xticks./[1, FR.plot, FR.plot])});
ylabel('Y°', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca,'linewidth',line_width);
set(h, 'FontName', 'Gill Sans MT', 'FontSize', 15, 'FontWeight', 'bold');
box off; hold off;

% clean data zoom in
subplot(2,8,15:16); hold all;
xticks = [1, numel(zoom_win)];
plot(data_ET{cSub}{cRun}.samples_clean(zoom_win,3), 'Color', subj_colors(cSub,:), 'LineWidth', line_width);
h = title('Zoom in on blinks'); set(h, 'FontName', 'Gill Sans MT', 'FontSize', 15, 'FontWeight', 'bold');
set(gca, 'XTick', xticks); set(gca, 'XTickLabels', {round(xticks./[1 FR.plot])});
set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([1, numel(zoom_win)]); ylim([-2 2]); set(gca,'linewidth',line_width)
box off; hold off; set(gcf,'color','w');
% NSD_ET_savefig(gcf, sprintf('preproc_sub%d_run%d', cSub, cRun), outpath);


% ---------------------------------------------------------------------- %
%                  Histograms + Gaussian mixture model
% ---------------------------------------------------------------------- %
% plot 2D histograms of gaze positions, a contour indicating 90% of a
% multivariate normal distribution (mnd) fit to the gaze positions, as well
% as a circle containing 90% of the data

% define bins for histogram
bins   = linspace(-3, 3, 81); % min, max, nbins
ggrid      = linspace(bins(1), bins(end), numel(bins)-1);

% threshold for mnd ellipse and circle (defines the contours)
con_level = 0.90;

% for each subject
h1 = figure('Position',[50,50,950,450]);
for cSub = 1:numel(data_ET)
    
    % compute gmdistribution for across-sessions average
    xx = cell2mat(arrayfun(@(x) data_ET{cSub}{x}.samples_clean(:,2), 1:numel(data_ET{cSub}), 'uni', 0)');
    yy = cell2mat(arrayfun(@(x) data_ET{cSub}{x}.samples_clean(:,3), 1:numel(data_ET{cSub}), 'uni', 0)');
    
    % remove NaNs
    nan_id = isnan(sum([xx,yy],2));
    xx = xx(~nan_id); yy = yy(~nan_id);
    
    % compute histogram and gmdist contour
    [ptsxx, ptsyy] = meshgrid(ggrid, ggrid);
    D = hist2d(xx,yy,bins,bins);
    gmfit = gmdistribution.fit([xx(:) yy(:)],1,'Replicates', 10);
    vals = mahal(gmfit,[ptsxx(:) ptsyy(:)]);
    if isvector(vals); vals = vals(:); end
    xsize = size(vals);
    n = sqrt(xsize(1));
    vals = flipud(reshape(vals,[n n xsize(2:end)]));
    
    % compute circle
    EE = cell2mat(arrayfun(@(x) data_ET{cSub}{x}.euclDist, 1:numel(data_ET{cSub}), 'uni', 0)');
    EE = EE(~isnan(EE));
    thr = 0:0.001:6;
    for n = 1:numel(thr); EEthr(n) = sum(EE<thr(n))./numel(EE); end
    [~, idx] = min(abs(EEthr - con_level));
    u = zeros(size(ptsxx));
    u((ptsxx.^2+ptsyy.^2)<thr(idx)^2)=1;
    
    % plot histogram
    figure(h1)
    subplot(2,4,cSub); hold all;
    imagesc(bins, bins, D);
    
    % plot mnd ellipse thresholded at 90%
    [~, hh] = contour(ggrid, ggrid, vals, repmat(-2*log(1-con_level), [1 2]));
    set(hh, 'Color', subj_colors(5,:)); set(hh, 'LineWidth', line_width);
    
    % plot circle containing 90% of data
    [cc, hh] = contour(ggrid,ggrid,u);
    set(hh, 'Color', [0 1 0]); set(hh, 'LineWidth', line_width);
    colormap(hot); set(gcf,'color','w');
    set(gca, 'XTick', [bins(1) 0 bins(end)]);
    set(gca, 'XTickLabels', {sprintf('%.0f°', bins(1)), '0°',sprintf('%.0f°', bins(end))});
    set(gca, 'YTick', [bins(1) 0 bins(end)]);
    set(gca, 'YTickLabels', {sprintf('%.0f°', bins(1)), '0°',sprintf('%.0f°', bins(end))});
    title(sprintf('Participant %d - %d runs', cSub, numel(data_ET{cSub})));
    
    % figure settings
    colormap(hot); set(gcf,'color','w');
    set(gca, 'XTick', [bins(1) 0 bins(end)]);
    set(gca, 'XTickLabels', {sprintf('%.0f°', bins(1)), '0°',sprintf('%.0f°', bins(end))});
    set(gca, 'YTick', [bins(1) 0 bins(end)]);
    set(gca, 'YTickLabels', {sprintf('%.0f°', bins(1)), '0°',sprintf('%.0f°', bins(end))});
    title(sprintf('Participant %d - %d runs', cSub, numel(data_ET{cSub})));
    xlim([bins(1), bins(end)]); ylim([bins(1), bins(end)]);
    hold off;
    
    % % sanity check: scatter plot of XYs
    % figure; plot(xx, yy, '.');
end
% NSD_ET_savefig(h1, 'hist2D_mnd_circ', outpath);

% ---------------------------------------------------------------------- %
%         Fraction of time spent below fixed deviation thresholds
% ---------------------------------------------------------------------- %
% Get euclidean error for subjects and runs
EuclError = arrayfun(@(cSub) cell2mat(arrayfun(@(cRun) ...
    data_ET{cSub}{cRun}.euclDist, ...
    1:numel(data_ET{cSub}), 'uni', 0)'), 1:numel(data_ET), 'uni', 0);

% remove NaNs
EuclError = arrayfun(@(x) EuclError{x}(~isnan(EuclError{x})), 1:8, 'uni', 0);

% compute fraction of time spent below 1deg
t_EuclError = arrayfun(@(x) sum(EuclError{x}<=1)/numel(EuclError{x}), 1:8, 'uni', 1);

% figure settings
figure('Position', [600 100 200 500]);
dot_alpha = 0.8;
dot_sz    = 75;
box_alpha = 0.5;
box_width = .9;

% plot results for a threshold of 1deg
h1 = boxplot(t_EuclError, 'outliersize', .1, 'width', box_width); hold all;
set(gca,'linewidth',line_width);
scatter(1.5 + rand(numel(t_EuclError),1)/1.1, t_EuclError, ...
    dot_sz, subj_colors, 'filled', 'MarkerFaceAlpha', dot_alpha, 'MarkerEdgeColor', 'k', 'LineWidth',.01);
h2 = findobj(gca,'Tag','Box');set(gcf,'color','w');
arrayfun(@(x) patch(get(h2(x),'XData'),get(h2(x),'YData'), [0.65, 0.65, 0.65], 'FaceAlpha',box_alpha), 1:length(h2), 'uni', 0);
h3 = boxplot(t_EuclError,ones(1,numel(t_EuclError)), 'Color', 'k', 'symbol','', 'width', box_width);
set(h3, 'LineWidth', line_width); delete(h1);
h4 = findobj(gca, 'Tag', 'Lower Whisker'); set(h4,'LineStyle','-');
h4 = findobj(gca, 'Tag', 'Upper Whisker'); set(h4,'LineStyle','-');
ylabel('Fraction of time spent deviating less than 1° from fixation', 'FontSize', 12, 'FontWeight', 'bold');
set(gcf,'color','w');  set(gca, 'XTick', [])
set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
ylim([0, max(t_EuclError)+0.1]); xlim([0 2.5]);
box off; hold off;
% NSD_ET_savefig(gcf, 'timewithin1deg_allSubs', outpath);

% recompute results for all other thresholds
thrs = [0.1:0.1:2];
t_EuclError2 = cell2mat(arrayfun(@(y) arrayfun(@(x) sum(EuclError{x}<y)/numel(EuclError{x}), 1:numel(data_ET), 'uni', 1), thrs, 'uni', 0)');

% plot results for all other thresholds
h = figure('Position', [100 100 600 300]); hold all;
arrayfun(@(x) plot(t_EuclError2(:,x), 'LineWidth', 4, 'Color', subj_colors(x,:)), 1:size(t_EuclError2,2), 'uni', 0);
set(gca, 'XTick', thrs(2:2:end).*10);
thr_labels = (strsplit(sprintf(' %.1f°', thrs(2:2:end)), ' '));
set(gca, 'XTickLabels', thr_labels(2:end), 'FontSize', 12, 'FontWeight', 'bold')
xlabel('Threshold', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Fraction of time spent below threshold', 'FontSize', 12, 'FontWeight', 'bold');
box off; set(gca,'linewidth',line_width);
line([numel(thrs)./2,numel(thrs)./2],0:1,'linewidth',line_width, 'LineStyle', '--', 'Color', [0 0 0]);
legend_labels = strsplit(sprintf('  subject %d', 1:8), '  ');
legend(legend_labels(2:end), 'Location', 'best'); legend('boxoff')
hold off;
% NSD_ET_savefig(gcf, 'timebelow_threshold', outpath);

% ---------------------------------------------------------------------- %
%               plot time series of two examplary subjects
% ---------------------------------------------------------------------- %

% settings
plotSubs = [1,5]; % S01 & S05
plotRuns = [7,10]; % scanning runs
tplot = [1:4500]+2000; % select data

for cSub = 1:numel(plotSubs)
    % mk figure
    figure('units','normalized','outerposition',[0 0 0.85 0.5]);
    set(gcf,'color','w');
    
    % plot 2D scatter plot of gaze positions (downsampled)
    subplot(2,8,[1:2, 9:10]);
    plotData = downsample(data_ET{plotSubs(cSub)}{plotRuns(cSub)}.samples_clean,FR.new/FR.plot);
    plot(plotData(:,2), plotData(:,3), '.','Color', [0 .55 0.7], 'MarkerSize', 12);
    ylabel('Y°', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('X°', 'FontSize', 12, 'FontWeight', 'bold');
    xlim([-3 3]); ylim([-3 3]);
    set(gca, 'XTick', [-3 0 3]);
    set(gca, 'YTick', [-3 0 3]);
    set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
    set(gca,'linewidth',line_width);
    box off;
    
    % gaze position over time
    subplot(2,8,[3:8]);
    plot(plotData(tplot,3), 'Color', [0 .55 0.7], 'linewidth', line_width);
    h = title('Gaze position');
    set(gca,'linewidth',line_width);
    set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
    set(gca, 'XTick', []);
    set(gca, 'YTick', [-3 0 3]);
    ylabel('Y°', 'FontSize', 12, 'FontWeight', 'bold');
    set(h, 'FontName', 'Gill Sans MT', 'FontSize', 15, 'FontWeight', 'bold');
    ylim([-3 3]);
    box off;
    
    % pupil over time
    subplot(2,8,[11:16]);
    plot(plotData(tplot,4)./max(plotData(:,4)), 'Color', [0 .55 0.7], 'linewidth', line_width);
    h = title('Pupil area');
    set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
    set(gca, 'YTick', [-3 0 3]);
    ylabel('Pupil area', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time', 'FontSize', 12, 'FontWeight', 'bold');
    set(h, 'FontName', 'Gill Sans MT', 'FontSize', 15, 'FontWeight', 'bold');
    
    % plot blinks
    plotData_blinks = downsample(data_ET{plotSubs(cSub)}{plotRuns(cSub)}.samples_blinks,2);
    hold on; plot(plotData_blinks(tplot,4)./max(plotData(:,4)), 'Color', [0.8 .4 0], 'linewidth', line_width);
    ylim([0 1]); set(gca, 'YTick', [0 1]);
    set(gca,'linewidth',line_width);
    box off;
    
    % save figure
    %     NSD_ET_savefig(gcf, sprintf('example_sub%d_run%d', plotSubs(cSub), plotRuns(cSub)), outpath);
end