function nsd_et_preproc_figures(data, plotset)
% ----------------------------------------------------------------------- %
%  Diagnostic figures of the cleaned NSD eye-tracking data. MN, July 2021
% ----------------------------------------------------------------------- %

if strcmp(plotset.flag, 'pupil')
    
    % plot pupil area over time
    figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w');
    subplot(2,1,1); plot(data{plotset.cRun}.samples(:,4), 'Color', [0 .55 0.7], 'LineWidth', 1.5);
    h = title(sprintf('Raw: participant %d - run %d', plotset.cSub, plotset.cRun));
    set(gca, 'XTick', []); set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
    xlim([0, numel(data{plotset.cRun}.samples(:,4))]);  box off; ylabel('Pupil area', 'FontSize', 12, 'FontWeight', 'bold');
    set(h, 'FontName', 'Gill Sans MT', 'FontSize', 15, 'FontWeight', 'bold');
    subplot(2,1,2); plot(data{plotset.cRun}.samples_clean(:,4), 'Color', [0 .55 0.7], 'LineWidth', 1.5);
    hold on; plot(data{plotset.cRun}.samples_blinks(:,4), 'Color', [0.8 .4 0], 'LineWidth', 1.5);
    xlabel('seconds', 'FontName', 'Gill Sans MT', 'FontSize', 12, 'FontWeight', 'bold');
    box off; h = title('Preprocessed'); set(h, 'FontName', 'Gill Sans MT', 'FontSize', 15, 'FontWeight', 'bold');
    xlim([0, numel(data{plotset.cRun}.samples_clean(:,4))]); ylabel('Pupil area', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'XTick', [1 numel(data{plotset.cRun}.samples_clean(:,4))/2, numel(data{plotset.cRun}.samples_clean(:,4))])
    set(gca, 'XTickLabel', [0 round((numel(data{plotset.cRun}.samples_clean(:,4))/plotset.newFR)/2) round((numel(data{plotset.cRun}.samples_clean(:,4))/plotset.newFR))])
    
    % save fig
    if exist(plotset.output)~=7; mkdir(plotset.output); end
    outname = strcat('pupil_', plotset.subname, data{plotset.cRun}.filename(8:end-12), '.jpg');
    saveas(gcf,fullfile(plotset.output, outname{1}));
    
elseif strcmp (plotset.flag, 'gazeposition')
    
    % plot 2D scatter plot of gaze position
    figure('units','normalized','outerposition',[0 0 1 0.75]); set(gcf,'color','w');
    subplot(2,5,1); plot(data{plotset.cRun}.samples_clean(:,2), data{plotset.cRun}.samples_clean(:,3), '.','Color', [0 .55 0.7]);
    ylabel('Y°', 'FontSize', 12, 'FontWeight', 'bold'); xlabel('X°', 'FontSize', 12, 'FontWeight', 'bold');
    xlim([-5 5]); ylim([-5 5]);
    
    % plot horizontal gaze position over time
    subplot(2,5,2:5); plot(data{plotset.cRun}.samples_clean(:,2), 'Color', [0 .55 0.7], 'LineWidth', 1.5);
    h = title(sprintf('Gaze postion: participant %d - run %d', plotset.cSub, plotset.cRun));
    set(gca, 'XTick', []); set(gca, 'FontName', 'Gill Sans MT', 'FontSize', 12);
    xlim([0, numel(data{plotset.cRun}.samples_clean(:,2))]); ylim([-5 5]);
    box off; ylabel('X°', 'FontSize', 12, 'FontWeight', 'bold');
    set(h, 'FontName', 'Gill Sans MT', 'FontSize', 15, 'FontWeight', 'bold');
    tmp = data{plotset.cRun}.samples_blinks(:,3); tmp(~isnan(tmp)) = 0;
    hold on; plot(1:numel(tmp), double(tmp), 'Marker','.', 'Color', [0.8 .4 0],'MarkerSize', 15);
    
    % plot Euclidean deviation from gaze center
    subplot(2,5,6);
    hist(data{plotset.cRun}.euclDist(data{plotset.cRun}.euclDist<10), 50);
    h = findobj(gca,'Type','patch'); h.FaceColor = [0 .55 0.7];
    ylabel('n samples', 'FontSize', 12, 'FontWeight', 'bold'); xlabel('EuclDist from median', 'FontSize', 12, 'FontWeight', 'bold')
    
    % plot vertical gaze position over time
    subplot(2,5,7:10); plot(data{plotset.cRun}.samples_clean(:,3), 'Color', [0 .55 0.7], 'LineWidth', 1.5);
    xlabel('seconds', 'FontName', 'Gill Sans MT', 'FontSize', 12, 'FontWeight', 'bold');
    box off; xlim([0, numel(data{plotset.cRun}.samples_clean(:,3))]);  ylim([-5 5]); ylabel('Y°', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'XTick', [1 numel(data{plotset.cRun}.samples_clean(:,3))/2, numel(data{plotset.cRun}.samples_clean(:,3))])
    set(gca, 'XTickLabel', [0 round((numel(data{plotset.cRun}.samples_clean(:,3))/plotset.newFR)/2) round((numel(data{plotset.cRun}.samples_clean(:,3))/plotset.newFR))])
    tmp = data{plotset.cRun}.samples_blinks(:,3); tmp(~isnan(tmp)) = 0;
    hold on; plot(1:numel(tmp), double(tmp), 'Marker','.', 'Color', [0.8 .4 0],'MarkerSize', 15);
    
    % save fig
    if exist(plotset.output)~=7; mkdir(plotset.output); end
    outname = strcat('XY_', plotset.subname, data{plotset.cRun}.filename(8:end-12), '.jpg');
    saveas(gcf,fullfile(plotset.output, outname{1}));end
close all
end