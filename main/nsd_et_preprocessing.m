% ----------------------------------------------------------------------- %
% Preprocessing of the NSD eye tracking data. M.N. July, 2021
% The following preprocessing steps are performed.
%
% EDF to ASCII conversion (requires edf2asci.exe on windows)
% Blink removal
% Temporal alignment to fMRI session
% Pixel to visual-degree conversionå
% Tracking-noise removal
% Detrending & median-centering
% Downsampling
% Smoothing
%
% The preprocessed data are then saved as .mat and .csv files.

% ----------------------------------------------------------------------- %

% Some settings
clearvars
do.plot      = 1; % want figures? yes = 1, no = 0
do.save      = 1; % save output?  yes = 1, no = 0
do.convert   = 0; % edf2asci conversion? yes = 1, no = 0
FR.old       = 2000; % Eyelink frame rate
FR.new       = 100;  % desired frame rate after downsampling
outputFolder = '/System/Volumes/Data/misc/data61/naum2/NSD/analysis_output';

% Find data
path2data = fullfile('/System/Volumes/Data/misc/data61/naum2/NSD/data_raw');
subjnames = dir(fullfile(path2data, 's*')); subjnames = {subjnames.name}';
path2data = fullfile(path2data, subjnames, 'eyedata');
addpath(genpath('/System/Volumes/Data/misc/data61/naum2/NSD/Scripts'))

% ----------------------------------------------------------------------- %
%                           Start preprocessing
% ----------------------------------------------------------------------- %

fprintf('\n Start processing')
for cSub = 7:8 %numel(path2data)
    fprintf('\n    Participant %d', cSub);
    
    % ---- EDF to ASCI conversion & splitting messages from samples ---- %
    Files2convert = dir(fullfile(path2data{cSub}, '*_*.edf'));
    Files2convert = fullfile(path2data{cSub}, {Files2convert(:).name}');
    
    % Convert all .edf files into .ascii format
    if do.convert == 1
        % system('set PATH=%PATH%;X:\Matthias\Scripts\EyeLinkFiles\')
        for cRun  = 1:size(Files2convert,1)
            fprintf('\n Converting file: %s', Files2convert{cRun});
            
            % Select edf file
            cFile = Files2convert{cRun};
            cd('X:\Matthias\Scripts\EyeLinkFiles\');
            
            % ........................................................... %
            % Convert file into .asc giving only sampling data
            system(sprintf('edf2asc.exe -miss NaN -ne -nflags %s', cFile));
            
            % Rename file
            oldName = sprintf('%s.asc', cFile(1:end-4));
            newName = sprintf('%s_samples.asc', cFile(1:end-4));
            movefile(oldName,newName)
            % ........................................................... %
            
            % ........................................................... %
            % Convert file into .asc giving only eye link messages
            system(sprintf('edf2asc.exe -miss NaN -neye -ns -nflags %s', cFile));
            
            % Rename file
            oldName = sprintf('%s.asc', cFile(1:end-4));
            newName = sprintf('%s_messages.asc', cFile(1:end-4));
            movefile(oldName,newName)
            % ........................................................... %
            
            % ........................................................... %
            % Convert file into .asc giving only noise labels (blinks, saccades...)
            system(sprintf('edf2asc.exe -miss NaN -e -nflags %s', cFile));
            
            % Rename file
            oldName = sprintf('%s.asc', cFile(1:end-4));
            newName = sprintf('%s_events.asc', cFile(1:end-4));
            movefile(oldName,newName)
            % ........................................................... %
        end % run loop
    end
    
    % ------------------------------------------------------------------- %
    %                              Data cleaning
    % ------------------------------------------------------------------- %
    data = cell(1,44);
    whichRuns = nsd_et_lookup; whichRuns = whichRuns{cSub}; % Run look-up table
    for cRun  = whichRuns
        
        % ------------------------ data import--------------------------- %
        fprintf('\n       Importing run %d', cRun)
        
        % samples
        settings.flag   = 'samples'; settings.cRun = cRun; settings.cSub = cSub;
        files2import    = dir(fullfile(path2data{cSub}, sprintf('eyedata*_%s*.asc', settings.flag)));
        data            = nsd_et_import(files2import,data,settings);
        data{cRun}.filename = files2import(cRun).name;
        
        % messages
        settings.flag   = 'messages';
        files2import    = dir(fullfile(path2data{cSub}, sprintf('eyedata*_%s*.asc', settings.flag)));
        data            = nsd_et_import(files2import,data,settings);
        
        % Events
        settings.flag   = 'events';
        files2import    = dir(fullfile(path2data{cSub}, sprintf('eyedata*_%s*.asc', settings.flag)));
        data            = nsd_et_import(files2import,data,settings);
        
        % ------------------------ blink removal------------------------- %
        % Note that we consider the eyelink-detected blinks & when the pupil
        % size is zero. We also increase the blink duration to remove the edge.
        fprintf(' - blink removal');
        addms2blink            = 50; % ms added to end of blink
        detected_blinks        = cell2mat(arrayfun(@(x) ...
            (data{cRun}.blinks.on(x)):(data{cRun}.blinks.off(x)+FR.old/1000*addms2blink), ...
            1:numel(data{cRun}.blinks.on), 'uni', 0));
        data{cRun}.blinks.bool = ismember(data{cRun}.samples(:,1),detected_blinks);
        
        % Make sure blinks are detected, even when eyelink detection failed
        data{cRun}.blinks.bool((data{cRun}.samples(:,4) == 0)) = 1;
        
        % adds another 100 ms to the beginning and ending of the blink)
        smth_kernel            = ones(1,FR.old/1000*200)./(FR.old/1000*200);
        data{cRun}.blinks.bool = (conv(data{cRun}.blinks.bool', smth_kernel, 'same')>0)';
        
        % ----------- align eyetracking data with imaging data ---------- %
        fprintf(' - align with fMRI');
        on      = find(data{cRun}.samples(:,1) == data{cRun}.messages.OnOff(1)); on = on(2);
        off     = find(data{cRun}.samples(:,1) == data{cRun}.messages.OnOff(2)); off = off(1);
        data{cRun}.samples     = data{cRun}.samples(on:off,:);
        data{cRun}.blinks.bool = data{cRun}.blinks.bool(on:off,:);
        
        % ------------------- Convert to visual degrees ----------------- %
        data{cRun}.samples(:,2:3) = data{cRun}.samples(:,2:3).*(8.4/714);
        
        % ---------------- split valid samples from blinks -------------- %
        data{cRun}.samples_clean  = nan(size(data{cRun}.samples));
        data{cRun}.samples_blinks = data{cRun}.samples_clean;
        data{cRun}.samples_clean(data{cRun}.blinks.bool~=1,:)  = data{cRun}.samples(data{cRun}.blinks.bool~=1,:);
        data{cRun}.samples_blinks(data{cRun}.blinks.bool==1,:) = data{cRun}.samples(data{cRun}.blinks.bool==1,:);
        
        % -------------------- remove tracking noise -------------------- %
        % Samples deviating more than 6 degrees from median gaze position
        % are removed (+- 250ms).
        fprintf(' - cut tracking noise');
        cutoff      = 6; % deg
        median_gaze = nanmedian(data{cRun}.samples_clean(:,2:3));
        euclDist    = sqrt((data{cRun}.samples_clean(:,2) - median_gaze(1)).^2 ...
            + (data{cRun}.samples_clean(:,3) - median_gaze(2)).^2);
        smthInMs    = 500; smth_kernel = ones(1,FR.old/1000*smthInMs)./(FR.old/1000*smthInMs);
        cut_samples = (conv(euclDist > cutoff', smth_kernel, 'same')>0)';
        data{cRun}.samples_clean(cut_samples,:) = NaN;
        
        % --------------------------- downsampling ---------------------- %
        fprintf(' - downsample');
        data{cRun}.samples_clean  = downsample(data{cRun}.samples_clean, FR.old/FR.new);
        data{cRun}.samples_blinks = downsample(data{cRun}.samples_blinks, FR.old/FR.new);
        %         data{cRun}.samples_raw_ds = downsample(data{cRun}.samples, FR.old/FR.new); % downsampled raw data for plotting
        
        % ------------------ detrending & median centering -------------- %
        fprintf(' - detrend & center');
        data{cRun}.samples_clean(:,2:3) = et_detrend(data{cRun}.samples_clean(:,2:3)')';
        median_center = nanmedian(data{cRun}.samples_clean(:,2:3));
        data{cRun}.samples_clean(:,2)   = data{cRun}.samples_clean(:,2)  - median_center(1);
        data{cRun}.samples_clean(:,3)   = data{cRun}.samples_clean(:,3)  - median_center(2);
        data{cRun}.samples_blinks(:,2)  = data{cRun}.samples_blinks(:,2) - median_center(1);
        data{cRun}.samples_blinks(:,3)  = data{cRun}.samples_blinks(:,3) - median_center(2);
        
        % ------------------------ temporal smoothing ------------------- %
        fprintf(' - smooth');
        smthInMs    = 50; % ms
        smth_kernel = ones(1,FR.new/1000*smthInMs)./(FR.new/1000*smthInMs);
        data{cRun}.samples_clean(:,2:4) = et_smooth(data{cRun}.samples_clean(:,2:4)', smth_kernel)';
        
        % ----------------------- Compute some metrics ------------------ %
        
        % recompute Euclidean distance on smoothed samples
        median_gaze = nanmedian(data{cRun}.samples_clean(:,2:3));
        data{cRun}.euclDist = sqrt((data{cRun}.samples_clean(:,2) - median_gaze(1)).^2 ...
            + (data{cRun}.samples_clean(:,3) - median_gaze(2)).^2);
        
        % how many valid samples remained?
        data{cRun}.valid_ratio = 100*(1-sum(any(isnan(data{cRun}.samples_clean(:,2:3))')')./size(data{cRun}.samples_clean,1));
        
        % -------------------------------- plot ------------------------- %
        if do.plot
            if data{cRun}.valid_ratio>33.333
                fprintf(' - plot');
                plotset.newFR   = FR.new;
                plotset.cSub    = cSub;
                plotset.subname = subjnames(cSub);
                plotset.cRun    = cRun;
                plotset.output  = outputFolder;
                plotset.flag    = 'gazeposition';
                nsd_et_preproc_figures(data, plotset)
                plotset.flag    = 'pupil';
                nsd_et_preproc_figures(data, plotset)
            end
        end
        
        % -------------------------- save csv file ---------------------- %
        if do.save
            if data{cRun}.valid_ratio>33.333
                fprintf(' - Save file');
                if exist(outputFolder)~=7; mkdir(outputFolder); end
                outname = strcat(data{cRun}.filename(1:8), subjnames(cSub), data{cRun}.filename(8:end-3), 'csv');
                writematrix(data{cRun}.samples_clean,fullfile(outputFolder, outname{1}))
            end
        end
    end
    
    % --------------------------- save mat file ------------------------- %
    if do.save
        fprintf(' - Save file');
        save(fullfile(outputFolder, sprintf('preproc_data_subj%d.mat', cSub)), 'data');
    end
    fprintf(' - done.');
    % ------------------------------------------------------------------- %
end

% potentially interesting for some users: split samples into individual TRs
% cID = 1:(1.6*FR.new); cTR = 0; XY = [];
% fprintf(' - TR split');     clearvars XY
% while cID(end)<size(data{cRun}.samples_clean,1)
%     cTR = cTR + 1;
%     XY.samples_ET{cTR}     = data{cRun}.samples_clean(cID, 2:3);
%     XY.samples_blinks{cTR} = ~isnan(data{cRun}.samples_blinks(cID, 4));
%     cID = cID + numel(cID);
% end
% save(fullfile(outputFolder, 'trajectory_ET.mat'), 'XY');