function data = nsd_et_import(files2importSamples,data, settings)
% ----------------------------------------------------------------------- %
% import all eye-tracking ascii files. MN, July 2021
% ----------------------------------------------------------------------- %

cRun = settings.cRun; cSub = settings.cSub;

if strcmp(settings.flag, 'samples')
    
    % read files
    data{cRun}.samples = textread(fullfile(files2importSamples(cRun).folder, files2importSamples(cRun).name));
    
elseif strcmp(settings.flag, 'messages')
    
    % read
    tmpImport = textscan(fopen(fullfile(files2importSamples(cRun).folder, files2importSamples(cRun).name)), '%s%s%s%s');
    tmpImport = arrayfun(@(x) tmpImport{x}(1:min(cellfun(@numel, tmpImport))), 1:numel(tmpImport), 'uni', 0);
    tmpImport = horzcat(tmpImport{1:end});
    
    % when did experiment start and end?
    data{cRun}.messages.OnOff = cellfun(@str2num, tmpImport(strcmp(tmpImport(:,3), 'SYNCTIME'), 2));
    
    % all other messages
    data{cRun}.messages.all   = tmpImport; % all messages
    
elseif strcmp(settings.flag, 'events')
    
    % read
    tmpImport = textscan(fopen(fullfile(files2importSamples(cRun).folder, files2importSamples(cRun).name)), '%s');
    tmpImport = arrayfun(@(x) tmpImport{x}(1:min(cellfun(@numel, tmpImport))), 1:numel(tmpImport), 'uni', 0);
    tmpImport = horzcat(tmpImport{1:end});
    
    % saccades and blinks
    saccs.on = cellfun(@str2num, tmpImport(find(strcmp(tmpImport, 'SSACC'))+2)); saccs.off = cellfun(@str2num, tmpImport(find(strcmp(tmpImport, 'ESACC'))+3));
    blinks.on = cellfun(@str2num, tmpImport(find(strcmp(tmpImport, 'SBLINK'))+2)); blinks.off = cellfun(@str2num, tmpImport(find(strcmp(tmpImport, 'EBLINK'))+3));
    
    % match onsets and offsets
    if numel(saccs.on) > numel(saccs.off); saccs.on(end) = []; elseif numel(saccs.on) < numel(saccs.off); saccs.off(1) = []; end
    if numel(blinks.on) > numel(blinks.off); blinks.on(end) = []; elseif numel(blinks.on) < numel(blinks.off); blinks.off(1) = []; end
    data{cRun}.saccs = saccs; data{cRun}.blinks =  blinks;
end

end