function EEG_analysis(inputDir, outputDir, settings)
%EEG_ANALYSIS Process cleaned EEG data — band power, spectral slope, and ERSP
%
%   EEG_analysis(inputDir, outputDir, settings)
%
%   Inputs:
%       inputDir  - Directory containing cleaned .set files
%       outputDir - Directory for output files
%       settings  - Struct with processing parameters
%
%   Output files (CSV primary, .mat for complex data):
%       bandpower_<roi>.csv   - Per-channel band power (long format)
%       spectral_slope_<roi>.csv - Aperiodic 1/f slope per channel/task
%       ersp_<roi>.mat        - Time-frequency data (too large for CSV)
%       ersp_<roi>_summary.csv - ERSP band summaries per task
%       analysis_log_<timestamp>.csv
%
%   Analyses:
%       1. Band Power    - Traditional frequency band power (Welch)
%       2. Spectral Slope- 1/f aperiodic exponent (FOOOF-style)
%       3. ERSP          - Event-related spectral perturbation

arguments
    inputDir    (1,1) string
    outputDir   (1,1) string
    settings    (1,1) struct = struct()
end

%% ========================================================================
%  Configuration with defaults
%  ========================================================================
defaults = struct();
defaults.task_names     = {'eyes_open','nback_0a','nback_1a',...
                           'nback_0b','nback_2a'};
defaults.task_duration  = 64;

% Frequency bands
defaults.bands = struct(...
    'delta',    [0.5,  4], ...
    'theta',    [4,    8], ...
    'alpha',    [8,   13], ...
    'beta',     [13,  30], ...
    'gamma',    [30,  50]);

% Spectral slope settings (1/f aperiodic)
defaults.slope_range    = [2, 40];      % Frequency range for slope fitting
defaults.slope_exclude  = [8, 13];      % Exclude alpha peak from fit

% ERSP settings
defaults.ersp_epoch     = [-0.5, 1.5];  % Epoch window (s)
defaults.ersp_baseline  = [-0.4, -0.1]; % Baseline window (s)
defaults.ersp_freqs     = [2, 50];      % Frequency range
defaults.ersp_cycles    = [3, 0.5];     % Wavelet cycles [min, factor]
defaults.epoch_thresh   = 100;          % µV threshold for rejection

% ROI definitions (empty = all channels)
defaults.rois = struct();

% Merge with provided settings
settings = mergeStructs(defaults, settings);

%% ========================================================================
%  Initialize
%  ========================================================================
fprintf('\n===== EEG ANALYSIS STARTED =====\n');
fprintf('Analyses: Band Power | Spectral Slope | ERSP\n');
eeglab nogui;

if ~exist(outputDir, 'dir'), mkdir(outputDir); end
fileList = dir(fullfile(inputDir, '**', '*.set'));
if isempty(fileList), error('No .set files found in: %s', inputDir); end
fprintf('Found %d files to process\n', length(fileList));

% ROI setup
roiNames = fieldnames(settings.rois);
outputTypes = ['allchan'; roiNames];

% Initialize accumulators
bp_cells = struct();
slope_cells = struct();
ersp_cells = struct();
for i = 1:length(outputTypes)
    bp_cells.(outputTypes{i}) = {};
    slope_cells.(outputTypes{i}) = {};
    ersp_cells.(outputTypes{i}) = {};
end

logEntries = {};

%% ========================================================================
%  Main processing loop
%  ========================================================================
nFiles = length(fileList);
for fileIdx = 1:nFiles
    filepath = fileList(fileIdx).folder;
    filename = fileList(fileIdx).name;
    
    fprintf('\n[%d/%d] Processing: %s\n', fileIdx, nFiles, filename);
    
    logEntry = struct();
    logEntry.Filename = filename;
    logEntry.SubjectID = '';
    logEntry.Group = '';
    logEntry.nChannels = 0;
    logEntry.Status = '';
    logEntry.Error = '';
    logEntry.Timestamp = char(datetime('now'));
    
    try
        %% Load data
        EEG = pop_loadset('filename', filename, 'filepath', filepath);
        [subjectID, groupName] = extractMetadata(filename, filepath);
        
        logEntry.SubjectID = subjectID;
        logEntry.Group = groupName;
        logEntry.nChannels = EEG.nbchan;
        
        fprintf('  Subject: %s | Group: %s | %d ch | %.1f s\n', ...
                subjectID, groupName, EEG.nbchan, EEG.pnts/EEG.srate);
        
        %% Validate task markers
        [taskMarkerIdx, ~] = findTaskMarkers(EEG, settings.task_names);
        
        if any(taskMarkerIdx == 0)
            missingTasks = settings.task_names(taskMarkerIdx == 0);
            warning('Missing task markers: %s. Skipping.', strjoin(missingTasks, ', '));
            logEntry.Status = 'Skipped';
            logEntry.Error = sprintf('Missing markers: %s', strjoin(missingTasks, ', '));
            logEntries{end+1} = logEntry; %#ok<AGROW>
            continue
        end
        
        %% 1. BAND POWER — per-channel
        fprintf('  → Band power...\n');
        bpAllchan = extractBandPower(EEG, settings, taskMarkerIdx, subjectID, groupName);
        bp_cells.allchan{end+1} = bpAllchan;
        
        for r = 1:length(roiNames)
            roiName = roiNames{r};
            roiChans = settings.rois.(roiName);
            bpROI = extractBandPower(EEG, settings, taskMarkerIdx, ...
                                      subjectID, groupName, roiChans);
            bp_cells.(roiName){end+1} = bpROI;
        end
        
        %% 2. SPECTRAL SLOPE — per-channel
        fprintf('  → Spectral slope (1/f)...\n');
        slopeAllchan = extractSpectralSlope(EEG, settings, taskMarkerIdx, subjectID, groupName);
        slope_cells.allchan{end+1} = slopeAllchan;
        
        for r = 1:length(roiNames)
            roiName = roiNames{r};
            roiChans = settings.rois.(roiName);
            slopeROI = extractSpectralSlope(EEG, settings, taskMarkerIdx, ...
                                             subjectID, groupName, roiChans);
            slope_cells.(roiName){end+1} = slopeROI;
        end
        
        %% 3. ERSP — time-frequency
        fprintf('  → ERSP (time-frequency)...\n');
        erspAllchan = extractERSP(EEG, settings, subjectID, groupName);
        ersp_cells.allchan{end+1} = erspAllchan;
        
        for r = 1:length(roiNames)
            roiName = roiNames{r};
            roiChans = settings.rois.(roiName);
            erspROI = extractERSP(EEG, settings, subjectID, groupName, roiChans);
            ersp_cells.(roiName){end+1} = erspROI;
        end
        
        logEntry.Status = 'Completed';
        logEntries{end+1} = logEntry; %#ok<AGROW>
        fprintf('  ✓ Completed\n');
        
    catch ME
        warning('Error processing %s: %s', filename, ME.message);
        logEntry.Status = 'Failed';
        logEntry.Error = ME.message;
        logEntries{end+1} = logEntry; %#ok<AGROW>
    end
end

%% ========================================================================
%  Save outputs
%  ========================================================================
fprintf('\n→ Writing output files...\n');

nCompleted = sum(cellfun(@(x) strcmp(x.Status, 'Completed'), logEntries));
fprintf('  Processed: %d/%d subjects\n', nCompleted, nFiles);

for i = 1:length(outputTypes)
    outType = outputTypes{i};
    
    % --- Band Power (CSV only) ---
    if ~isempty(bp_cells.(outType)) && ~all(cellfun(@isempty, bp_cells.(outType)))
        bp_table = vertcat(bp_cells.(outType){:});
        bpCsv = fullfile(outputDir, sprintf('bandpower_%s.csv', outType));
        writetable(bp_table, bpCsv);
        fprintf('  ✓ bandpower_%s.csv (%d rows)\n', outType, height(bp_table));
    end
    
    % --- Spectral Slope (CSV only) ---
    if ~isempty(slope_cells.(outType)) && ~all(cellfun(@isempty, slope_cells.(outType)))
        slope_table = vertcat(slope_cells.(outType){:});
        slopeCsv = fullfile(outputDir, sprintf('spectral_slope_%s.csv', outType));
        writetable(slope_table, slopeCsv);
        fprintf('  ✓ spectral_slope_%s.csv (%d rows)\n', outType, height(slope_table));
    end
    
    % --- ERSP (.mat for full TF data, CSV for summary) ---
    if ~isempty(ersp_cells.(outType)) && ~all(cellfun(@isempty, ersp_cells.(outType)))
        erspList = ersp_cells.(outType);
        nonEmpty = ~cellfun(@isempty, erspList);
        
        if any(nonEmpty)
            erspData = vertcat(erspList{nonEmpty});
            
            % Save full TF data to .mat (too large for CSV)
            erspFile = fullfile(outputDir, sprintf('ersp_%s.mat', outType));
            settingsOut = settings; %#ok<NASGU>
            save(erspFile, 'erspData', 'settingsOut', '-v7.3');
            
            % Create summary CSV (mean power per band per time window)
            summaryTable = createERSPSummary(erspData, settings);
            summaryCsv = fullfile(outputDir, sprintf('ersp_%s_summary.csv', outType));
            writetable(summaryTable, summaryCsv);
            
            fprintf('  ✓ ersp_%s.mat + summary.csv (%d subjects)\n', ...
                    outType, length(unique({erspData.SubjectID})));
        end
    end
end

% Write processing log
if ~isempty(logEntries)
    logTable = struct2table([logEntries{:}]);
    logFilename = sprintf('analysis_log_%s.csv', ...
                          datetime('now', 'Format', 'yyyy-MM-dd_HH-mm'));
    writetable(logTable, fullfile(outputDir, logFilename));
    fprintf('  ✓ %s\n', logFilename);
end

fprintf('\n===== EEG ANALYSIS COMPLETE =====\n');
end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function merged = mergeStructs(defaults, overrides)
    merged = defaults;
    if isempty(overrides), return; end
    fields = fieldnames(overrides);
    for i = 1:length(fields)
        if ~isempty(overrides.(fields{i}))
            merged.(fields{i}) = overrides.(fields{i});
        end
    end
end

function [subjectID, groupName] = extractMetadata(filename, filepath)
    tokens = regexp(filename, '(TBI\d{3}|TRE\d{3})', 'match');
    if isempty(tokens)
        tokens = regexp(filepath, '(TBI\d{3}|TRE\d{3})', 'match');
    end
    subjectID = 'UNKNOWN';
    if ~isempty(tokens), subjectID = tokens{1}; end
    
    groupPatterns = {'low_performers', 'normal_performers', 'healthy_controls', ...
                     'impaired', 'maintained', 'healthy'};
    groupName = 'unknown';
    for i = 1:length(groupPatterns)
        if contains(filepath, groupPatterns{i})
            groupName = groupPatterns{i};
            break;
        end
    end
    
    groupMap = containers.Map(...
        {'impaired', 'maintained', 'healthy'}, ...
        {'low_performers', 'normal_performers', 'healthy_controls'});
    if groupMap.isKey(groupName)
        groupName = groupMap(groupName);
    end
end

function [markerIdx, latencies] = findTaskMarkers(EEG, taskNames)
    nTasks = length(taskNames);
    markerIdx = zeros(1, nTasks);
    latencies = zeros(1, nTasks);
    
    if isempty(EEG.event), return; end
    eventTypes = {EEG.event.type};
    
    for t = 1:nTasks
        taskName = taskNames{t};
        idx = find(strcmp(eventTypes, taskName), 1);
        if isempty(idx)
            idx = find(strcmp(eventTypes, [taskName '_stim']), 1);
        end
        if isempty(idx)
            idx = find(strcmp(eventTypes, [taskName '_task']), 1);
        end
        if ~isempty(idx)
            markerIdx(t) = idx;
            latencies(t) = EEG.event(idx).latency;
        end
    end
end

%% ========================================================================
%  BAND POWER — PER-CHANNEL
%  ========================================================================

function bpTable = extractBandPower(EEG, settings, taskMarkerIdx, subjectID, groupName, channelSubset)
arguments
    EEG             (1,1) struct
    settings        (1,1) struct
    taskMarkerIdx   (1,:) double
    subjectID       (1,1) string
    groupName       (1,1) string
    channelSubset   cell = {}
end

[chanIdx, chanLabels] = resolveChannels(EEG, channelSubset);

bandNames = fieldnames(settings.bands);
nBands = length(bandNames);
nChans = length(chanIdx);
nTasks = length(settings.task_names);

nRows = nChans * nTasks;
rowIdx = 0;

absColNames = cellfun(@(b) [b '_abs'], bandNames, 'UniformOutput', false);
relColNames = cellfun(@(b) [b '_rel'], bandNames, 'UniformOutput', false);
logColNames = cellfun(@(b) [b '_log'], bandNames, 'UniformOutput', false);

absData = NaN(nRows, nBands);
relData = NaN(nRows, nBands);
logData = NaN(nRows, nBands);
subjectIDs = repmat({char(subjectID)}, nRows, 1);
groupNames = repmat({char(groupName)}, nRows, 1);
channels = cell(nRows, 1);
tasks = cell(nRows, 1);

for t = 1:nTasks
    if taskMarkerIdx(t) == 0, continue; end
    
    taskName = matlab.lang.makeValidName(settings.task_names{t});
    taskStart = round(EEG.event(taskMarkerIdx(t)).latency);
    taskEnd = min(EEG.pnts, taskStart + round(settings.task_duration * EEG.srate));
    taskSegment = EEG.data(chanIdx, taskStart:taskEnd);
    
    for ch = 1:nChans
        rowIdx = rowIdx + 1;
        
        winLen = min(2 * EEG.srate, size(taskSegment, 2));
        [pxx, freqs] = pwelch(taskSegment(ch, :), winLen, winLen/2, [], EEG.srate);
        totalPower = sum(pxx);
        
        for b = 1:nBands
            freqRange = settings.bands.(bandNames{b});
            freqIdx = (freqs >= freqRange(1)) & (freqs <= freqRange(2));
            absPow = sum(pxx(freqIdx));
            absData(rowIdx, b) = absPow;
            relData(rowIdx, b) = absPow / totalPower;
            logData(rowIdx, b) = 10 * log10(absPow + eps);
        end
        
        channels{rowIdx} = chanLabels{ch};
        tasks{rowIdx} = taskName;
    end
end

validRows = 1:rowIdx;
bpTable = table(subjectIDs(validRows), groupNames(validRows), ...
                channels(validRows), tasks(validRows), ...
                'VariableNames', {'SubjectID', 'Group', 'Channel', 'Task'});

for b = 1:nBands
    bpTable.(absColNames{b}) = absData(validRows, b);
    bpTable.(relColNames{b}) = relData(validRows, b);
    bpTable.(logColNames{b}) = logData(validRows, b);
end
end

%% ========================================================================
%  SPECTRAL SLOPE — 1/f APERIODIC COMPONENT
%  ========================================================================

function slopeTable = extractSpectralSlope(EEG, settings, taskMarkerIdx, subjectID, groupName, channelSubset)
%EXTRACTSPECTRALSLOPE Compute 1/f aperiodic slope (FOOOF-style)
%
%   Fits a line to log-log power spectrum, excluding oscillatory peaks.
%   Returns slope (exponent) and offset per channel per task.

arguments
    EEG             (1,1) struct
    settings        (1,1) struct
    taskMarkerIdx   (1,:) double
    subjectID       (1,1) string
    groupName       (1,1) string
    channelSubset   cell = {}
end

[chanIdx, chanLabels] = resolveChannels(EEG, channelSubset);

nChans = length(chanIdx);
nTasks = length(settings.task_names);
nRows = nChans * nTasks;
rowIdx = 0;

% Pre-allocate
subjectIDs = repmat({char(subjectID)}, nRows, 1);
groupNames = repmat({char(groupName)}, nRows, 1);
channels = cell(nRows, 1);
tasks = cell(nRows, 1);
slopes = NaN(nRows, 1);
offsets = NaN(nRows, 1);
rsquared = NaN(nRows, 1);

slopeRange = settings.slope_range;
excludeRange = settings.slope_exclude;

for t = 1:nTasks
    if taskMarkerIdx(t) == 0, continue; end
    
    taskName = matlab.lang.makeValidName(settings.task_names{t});
    taskStart = round(EEG.event(taskMarkerIdx(t)).latency);
    taskEnd = min(EEG.pnts, taskStart + round(settings.task_duration * EEG.srate));
    taskSegment = EEG.data(chanIdx, taskStart:taskEnd);
    
    for ch = 1:nChans
        rowIdx = rowIdx + 1;
        
        % Compute PSD
        winLen = min(4 * EEG.srate, size(taskSegment, 2));
        [pxx, freqs] = pwelch(taskSegment(ch, :), winLen, winLen/2, [], EEG.srate);
        
        % Select frequency range for fitting
        fitIdx = (freqs >= slopeRange(1)) & (freqs <= slopeRange(2));
        
        % Exclude oscillatory peaks (e.g., alpha)
        if ~isempty(excludeRange)
            excludeIdx = (freqs >= excludeRange(1)) & (freqs <= excludeRange(2));
            fitIdx = fitIdx & ~excludeIdx;
        end
        
        if sum(fitIdx) < 5
            % Not enough points for reliable fit
            channels{rowIdx} = chanLabels{ch};
            tasks{rowIdx} = taskName;
            continue;
        end
        
        % Log-log transform
        logFreqs = log10(freqs(fitIdx));
        logPower = log10(pxx(fitIdx));
        
        % Robust linear fit
        try
            [p, S] = polyfit(logFreqs, logPower, 1);
            slopes(rowIdx) = -p(1);  % Negative because slope is typically negative
            offsets(rowIdx) = p(2);
            
            % R-squared
            yfit = polyval(p, logFreqs);
            SSres = sum((logPower - yfit).^2);
            SStot = sum((logPower - mean(logPower)).^2);
            rsquared(rowIdx) = 1 - (SSres / SStot);
        catch
            % Fit failed
        end
        
        channels{rowIdx} = chanLabels{ch};
        tasks{rowIdx} = taskName;
    end
end

validRows = 1:rowIdx;
slopeTable = table(...
    subjectIDs(validRows), groupNames(validRows), ...
    channels(validRows), tasks(validRows), ...
    slopes(validRows), offsets(validRows), rsquared(validRows), ...
    'VariableNames', {'SubjectID', 'Group', 'Channel', 'Task', ...
                      'Slope', 'Offset', 'R2'});
end

%% ========================================================================
%  ERSP — EVENT-RELATED SPECTRAL PERTURBATION
%  ========================================================================

function erspData = extractERSP(EEG, settings, subjectID, groupName, channelSubset)
%EXTRACTERSP Compute event-related spectral perturbation
%
%   Time-frequency decomposition around stimulus events, baseline normalized.
%   Returns struct array with TF matrices per task.

arguments
    EEG             (1,1) struct
    settings        (1,1) struct
    subjectID       (1,1) string
    groupName       (1,1) string
    channelSubset   cell = {}
end

[chanIdx, ~] = resolveChannels(EEG, channelSubset);

% Get all stimulus events
if isempty(EEG.event)
    erspData = struct([]);
    return;
end

eventTypes = {EEG.event.type};
stimEvents = find(contains(eventTypes, '_stim'));

if isempty(stimEvents)
    erspData = struct([]);
    return;
end

% ERSP parameters
epochWin = settings.ersp_epoch;
baseWin = settings.ersp_baseline;
freqRange = settings.ersp_freqs;
cycles = settings.ersp_cycles;

% Time-frequency parameters
nFreqs = 40;
freqs = linspace(freqRange(1), freqRange(2), nFreqs);

epochSamples = round(epochWin * EEG.srate);
nTimePts = epochSamples(2) - epochSamples(1) + 1;
times = linspace(epochWin(1)*1000, epochWin(2)*1000, nTimePts);

% Baseline indices
baseIdx = (times >= baseWin(1)*1000) & (times <= baseWin(2)*1000);

% Group events by task
taskNames = settings.task_names(~strcmp(settings.task_names, 'eyes_open'));
erspData = struct([]);

for t = 1:length(taskNames)
    taskName = taskNames{t};
    taskStimPattern = [taskName '_stim'];
    taskStimIdx = find(contains(eventTypes, taskStimPattern));
    
    if isempty(taskStimIdx)
        continue;
    end
    
    % Collect valid epochs
    epochs = [];
    for s = 1:length(taskStimIdx)
        evtIdx = taskStimIdx(s);
        lat = round(EEG.event(evtIdx).latency);
        
        epochStart = lat + epochSamples(1);
        epochEnd = lat + epochSamples(2);
        
        if epochStart < 1 || epochEnd > EEG.pnts
            continue;
        end
        
        % Extract and average across ROI channels
        epochData = mean(EEG.data(chanIdx, epochStart:epochEnd), 1);
        
        % Artifact rejection
        if max(abs(epochData)) > settings.epoch_thresh
            continue;
        end
        
        epochs = cat(1, epochs, epochData);
    end
    
    if size(epochs, 1) < 3
        continue;  % Need at least 3 epochs
    end
    
    % Compute time-frequency decomposition via wavelets
    tfMatrix = computeWaveletTF(epochs, EEG.srate, freqs, cycles);
    
    % Convert to power (magnitude squared)
    tfPower = abs(tfMatrix).^2;
    
    % Average across trials
    meanTF = squeeze(mean(tfPower, 1));  % [freqs × time]
    
    % Baseline normalization (dB)
    baseline = mean(meanTF(:, baseIdx), 2);
    erspDB = 10 * log10(bsxfun(@rdivide, meanTF, baseline));
    
    % Store
    idx = length(erspData) + 1;
    erspData(idx).SubjectID = char(subjectID);
    erspData(idx).Group = char(groupName);
    erspData(idx).Task = matlab.lang.makeValidName(taskName);
    erspData(idx).ERSP = erspDB;
    erspData(idx).Freqs = freqs;
    erspData(idx).Times = times;
    erspData(idx).nTrials = size(epochs, 1);
end

if isempty(erspData)
    erspData = struct([]);
end
end

function tfMatrix = computeWaveletTF(epochs, srate, freqs, cycles)
%COMPUTEWAVELETTF Morlet wavelet time-frequency decomposition
%
%   epochs: [nTrials × nTimepoints]
%   Returns: [nTrials × nFreqs × nTimepoints]

[nTrials, nPts] = size(epochs);
nFreqs = length(freqs);

tfMatrix = zeros(nTrials, nFreqs, nPts);

for fi = 1:nFreqs
    freq = freqs(fi);
    
    % Adaptive number of cycles
    if length(cycles) == 2
        nCycles = cycles(1) + (cycles(2) * freq);
        nCycles = min(nCycles, freq);  % Cap at frequency value
    else
        nCycles = cycles(1);
    end
    
    % Create Morlet wavelet
    sigma_t = nCycles / (2 * pi * freq);
    t = -3*sigma_t : 1/srate : 3*sigma_t;
    wavelet = exp(2*1i*pi*freq*t) .* exp(-t.^2 / (2*sigma_t^2));
    wavelet = wavelet / sum(abs(wavelet));  % Normalize
    
    % Convolve each trial
    for tr = 1:nTrials
        convResult = conv(epochs(tr,:), wavelet, 'same');
        tfMatrix(tr, fi, :) = convResult;
    end
end
end

function summaryTable = createERSPSummary(erspData, settings)
%CREATEERSPSUMMARY Create summary table with mean ERSP per band/time window

if isempty(erspData)
    summaryTable = table();
    return;
end

bandNames = fieldnames(settings.bands);
timeWindows = struct('early', [0, 250], 'mid', [250, 500], 'late', [500, 1000]);
windowNames = fieldnames(timeWindows);

rows = {};
for i = 1:length(erspData)
    subj = erspData(i);
    freqs = subj.Freqs;
    times = subj.Times;
    ersp = subj.ERSP;
    
    for b = 1:length(bandNames)
        band = bandNames{b};
        bandRange = settings.bands.(band);
        freqIdx = (freqs >= bandRange(1)) & (freqs <= bandRange(2));
        
        for w = 1:length(windowNames)
            winName = windowNames{w};
            winRange = timeWindows.(winName);
            timeIdx = (times >= winRange(1)) & (times <= winRange(2));
            
            % Mean power in this band/window
            meanPow = mean(ersp(freqIdx, timeIdx), 'all', 'omitnan');
            
            row = {subj.SubjectID, subj.Group, subj.Task, ...
                   band, winName, meanPow, subj.nTrials};
            rows = [rows; row]; %#ok<AGROW>
        end
    end
end

summaryTable = cell2table(rows, 'VariableNames', ...
    {'SubjectID', 'Group', 'Task', 'Band', 'TimeWindow', 'MeanPower_dB', 'nTrials'});
end

%% ========================================================================
%  CHANNEL RESOLUTION HELPER
%  ========================================================================

function [chanIdx, chanLabels] = resolveChannels(EEG, channelSubset)
%RESOLVECHANNELS Get channel indices and labels from subset or all

if isempty(channelSubset)
    chanIdx = 1:EEG.nbchan;
    chanLabels = {EEG.chanlocs.labels};
else
    allLabels = {EEG.chanlocs.labels};
    chanIdx = [];
    chanLabels = {};
    for i = 1:length(channelSubset)
        idx = find(strcmpi(allLabels, channelSubset{i}), 1);
        if ~isempty(idx)
            chanIdx(end+1) = idx; %#ok<AGROW>
            chanLabels{end+1} = allLabels{idx}; %#ok<AGROW>
        end
    end
    if isempty(chanIdx)
        error('No valid channels found in channelSubset');
    end
end
end
