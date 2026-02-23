function EEG_analysis(inputDir, outputDir, settings)
%EEG_ANALYSIS Process cleaned EEG data — preserve per-channel and per-trial data
%
%   EEG_analysis(inputDir, outputDir, settings)
%
%   Inputs:
%       inputDir  - Directory containing cleaned .set files
%       outputDir - Directory for output .mat files
%       settings  - Struct with processing parameters
%
%   Expected settings fields:
%       settings.task_names   - Cell array of task marker names
%       settings.task_duration- Duration of each task in seconds
%       settings.bands        - Struct of frequency bands
%       settings.erp_epoch    - [start end] in seconds
%       settings.erp_baseline - [start end] in seconds
%       settings.p300_window  - [start end] in seconds
%       settings.epoch_thresh - µV threshold for epoch rejection
%       settings.rois         - Struct of ROI definitions
%
%   Output .mat files (per ROI):
%       bandpower_<roi>.mat  - Per-channel band power (long format)
%       erp_<roi>.mat        - Per-trial ERP metrics + waveforms
%       analysis_log_<timestamp>.csv
%
%   Data structure — bandpower .mat:
%       bp_table: table with one row per Subject × Channel × Task
%           Columns: SubjectID, Group, Channel, Task,
%                    <band>_abs, <band>_rel, <band>_log for each band
%       settings: processing settings used
%
%   Data structure — erp .mat:
%       erp_table: table with one row per Subject × Task × Trial
%           Columns: SubjectID, Group, Task, Trial, Accuracy,
%                    PeakAmp, PeakLat, MeanAmp, CentLat
%       waveform_matrix: [nTotalTrials × nTimepoints] double
%       waveform_info:   table indexing into waveform_matrix
%           Columns: SubjectID, Group, Task, Trial, Accuracy
%       times: [1 × nTimepoints] time vector in ms
%       settings: processing settings used

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
    'gamma',    [30,  50], ...
    'broadband',[0.5, 55]);

% ERP settings
defaults.erp_epoch      = [-0.2, 0.8];
defaults.erp_baseline   = [-0.2, 0];
defaults.p300_window    = [0.25, 0.5];
defaults.epoch_thresh   = 100;

% ROI definitions
defaults.rois = struct();

% Merge with provided settings
settings = mergeStructs(defaults, settings);

%% ========================================================================
%  Initialize
%  ========================================================================
fprintf('\n===== EEG ANALYSIS STARTED =====\n');
eeglab nogui;

%-- Initial file validation
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
fileList = dir(fullfile(inputDir, '**', '*.set'));
if isempty(fileList), error('No .set files found in: %s', inputDir); end
fprintf('Found %d files to process\n', length(fileList));

%-- Initial ROI validation
roiNames = fieldnames(settings.rois);
outputTypes = ['allchan'; roiNames];

%-- Initialize cell accumulators (one cell per subject, vertcat at end)
bp_cells = struct();
erp_cells = struct();
wf_cells = struct();   % waveform data
for i = 1:length(outputTypes)
    bp_cells.(outputTypes{i}) = {};
    erp_cells.(outputTypes{i}) = {};
    wf_cells.(outputTypes{i}) = {};
end

%-- Initialize processing log
logEntries = {};

%% ========================================================================
%  Main processing loop
%  ========================================================================
nFiles = length(fileList);
for fileIdx = 1:nFiles
    filepath = fileList(fileIdx).folder;
    filename = fileList(fileIdx).name;

    fprintf('\n[%d/%d] Processing: %s\n', fileIdx, nFiles, filename);

    %-- Initialize log entry
    logEntry = struct();
    logEntry.Filename = filename;
    logEntry.Filepath = filepath;
    logEntry.SubjectID = '';
    logEntry.Group = '';
    logEntry.nChannels = 0;
    logEntry.Duration = 0;
    logEntry.Status = '';
    logEntry.Error = '';
    logEntry.Timestamp = char(datetime('now'));

    try
        %% Load data
        EEG = pop_loadset('filename', filename, 'filepath', filepath);

        % Extract metadata from path/filename
        [subjectID, groupName] = extractMetadata(filename, filepath);

        logEntry.SubjectID = subjectID;
        logEntry.Group = groupName;
        logEntry.nChannels = EEG.nbchan;
        logEntry.Duration = EEG.pnts / EEG.srate;

        fprintf('  Subject: %s | Group: %s | %d ch | %.1f s\n', ...
                subjectID, groupName, EEG.nbchan, EEG.pnts/EEG.srate);

        %% Validate and locate task markers
        [taskMarkerIdx, ~] = findTaskMarkers(EEG, settings.task_names);

        if any(taskMarkerIdx == 0)
            missingTasks = settings.task_names(taskMarkerIdx == 0);
            warning('Missing task markers: %s. Skipping.', strjoin(missingTasks, ', '));
            logEntry.Status = 'Skipped';
            logEntry.Error = sprintf('Missing markers: %s', strjoin(missingTasks, ', '));
            logEntries{end+1} = logEntry; %#ok<AGROW>
            continue
        end

        %% Band Power Analysis — per-channel, no averaging
        fprintf('  → Computing band power (per-channel)...\n');

        % All channels
        bpAllchan = extractBandPower(EEG, settings, taskMarkerIdx, subjectID, groupName);
        bp_cells.allchan{end+1} = bpAllchan;

        % ROIs
        for r = 1:length(roiNames)
            roiName = roiNames{r};
            roiChans = settings.rois.(roiName);
            bpROI = extractBandPower(EEG, settings, taskMarkerIdx, ...
                                      subjectID, groupName, roiChans);
            bp_cells.(roiName){end+1} = bpROI;
        end

        %% ERP Analysis — per-trial, no averaging
        fprintf('  → Computing ERPs (per-trial)...\n');

        % All channels
        [erpAllchan, wfAllchan] = extractERP(EEG, settings, taskMarkerIdx, ...
                                              subjectID, groupName);
        erp_cells.allchan{end+1} = erpAllchan;
        wf_cells.allchan{end+1} = wfAllchan;

        % ROIs
        for r = 1:length(roiNames)
            roiName = roiNames{r};
            roiChans = settings.rois.(roiName);
            [erpROI, wfROI] = extractERP(EEG, settings, taskMarkerIdx, ...
                                          subjectID, groupName, roiChans);
            erp_cells.(roiName){end+1} = erpROI;
            wf_cells.(roiName){end+1} = wfROI;
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

% Compute ERP time vector once (needed for .mat output)
epochSamples = round(settings.erp_epoch * 256);  % placeholder srate
% Actually recompute from any successfully loaded EEG — use settings
nEpochPts = round(diff(settings.erp_epoch) * 256) + 1;  % approximate
% The exact time vector is embedded in the waveform cells; we'll take
% it from the settings directly
erpTimes = linspace(settings.erp_epoch(1)*1000, settings.erp_epoch(2)*1000, ...
                     round(diff(settings.erp_epoch) * 256) + 1);
% NOTE: actual nTimepoints comes from the data; we'll get it from
% the waveform matrix after concatenation

for i = 1:length(outputTypes)
    outType = outputTypes{i};

    % --- Bandpower .mat ---
    if ~isempty(bp_cells.(outType))
        bp_table = vertcat(bp_cells.(outType){:});
        settingsOut = settings; %#ok<NASGU>

        bpFile = fullfile(outputDir, sprintf('bandpower_%s.mat', outType));
        save(bpFile, 'bp_table', 'settingsOut', '-v7.3');

        bpCsv = fullfile(outputDir, sprintf('bandpower_%s.csv', outType));
        writetable(bp_table, bpCsv);

        fprintf('  ✓ bandpower_%s (.mat + .csv, %d rows, %d subjects)\n', ...
                outType, height(bp_table), length(unique(bp_table.SubjectID)));
    end

    % --- ERP ---
    if ~isempty(erp_cells.(outType))
        erp_table = vertcat(erp_cells.(outType){:});

        % Concatenate waveform matrices
        wfList = wf_cells.(outType);
        nonEmpty = ~cellfun(@isempty, wfList);
        if any(nonEmpty)
            wfData = vertcat(wfList{nonEmpty});
            waveform_matrix = vertcat(wfData.waveform); %#ok<NASGU>
            waveform_info = removevars(struct2table(wfData), {'waveform','times'}); %#ok<NASGU>
            times = wfData(1).times; %#ok<NASGU>
        else
            waveform_matrix = []; %#ok<NASGU>
            waveform_info = table(); %#ok<NASGU>
            times = erpTimes; %#ok<NASGU>
        end

        settingsOut = settings; %#ok<NASGU>
        erpFile = fullfile(outputDir, sprintf('erp_%s.mat', outType));
        save(erpFile, 'erp_table', 'waveform_matrix', 'waveform_info', ...
             'times', 'settingsOut', '-v7.3');

        % ERP metrics CSV
        erpCsv = fullfile(outputDir, sprintf('erp_%s.csv', outType));
        writetable(erp_table, erpCsv);

        % Waveform CSV: info columns + time-point columns
        if ~isempty(waveform_matrix)
            timeColNames = arrayfun(@(t) sprintf('t%.1f', t), times, 'UniformOutput', false);
            wfTable = [waveform_info, ...
                       array2table(waveform_matrix, 'VariableNames', timeColNames)];
            wfCsv = fullfile(outputDir, sprintf('waveforms_%s.csv', outType));
            writetable(wfTable, wfCsv);
            fprintf('  ✓ erp_%s (.mat + .csv + waveforms, %d trials, %d subjects)\n', ...
                    outType, height(erp_table), length(unique(erp_table.SubjectID)));
        else
            fprintf('  ✓ erp_%s (.mat + .csv, %d trials, %d subjects)\n', ...
                    outType, height(erp_table), length(unique(erp_table.SubjectID)));
        end
    end
end

% Write processing log
if ~isempty(logEntries)
    logTable = struct2table([logEntries{:}]);
    logFilename = sprintf('analysis_log_%s.csv', ...
                          datetime('now', 'Format', 'yyyy-MM-dd_HH-mm'));
    logFile = fullfile(outputDir, logFilename);
    writetable(logTable, logFile);
    fprintf('  ✓ %s\n', logFilename);
end

fprintf('\n===== EEG ANALYSIS COMPLETE =====\n');
end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function merged = mergeStructs(defaults, overrides)
%MERGESTRUCTS Merge two structs, with overrides taking precedence
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
%EXTRACTMETADATA Extract subject ID and group from filename/path

    % Extract subject ID (TBI### or TRE###)
    tokens = regexp(filename, '(TBI\d{3}|TRE\d{3})', 'match');
    if isempty(tokens)
        tokens = regexp(filepath, '(TBI\d{3}|TRE\d{3})', 'match');
    end

    if isempty(tokens)
        subjectID = 'UNKNOWN';
    else
        subjectID = tokens{1};
    end

    % Extract group from path
    groupPatterns = {'low_performers', 'normal_performers', 'healthy_controls', ...
                     'impaired', 'maintained', 'healthy'};
    groupName = 'unknown';

    for i = 1:length(groupPatterns)
        if contains(filepath, groupPatterns{i})
            groupName = groupPatterns{i};
            break;
        end
    end

    % Normalize group names
    groupMap = containers.Map(...
        {'impaired', 'maintained', 'healthy'}, ...
        {'low_performers', 'normal_performers', 'healthy_controls'});
    if groupMap.isKey(groupName)
        groupName = groupMap(groupName);
    end
end

function [markerIdx, latencies] = findTaskMarkers(EEG, taskNames)
%FINDTASKMARKERS Locate task markers in EEG event structure

    nTasks = length(taskNames);
    markerIdx = zeros(1, nTasks);
    latencies = zeros(1, nTasks);

    if isempty(EEG.event)
        return;
    end

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
        else
            warning('No marker found for task: %s', taskName);
        end
    end
end

%% ========================================================================
%  BAND POWER — PER-CHANNEL (NO AVERAGING)
%  ========================================================================

function bpTable = extractBandPower(EEG, settings, taskMarkerIdx, subjectID, groupName, channelSubset)
%EXTRACTBANDPOWER Compute band power per channel per task — no averaging
%
%   Returns a TALL table with one row per channel × task.
%   Each row has absolute, relative, and log power for all bands.
%
%   Output columns:
%       SubjectID, Group, Channel, Task,
%       delta_abs, theta_abs, alpha_abs, beta_abs, gamma_abs, broadband_abs,
%       delta_rel, theta_rel, ..., delta_log, theta_log, ...

arguments
    EEG             (1,1) struct
    settings        (1,1) struct
    taskMarkerIdx   (1,:) double
    subjectID       (1,1) string
    groupName       (1,1) string
    channelSubset   cell = {}
end

%% Resolve channel indices
if isempty(channelSubset)
    chanIdx = 1:EEG.nbchan;
    chanLabels = {EEG.chanlocs.labels};
else
    allLabels = {EEG.chanlocs.labels};
    chanIdx = [];
    chanLabels = {};
    for i = 1:length(channelSubset)
        idx = find(strcmpi(allLabels, channelSubset{i}), 1);
        if isempty(idx)
            warning('Channel "%s" not found, skipping', channelSubset{i});
        else
            chanIdx(end+1) = idx; %#ok<AGROW>
            chanLabels{end+1} = allLabels{idx}; %#ok<AGROW>
        end
    end
    if isempty(chanIdx)
        error('No valid channels found in channelSubset');
    end
end

%% Setup
bandNames = fieldnames(settings.bands);
nBands = length(bandNames);
nChans = length(chanIdx);
nTasks = length(settings.task_names);

% Pre-allocate cell array for rows (nChans × nTasks rows total)
nRows = nChans * nTasks;
rowIdx = 0;

% Column names for band powers
absColNames = cellfun(@(b) [b '_abs'], bandNames, 'UniformOutput', false);
relColNames = cellfun(@(b) [b '_rel'], bandNames, 'UniformOutput', false);
logColNames = cellfun(@(b) [b '_log'], bandNames, 'UniformOutput', false);

% Pre-allocate numeric arrays
absData = NaN(nRows, nBands);
relData = NaN(nRows, nBands);
logData = NaN(nRows, nBands);
subjectIDs = repmat({char(subjectID)}, nRows, 1);
groupNames = repmat({char(groupName)}, nRows, 1);
channels = cell(nRows, 1);
tasks = cell(nRows, 1);

%% Compute power per channel per task
for t = 1:nTasks
    if taskMarkerIdx(t) == 0
        continue
    end

    taskName = matlab.lang.makeValidName(settings.task_names{t});
    taskStart = round(EEG.event(taskMarkerIdx(t)).latency);
    taskEnd = min(EEG.pnts, taskStart + round(settings.task_duration * EEG.srate));
    taskSegment = EEG.data(chanIdx, taskStart:taskEnd);

    for ch = 1:nChans
        rowIdx = rowIdx + 1;

        % PSD via Welch
        winLen = min(2 * EEG.srate, size(taskSegment, 2));
        [pxx, freqs] = pwelch(taskSegment(ch, :), winLen, winLen/2, [], EEG.srate);
        totalPower = sum(pxx);

        % Band powers
        for b = 1:nBands
            freqRange = settings.bands.(bandNames{b});
            if iscell(freqRange), freqRange = [freqRange{:}]; end
            freqIdx = (freqs >= freqRange(1)) & (freqs <= freqRange(2));

            absPow = sum(pxx(freqIdx));
            absData(rowIdx, b) = absPow;
            relData(rowIdx, b) = absPow / totalPower;
            logData(rowIdx, b) = 10 * log10(absPow + eps);
        end

        % Metadata
        channels{rowIdx} = chanLabels{ch};
        tasks{rowIdx} = taskName;
    end
end

% Trim unused rows (from skipped tasks)
validRows = 1:rowIdx;
absData = absData(validRows, :);
relData = relData(validRows, :);
logData = logData(validRows, :);

%% Build output table
bpTable = table(subjectIDs(validRows), groupNames(validRows), ...
                channels(validRows), tasks(validRows), ...
                'VariableNames', {'SubjectID', 'Group', 'Channel', 'Task'});

% Add band columns
for b = 1:nBands
    bpTable.(absColNames{b}) = absData(:, b);
    bpTable.(relColNames{b}) = relData(:, b);
    bpTable.(logColNames{b}) = logData(:, b);
end

% Convert categorical columns
bpTable.SubjectID = categorical(bpTable.SubjectID);
bpTable.Group = categorical(bpTable.Group);
bpTable.Channel = categorical(bpTable.Channel);
bpTable.Task = categorical(bpTable.Task);
end

%% ========================================================================
%  ERP — PER-TRIAL (NO AVERAGING)
%  ========================================================================

function [erpTable, wfData] = extractERP(EEG, settings, taskMarkerIdx, subjectID, groupName, channelSubset)
%EXTRACTERP Compute per-trial ERP metrics — no averaging across trials
%
%   Returns:
%       erpTable - Tall table with one row per trial × task
%           Columns: SubjectID, Group, Task, Trial, Accuracy,
%                    PeakAmp, PeakLat, MeanAmp, CentLat
%       wfData   - Struct array for waveform storage
%           Fields: SubjectID, Group, Task, Trial, Accuracy,
%                   waveform [1×nPts], times [1×nPts]

arguments
    EEG             (1,1) struct
    settings        (1,1) struct
    taskMarkerIdx   (1,:) double
    subjectID       (1,1) string
    groupName       (1,1) string
    channelSubset   cell = {}
end

%% Resolve channel indices
if isempty(channelSubset)
    chanIdx = 1:EEG.nbchan;
else
    allLabels = {EEG.chanlocs.labels};
    chanIdx = [];
    for i = 1:length(channelSubset)
        idx = find(strcmpi(allLabels, channelSubset{i}), 1);
        if isempty(idx)
            warning('Channel "%s" not found, skipping', channelSubset{i});
        else
            chanIdx(end+1) = idx; %#ok<AGROW>
        end
    end
    if isempty(chanIdx)
        error('No valid channels found in channelSubset');
    end
end

%% Setup
nTasks = length(settings.task_names);
epochSamples = round(settings.erp_epoch * EEG.srate);
baselineSamples = round(settings.erp_baseline * EEG.srate);
nEpochPts = epochSamples(2) - epochSamples(1) + 1;
times = linspace(settings.erp_epoch(1)*1000, settings.erp_epoch(2)*1000, nEpochPts);

% P300 window
p300Idx = (times >= settings.p300_window(1)*1000) & (times <= settings.p300_window(2)*1000);
p300Times = times(p300Idx);

% Event types
if isempty(EEG.event)
    eventTypes = {};
else
    eventTypes = {EEG.event.type};
end

%% Pre-allocate (generous estimate, trim later)
maxTrials = 500 * nTasks;  % upper bound
subjectIDs = cell(maxTrials, 1);
groupNames = cell(maxTrials, 1);
taskLabels = cell(maxTrials, 1);
trialNums  = zeros(maxTrials, 1);
accuracies = cell(maxTrials, 1);
peakAmps   = NaN(maxTrials, 1);
peakLats   = NaN(maxTrials, 1);
meanAmps   = NaN(maxTrials, 1);
centLats   = NaN(maxTrials, 1);
waveforms  = NaN(maxTrials, nEpochPts);

rowIdx = 0;

%% Extract per-trial data
for t = 1:nTasks
    if taskMarkerIdx(t) == 0
        continue
    end

    taskName = settings.task_names{t};
    taskLabel = matlab.lang.makeValidName(taskName);
    stimPattern = [taskName '_stim'];

    stimIdx = find(contains(eventTypes, stimPattern));
    if isempty(stimIdx)
        continue
    end

    trialCount = 0;

    for s = 1:length(stimIdx)
        evtIdx = stimIdx(s);
        lat = round(EEG.event(evtIdx).latency);

        epochStart = lat + epochSamples(1);
        epochEnd   = lat + epochSamples(2);

        % Bounds check
        if epochStart < 1 || epochEnd > EEG.pnts
            continue
        end

        % Extract epoch (selected channels)
        epochData = EEG.data(chanIdx, epochStart:epochEnd);

        % Baseline correction
        baselineStart = 1;
        baselineEnd = baselineStart + (baselineSamples(2) - baselineSamples(1));
        baselineEnd = min(baselineEnd, size(epochData, 2));
        baselineMean = mean(epochData(:, baselineStart:baselineEnd), 2);
        epochData = epochData - baselineMean;

        % Artifact rejection
        if max(abs(epochData(:))) > settings.epoch_thresh
            continue
        end

        % Average across channels for this single trial
        % (spatial averaging within ROI is appropriate; we preserve TRIALS)
        trialWaveform = mean(epochData, 1);  % [1 × nEpochPts]

        % Determine accuracy
        accuracy = 'incorrect';
        if isfield(EEG.event(evtIdx), 'expected') && isfield(EEG.event(evtIdx), 'response')
            expected = EEG.event(evtIdx).expected;
            response = EEG.event(evtIdx).response;

            if isnumeric(expected) && isnumeric(response)
                isCorrect = (expected == response);
            elseif ischar(expected) || isstring(expected)
                isCorrect = strcmpi(string(expected), string(response));
            else
                isCorrect = isequal(expected, response);
            end

            if isCorrect
                accuracy = 'correct';
            end
        end

        % Compute single-trial P300 metrics
        p300Data = trialWaveform(p300Idx);

        [peakAmp, maxI] = max(p300Data);
        peakLat = p300Times(maxI);
        meanAmp = mean(p300Data);

        % Center-of-mass latency
        weights = p300Data - min(p300Data) + eps;
        centLat = sum(p300Times .* weights) / sum(weights);

        % Store
        trialCount = trialCount + 1;
        rowIdx = rowIdx + 1;

        subjectIDs{rowIdx} = char(subjectID);
        groupNames{rowIdx} = char(groupName);
        taskLabels{rowIdx} = taskLabel;
        trialNums(rowIdx)  = trialCount;
        accuracies{rowIdx} = accuracy;
        peakAmps(rowIdx)   = peakAmp;
        peakLats(rowIdx)   = peakLat;
        meanAmps(rowIdx)   = meanAmp;
        centLats(rowIdx)   = centLat;
        waveforms(rowIdx, :) = trialWaveform;
    end
end

%% Trim and build output
validRows = 1:rowIdx;

if rowIdx == 0
    erpTable = table();
    wfData = struct([]);
    return
end

erpTable = table(...
    categorical(subjectIDs(validRows)), ...
    categorical(groupNames(validRows)), ...
    categorical(taskLabels(validRows)), ...
    trialNums(validRows), ...
    categorical(accuracies(validRows)), ...
    peakAmps(validRows), ...
    peakLats(validRows), ...
    meanAmps(validRows), ...
    centLats(validRows), ...
    'VariableNames', {'SubjectID', 'Group', 'Task', 'Trial', 'Accuracy', ...
                      'PeakAmp', 'PeakLat', 'MeanAmp', 'CentLat'});

% Waveform struct array (for .mat storage)
wfData = struct();
for r = 1:rowIdx
    wfData(r).SubjectID = subjectIDs{r};
    wfData(r).Group     = groupNames{r};
    wfData(r).Task      = taskLabels{r};
    wfData(r).Trial     = trialNums(r);
    wfData(r).Accuracy  = accuracies{r};
    wfData(r).waveform  = waveforms(r, :);
    wfData(r).times     = times;
end
wfData = wfData(:);  % column vector
end
