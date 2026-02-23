function EEG_cleanRawData(options)
%EEG_CLEANRAWDATA Preprocess raw EEG data using modified PREP pipeline

arguments
    options.InputDir        (1,1) string
    options.OutputDir       (1,1) string  = fullfile(pwd, "clean", "EEG")
    options.PerformanceDir  (1,1) string  = ""
    options.Hierarchy       (1,:) cell    = {'group', 'subject', 'visit'}
    options.ChannelLabels   (1,:) cell    = {'AFp1','AFp2','AFF5h',     ...
                                             'AFF1h','AFF2h','AFF6h',   ...
                                             'FFC5h','FFC3h','FFC4h',   ...
                                             'FFC6h','FCC5h','FCC3h',   ...
                                             'FCC4h','FCC6h','CCP3h','A1'}
    options.ChannelLayout   (1,1) string  = ["plugins/dipfit/" + ...
                                             "standard_BEM/elec/" + ...
                                             "standard_1005.elc"]
    options.BandpassFilter  (1,2) double  = [0.1 55]
    options.LineFrequency   (1,1) double  = 60
    options.TargetSrate     (1,1) double  = 256
    options.TaskNames       (1,:) cell    = {'eyes_open','nback_0a', ...
                                             'nback_1a', 'nback_0b', ...
                                             'nback_2a'}
    options.TaskDuration    (1,1) double  = 64
    options.PerfTaskMap      (1,1) struct  = struct('eyes_open', '')
    options.UseFixedInterval (1,1) logical = true
    options.FixedInterval    (1,1) double  = 2
    options.ASRBurstCrit    (1,1) double  = 20
    options.ICAThreshold    (1,1) double  = 0.8
    options.Overwrite       (1,1) logical = false
    options.FileExtension   (1,1) string  = "hdf5"
    options.SaveFigures     (1,1) logical = true
end

disp('___________________________________________________________________')
disp(' ')
disp('                     EEG PREPROCESSING STARTED                     ')
disp('___________________________________________________________________')

%% ========================================================================
%  STEP 1: Setup directories and scan for files, init. logging
%  ========================================================================
fprintf('[1/9]: Setting up directories and scanning for files...\n');
eeglab nogui; ALLEEG = [];
if ~exist(options.OutputDir, 'dir'), mkdir(options.OutputDir); end

leafDirs = getLeafs(options.InputDir);
fprintf('       → Found %d leaf directories\n', length(leafDirs));

% Build file list with metadata from hierarchy
fileList = struct('filepath', {}, 'filename', {}, 'relpath', {}, ...
                  'group', {}, 'subject', {}, 'visit', {},'outputpath',{});

for leaf = 1:length(leafDirs)
    files=dir(fullfile(leafDirs{leaf},['*.' char(options.FileExtension)]));
    if isempty(files), continue; end
    
    relPath = strrep(leafDirs{leaf}, options.InputDir, '');
    relPath = strip(relPath, 'left', filesep);
    pathParts = strsplit(relPath, filesep);
    pathParts = pathParts(~cellfun(@isempty, pathParts));
    
    % Map path parts to hierarchy levels
    metadata = struct('group', '', 'subject', '', 'visit', '');
    for h = 1:min(length(options.Hierarchy), length(pathParts))
        metadata.(options.Hierarchy{h}) = pathParts{h};
    end
    
    % Add each file to the list
    for f = 1:length(files)
        idx = length(fileList) + 1;
        fileList(idx).filepath = leafDirs{leaf};
        fileList(idx).filename = files(f).name;
        fileList(idx).relpath = relPath;
        fileList(idx).group = metadata.group;
        fileList(idx).subject = metadata.subject;
        fileList(idx).visit = metadata.visit;
        fileList(idx).outputpath = fullfile(options.OutputDir, relPath);
    end
end

fprintf('       → Found %d EEG files to process\n', length(fileList));
if isempty(fileList), warning('No EEG files found. Exiting.'); return; end

logFile = fullfile(options.OutputDir,sprintf('preprocessing_log_%s.csv',...
                   datetime('now', 'Format', 'yyyy-MM-dd_HH-mm')));
logEntries = {};

%% ========================================================================
%  MAIN PROCESSING LOOP
%  ========================================================================
nFiles = length(fileList);
nTasks = length(options.TaskNames);

for fileIdx = 1:nFiles
    logEntry = initLogEntry(fileList(fileIdx));
    fprintf('\n========================================\n');
    fprintf('[%d/%d] %s\n', fileIdx, nFiles, fileList(fileIdx).filename);
    fprintf('        Group: %s | Subject: %s | Visit: %s\n', ...
            fileList(fileIdx).group, fileList(fileIdx).subject, ...
            fileList(fileIdx).visit);
    fprintf('========================================\n');
    
    [~, baseName, ~] = fileparts(fileList(fileIdx).filename);
    outputFile = fullfile(fileList(fileIdx).outputpath, [baseName '.set']);
    
    if exist(outputFile, 'file') && ~options.Overwrite
        response = input(sprintf(['File exists: %s\nOverwrite? ' ...
                                  '[y]/n/f(abort all): '],outputFile),'s');
        if strcmpi(response, 'f')
            fprintf('Aborting all processing.\n');
            break;
        elseif ~strcmpi(response, 'y') && ~isempty(response)
            fprintf('       → Skipping (user chose not to overwrite)\n');
            logEntry.Status = 'Skipped';
            logEntry.Notes = 'Output exists - no overwrite';
            logEntries{end+1} = logEntry; %#ok<AGROW>
            continue;
        end
    end
    
    udir = fileList(fileIdx).outputpath;
    if ~exist(udir, 'dir'), mkdir(udir); end
    
    try
        %% STEP 2: Load and configure
        fprintf('[2/9]: Loading and configuring...\n');
        
        EEG = pop_loadhdf5('filename', fileList(fileIdx).filename,    ...
                           'filepath',char(fileList(fileIdx).filepath), ...
                           'rejectchans', []);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
        
        %-- Log original
        EEG_raw = EEG;
        logEntry.OrigChannels = EEG.nbchan;
        logEntry.OrigSrate = EEG.srate;
        logEntry.OrigDuration = EEG.pnts / EEG.srate;
        
        % Assign channel labels and locations
        [EEG.chanlocs(1:length(options.ChannelLabels)).labels] = ...
            deal(options.ChannelLabels{:});
        EEG = pop_chanedit(EEG, 'lookup', options.ChannelLayout);
        EEG_original_locs = EEG.chanlocs;
        
        %-- Reset original copy's labels
        lb =arrayfun(@num2str, 1:size(EEG.data,1), 'UniformOutput', false);
        [EEG_raw.chanlocs(1:size(EEG_raw.data,1)).labels] = deal(lb{:});
        
        % Resample
        if EEG.srate > options.TargetSrate
            EEG = pop_resample(EEG, options.TargetSrate);
            EEG_raw = pop_resample(EEG_raw, options.TargetSrate);
            logEntry.ResampledTo = options.TargetSrate;
        else
            logEntry.ResampledTo = EEG.srate;
        end
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        fprintf('       ✓ Data loaded: %d ch, %.1f Hz, %.1f s\n', ...
                EEG.nbchan, EEG.srate, EEG.pnts/EEG.srate);
        
        %% STEP 3: Marker validation (interactive UI)
        fprintf('[3/9]: Validating markers...\n');
        
        % Log original markers
        logEntry.OrigMarkers = length(EEG.event);
        if ~isempty(EEG.event)
            logEntry.OrigMarkerLatencies = num2str([EEG.event.latency]...
                                                    /EEG.srate, '%.1f ');
        else
            logEntry.OrigMarkerLatencies = '[]';
        end
        
        % Check if markers need correction
        if isempty(EEG.event) || length(EEG.event) ~= nTasks
            fprintf('       ✗ Found %d markers, expected %d\n', ...
                    length(EEG.event), nTasks);
            fprintf('       ☷ Opening marker editor...\n');
            EEG = customMarkerEditor(EEG, options.TaskNames);
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            % Set duration for all events
            [EEG.event.duration] = deal(options.TaskDuration * EEG.srate);
            
            if length(EEG.event) ~= nTasks
                warning('Record has %d markers post-edit.',...
                        length(EEG.event));
                response = input(['Continue anyway? Press "t" to ' ...
                                  'terminate [y]/n : '], 's');
                if ~strcmpi(response, 'y') && ~isempty(response)
                    if strcmpi(response, 't'), return, end
                    logEntry.Status = 'Skipped';
                    logEntry.Notes = 'Incorrect markers after editing';
                    logEntries{end+1} = logEntry; %#ok<AGROW>
                    continue;
                end
            else
                fprintf('       ✓ Markers corrected!\n');
            end
        else
            % Assign task names to existing markers
            [EEG.event(1:nTasks).type] = deal(options.TaskNames{:});
            [EEG.event.duration] = deal(options.TaskDuration * EEG.srate);
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            fprintf('       ✓ %d markers validated\n', length(EEG.event));
        end
        
        % [INDEV] Add intermediate markers every 2s for all tasks
        %-- load performance data
        perfData = [];
        if options.PerformanceDir ~= ""
            perfData = loadPerformanceFile(options.PerformanceDir, ...
                                           fileList(fileIdx).subject);
            if ~isempty(perfData)
                fprintf('       ✓ Performance file loaded: %d trials\n',...
                        height(perfData));
                logEntry.Notes = [logEntry.Notes ' PerfFile=loaded'];
            else
                fprintf('       ⚠ No performance file matching ID\n');
            end
        end
        
        fprintf('       → Adding intermediate task markers...\n');
        EEG = addIntermediateMarkers(EEG, options.TaskDuration, ...
                                     perfData, ...
                                     options.PerfTaskMap, ...
                                     options.UseFixedInterval, ...
                                     options.FixedInterval);
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %% STEP 4: Trim to task period
        fprintf('[4/9]: Trimming to task period...\n');
        
        if isempty(EEG.event)
            error('No markers found after validation.');
        end
        
        % Find task boundarie
        taskEvents = EEG.event(~contains({EEG.event.type}, '_stim') & ...
                       ~contains({EEG.event.type}, '_task'));
        if isempty(taskEvents)
            taskEvents = EEG.event;
        end
        
        first_task = taskEvents(1).latency;
        last_task = taskEvents(end).latency;
        buffer_samp = round(1 * EEG.srate);  % 1 second buffer
        task_dur = round(options.TaskDuration * EEG.srate);
        
        trim_start = max(1, first_task - buffer_samp);
        trim_end = min(EEG.pnts, last_task + task_dur + buffer_samp);
        
        % Perform trim
        EEG = pop_select(EEG, 'point',[round(trim_start) round(trim_end)]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        logEntry.TrimStart = trim_start / EEG.srate;
        logEntry.TrimEnd = trim_end / EEG.srate;
        logEntry.TrimmedDuration = EEG.pnts / EEG.srate;
        fprintf('       ✓ Trimmed: %.1f to %.1f s (%.1f s total)\n', ...
            logEntry.TrimStart,logEntry.TrimEnd,logEntry.TrimmedDuration);
        
        %% STEP 5: Line noise removal
        fprintf('[5/9]: Removing line noise...\n');
        lf = options.LineFrequency;
        lineNoiseFreqs = lf:lf:floor(EEG.srate/2);
        EEG = pop_cleanline(EEG,                        ...
            'Bandwidth',                2,              ...
            'ChanCompIndices',          1:EEG.nbchan,   ...
            'SignalType',               'Channels',     ...
            'ComputeSpectralPower',     false,          ...
            'LineFrequencies',          lineNoiseFreqs, ...
            'NormalizeSpectrum',        false,          ...
            'LineAlpha',                0.01,           ...
            'PaddingFactor',            2,              ...
            'PlotFigures',              false,          ...
            'ScanForLines',             true,           ...
            'SmoothingFactor',          100,            ...
            'VerbosityLevel',           0,              ...
            'SlidingWinLength',         4,              ...
            'SlidingWinStep',           1);
        fprintf('       ✓ Line noise removed\n');
        
        %% STEP 6: Robust re-referencing
        fprintf('[6/9]: Robust re-referencing...\n');
        
        % Detect and remove bad channels
        [EEG_clean, ~, ~, rem_chan] = clean_artifacts(EEG, ...
            'FlatlineCriterion', 5, ...
            'ChannelCriterion', 0.8, ...
            'LineNoiseCriterion', 4, ...
            'Highpass', [0.25 0.75], ...
            'BurstCriterion', 'off', ...
            'WindowCriterion', 'off');
        
        logEntry.BadChannels = sum(rem_chan);
        badChanNames = options.ChannelLabels(rem_chan);
        logEntry.BadChanNames = strjoin(badChanNames, ';');
        
        % Interpolate bad channels
        EEG = pop_interp(EEG_clean, EEG_original_locs, 'spherical');
        
        % Average reference
        EEG = pop_reref(EEG, []);
        fprintf('       ✓ Re-referenced (interpolated %d ch: %s)\n', ...
                sum(rem_chan), logEntry.BadChanNames);
        
        %% STEP 7: Filtering and ASR
        fprintf('[7/9]: Bandpass filtering and ASR...\n');
        
        % Create 1Hz highpass copy for ICA
        EEGi = pop_eegfiltnew(EEG, 'locutoff', 1, 'hicutoff', []);
        
        % Bandpass filter main data
        EEG = pop_eegfiltnew(EEG, ...
            'locutoff', options.BandpassFilter(1), ...
            'hicutoff', options.BandpassFilter(2));
        
        % ASR burst correction
        [EEG, ~] = clean_artifacts(EEG, ...
            'FlatlineCriterion', 'off', ...
            'ChannelCriterion', 'off', ...
            'LineNoiseCriterion', 'off', ...
            'BurstCriterion', options.ASRBurstCrit, ...
            'WindowCriterion', 'off');
        fprintf('       ✓ Filtering and ASR complete\n');
        
        %% STEP 8: ICA and ICLabel
        fprintf('[8/9]: ICA decomposition and artifact rejection...\n');
        
        % Run ICA on 1Hz filtered copy
        EEGi = pop_runica(EEGi, 'icatype', 'runica', 'extended', 1, ...
                          'pca', EEG.nbchan - 1);
        
        % Transfer ICA weights
        EEG.icaweights = EEGi.icaweights;
        EEG.icasphere = EEGi.icasphere;
        EEG.icachansind = EEGi.icachansind;
        EEG.icawinv = [];
        EEG = eeg_checkset(EEG);
        clear EEGi;
        
        % Run ICLabel
        EEG = pop_iclabel(EEG, 'default');
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        % Reject artifact components
        reject = [];
        for i = 1:size(EEG.etc.ic_classification.ICLabel.classifications,1)
            iclab=EEG.etc.ic_classification.ICLabel.classifications(i, :);
            prob_artifact = sum(iclab(2:6));
            if prob_artifact > options.ICAThreshold
                reject = [reject i]; %#ok<AGROW>
            end
        end
        
        nComponents = size(EEG.icaweights, 1);
        EEG = pop_subcomp(EEG, reject, 0);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,...
                                               'gui', 'off');
        
        logEntry.ICAComponents = nComponents;
        logEntry.ICARejected = length(reject);
        logEntry.ICARejectedIdx = mat2str(reject);
        
        fprintf('       ✓ Rejected %d/%d ICA components\n', ...
                length(reject), nComponents);
        
        %% STEP 9: Visualization and validation
        fprintf('[9/9]: Generating validation visualization...\n');
        
        visualizeRawVsCleaned(EEG_raw, EEG, options.ChannelLabels,      ...
                              round(trim_start),round(trim_end),        ...
                              rem_chan,baseName);
        
        if options.SaveFigures
            figPath = fullfile(fileList(fileIdx).outputpath, ...
                      [baseName '_validation.png']);
            saveas(gcf, figPath);
            fprintf('       ✓ Validation figure saved: %s\n', figPath);
        end
        
        response = input(['       Review visualization. ' ...
                         'Continue to save? [y]/n : '], 's');
        if strcmpi(response, 'n')
            close(gcf);
            logEntry.Status = 'Rejected';
            logEntry.Notes = 'User rejected after visual review';
            logEntries{end+1} = logEntry; %#ok<AGROW>
            continue;
        end
        close(gcf);
        
        %% Save
        output_path = char(fileList(fileIdx).outputpath);
        EEG = pop_saveset(EEG, 'filename', baseName, ...
                               'filepath', output_path,...
                               'savemode', 'onefile');
        
        logEntry.Status = 'Completed';
        logEntry.FinalChannels = EEG.nbchan;
        logEntry.FinalDuration = EEG.pnts / EEG.srate;
        fprintf('       ✓ Saved: %s\n', outputFile);
        
    catch ME
        logEntry.Status = 'Failed';
        logEntry.Notes = strrep(ME.message, ',', ';');
        fprintf('       ✗ ERROR: %s\n', ME.message);
        fprintf('       Stack trace:\n');
        for k = 1:length(ME.stack)
            fprintf('         %s (line %d)\n', ...
                    ME.stack(k).name, ME.stack(k).line);
        end
    end
    logEntries{end+1} = logEntry; %#ok<AGROW>
end

%% ========================================================================
%  Finalize
%  ========================================================================
% Write log
if ~isempty(logEntries)
    logTable = struct2table([logEntries{:}]);
    writetable(logTable, logFile);
    fprintf('\n✓ Log saved: %s\n', logFile);
end

disp('___________________________________________________________________')
fprintf('\n                   EEG PREPROCESSING COMPLETE\n')
disp('___________________________________________________________________')

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function logEntry = initLogEntry(fileInfo)
%INITLOGENTRY Initialize a log entry structure
    logEntry.Timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    logEntry.Filename = fileInfo.filename;
    logEntry.Group = fileInfo.group;
    logEntry.Subject = fileInfo.subject;
    logEntry.Visit = fileInfo.visit;
    logEntry.Status = 'Started';
    logEntry.OrigChannels = 0;
    logEntry.OrigSrate = 0;
    logEntry.OrigDuration = 0;
    logEntry.UnitScaling = '';
    logEntry.ResampledTo = 0;
    logEntry.OrigMarkers = 0;
    logEntry.OrigMarkerLatencies = '';
    logEntry.TrimStart = 0;
    logEntry.TrimEnd = 0;
    logEntry.TrimmedDuration = 0;
    logEntry.BadChannels = 0;
    logEntry.BadChanNames = '';
    logEntry.ICAComponents = 0;
    logEntry.ICARejected = 0;
    logEntry.ICARejectedIdx = '';
    logEntry.FinalChannels = 0;
    logEntry.FinalDuration = 0;
    logEntry.Notes = '';
end

function perfData = loadPerformanceFile(perfDir, subjectID)
%LOADPERFORMANCEFILE Load performance CSV for a given subject
%
%   Searches for *_COG.csv files containing subject ID
%   Returns table or [] if no file found

    perfData = [];
    
    if perfDir == "" || ~exist(perfDir, 'dir')
        return;
    end
    
    % Hardcoded file pattern
    searchPattern = fullfile(perfDir, '**', '*_COG.csv');
    files = dir(searchPattern);
    
    if isempty(files)
        return;
    end
    
    % Find file containing subject ID
    matchIdx = [];
    for f = 1:length(files)
        if contains(files(f).name, subjectID, 'IgnoreCase', true)
            matchIdx = f;
            break;
        end
    end
    
    if isempty(matchIdx)
        % Try more flexible matching
        for f = 1:length(files)
            [~, fname, ~] = fileparts(files(f).name);
            tokens = regexp(fname, '([A-Z]{2,3}\d{3})', 'match');
            if any(strcmpi(tokens, subjectID))
                matchIdx = f;
                break;
            end
        end
    end
    
    if isempty(matchIdx)
        return;
    end
    
    % Load the CSV
    try
        perfFile = fullfile(files(matchIdx).folder, files(matchIdx).name);
        perfData = readtable(perfFile, 'VariableNamingRule', 'preserve');
        
        % Validate expected columns
        requiredCols = {'StimulusType', 'StimOffset'};
        if ~all(ismember(requiredCols, perfData.Properties.VariableNames))
            warning('Performance file missing required columns. Ignoring.');
            perfData = [];
        end
    catch ME
        perfData = [];
    end
end

function EEG = addIntermediateMarkers(EEG, taskDuration, ...
                                      perfData, perfTaskMap, ...
                                      useFixedInterval, fixedInterval)
%ADDINTERMEDIATEMARKERS Add stimulus markers within each task
%
%   If perfData is provided and task is not mapped to '', uses actual 
%   stimulus onset times from the performance file. Otherwise uses fixed
%   interval markers.

    if isempty(EEG.event), return; end
    mainTaskIdx = find(~contains({EEG.event.type}, '_stim'));
    nMarkersAdded = 0;
    
    for t = 1:length(mainTaskIdx)
        idx = mainTaskIdx(t);
        taskStart = EEG.event(idx).latency;
        taskName = EEG.event(idx).type;
        
        % Check if this task should use fixed intervals (mapped to '')
        fieldName = matlab.lang.makeValidName(taskName);
        useFixed = false;
        if isstruct(perfTaskMap) && isfield(perfTaskMap, fieldName)
            mappedValue = perfTaskMap.(fieldName);
            if isempty(mappedValue) || strcmp(mappedValue, '')
                useFixed = true;
            end
        end
        
        % Check if we have performance data for this task
        hasPerf = false;
        taskTrials = [];
        if ~useFixed && ~isempty(perfData)
            % Use task name directly to find trials in performance file
            taskTrials=perfData(strcmp(perfData.StimulusType, taskName),:);
            if ~isempty(taskTrials) && height(taskTrials) > 0
                hasPerf = true;
            end
        end
        
        if hasPerf
            %-- STIMULUS FILE MARKERS -------------------------------------
            stimOffsets_ms = taskTrials.StimOffset;  % in milliseconds
            
            for s = 1:length(stimOffsets_ms)
                % Convert ms to samples
                stimOffset_samp=round((stimOffsets_ms(s)/1000)*EEG.srate);
                stimLatency = taskStart + stimOffset_samp;
                
                % Only add if within data bounds
                if stimLatency > 0 && stimLatency <= EEG.pnts
                    EEG.event(end+1).type = [taskName '_stim'];
                    EEG.event(end).latency = stimLatency;
                    
                    % Calculate duration to next stimulus (or end of task)
                    if s < length(stimOffsets_ms)
                        nextOffset_samp=round((stimOffsets_ms(s+1)/1000)...
                                              * EEG.srate);
                        EEG.event(end).duration = nextOffset_samp ...
                                                  - stimOffset_samp;
                    else
                        taskEnd_samp = taskDuration * EEG.srate;
                        EEG.event(end).duration = taskEnd_samp ...
                                                  - stimOffset_samp;
                    end
                    
                    EEG.event(end).position = s;
                    if ismember('Stimulus', ...
                                taskTrials.Properties.VariableNames)
                        EEG.event(end).stimulus = taskTrials.Stimulus{s};
                    end
                    if ismember('ExpectedResponse', ...
                                taskTrials.Properties.VariableNames)
                        EEG.event(end).expected = ...
                            taskTrials.ExpectedResponse(s);
                    end
                    if ismember('ActualResponse', ...
                                taskTrials.Properties.VariableNames)
                        resp = taskTrials.ActualResponse(s);
                        if iscell(resp)
                            EEG.event(end).response = resp{1};
                        else
                            EEG.event(end).response = resp;
                        end
                    end
                    if ismember('ReactionTime', ...
                                taskTrials.Properties.VariableNames)
                        EEG.event(end).rt = taskTrials.ReactionTime(s);
                    end
                    
                    nMarkersAdded = nMarkersAdded + 1;
                end
            end
        %-- FIXED INTERVAL MARKERS ----------------------------------------   
        elseif useFixedInterval
            numIntermediate = floor((taskDuration - fixedInterval) / ...
                                    fixedInterval);
            
            for i = 1:numIntermediate
                stimLatency = taskStart + (i * fixedInterval * EEG.srate);
                
                if stimLatency > 0 && stimLatency <= EEG.pnts
                    EEG.event(end+1).type = [taskName '_task'];
                    EEG.event(end).latency = stimLatency;
                    EEG.event(end).duration = fixedInterval * EEG.srate;
                    EEG.event(end).position = i;
                    nMarkersAdded = nMarkersAdded + 1;
                end
            end
        else
            fprintf('         [%s] No intermediate markers added\n', ...
                    taskName);
        end
    end
    
    [~, sortIdx] = sort([EEG.event.latency]);
    EEG.event = EEG.event(sortIdx);
    
    urevents = num2cell(1:length(EEG.event));
    [EEG.event.urevent] = deal(urevents{:});
    
    fprintf('       ✓ Total intermediate markers added: %d\n', ...
            nMarkersAdded);
end

%% ========================================================================
%  INTERACTIVE MARKER EDITOR
%  ========================================================================

function EEG = customMarkerEditor(EEG, taskNames)
%CUSTOMMARKEREDITOR Interactive GUI for editing EEG markers
    
    % Create figure
    fig = figure('Name', 'Marker Editor', 'Position', [50 50 1400 800], ...
                 'NumberTitle', 'off', 'MenuBar', 'none', 'Resize', 'on');
    
    % Ensure we have event structure
    if isempty(EEG.event)
        EEG.event = struct('type', {}, 'latency', {}, 'duration', {});
    end
    
    % GUI Data
    guiData.EEG = EEG;
    guiData.taskNames = taskNames;
    guiData.markerLines = [];
    guiData.nTasks = length(taskNames);
    
    % Main plot area (top 70%)
    guiData.ax = axes('Parent', fig, 'Position', [0.05 0.35 0.9 0.6]);
    
    % Plot EEG data (subset of channels for visibility)
    nChanToPlot = min(8, EEG.nbchan);
    timeVec = (0:EEG.pnts-1) / EEG.srate;
    
    % Plot with offset for visibility
    plotData = EEG.data(1:nChanToPlot, :);
    offset = max(abs(plotData(:))) * 2;
    if offset == 0, offset = 100; end  % Fallback for flat data
    
    hold(guiData.ax, 'on');
    for ch = 1:nChanToPlot
        plot(guiData.ax, timeVec, plotData(ch,:) + offset*(nChanToPlot-ch), 'b');
    end
    
    xlabel(guiData.ax, 'Time (seconds)');
    ylabel(guiData.ax, 'Channels (offset)');
    title(guiData.ax, sprintf('EEG Data - %s | Click to add markers', EEG.filename), ...
          'Interpreter', 'none');
    ylim(guiData.ax, [-offset, offset*nChanToPlot]);
    xlim(guiData.ax, [0, EEG.xmax]);
    grid(guiData.ax, 'on');
    
    % Control panel (bottom 30%)
    panelY = 0.02;
    panelHeight = 0.28;
    
    % Marker latency display and edit
    uicontrol('Style', 'text', 'String', 'Marker Latencies (seconds):', ...
              'Units', 'normalized', 'Position', [0.05 panelY+panelHeight-0.03 0.3 0.025], ...
              'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 10);
    
    % Create text fields and delete buttons for each marker
    guiData.textFields = [];
    guiData.deleteButtons = [];
    guiData.taskLabels = [];
    fieldWidth = 0.09;
    
    for i = 1:guiData.nTasks
        % Task name label
        guiData.taskLabels(i) = uicontrol('Style', 'text', ...
                  'String', sprintf('%s:', taskNames{i}), ...
                  'Units', 'normalized', ...
                  'Position', [0.05 + (i-1)*fieldWidth panelY+panelHeight-0.06 0.08 0.025], ...
                  'HorizontalAlignment', 'center', 'FontSize', 8);
        
        % Editable text field
        guiData.textFields(i) = uicontrol('Style', 'edit', ...
                  'String', '', ...
                  'Units', 'normalized', ...
                  'Position', [0.05 + (i-1)*fieldWidth panelY+panelHeight-0.09 0.08 0.03], ...
                  'HorizontalAlignment', 'center', ...
                  'FontSize', 10, ...
                  'Callback', @(src,evt) updateFromTextField(fig, i));
        
        % Delete button below text field
        guiData.deleteButtons(i) = uicontrol('Style', 'pushbutton', ...
                  'String', 'X', ...
                  'Units', 'normalized', ...
                  'Position', [0.05 + (i-1)*fieldWidth panelY+panelHeight-0.12 0.08 0.025], ...
                  'FontSize', 9, ...
                  'ForegroundColor', [0.8 0 0], ...
                  'FontWeight', 'bold', ...
                  'Callback', @(src,evt) deleteMarkerAtIndex(fig, i));
    end
    
    % Buttons row
    buttonY = panelY + panelHeight - 0.18;
    
    uicontrol('Style', 'pushbutton', 'String', '← Add 74s Before First', ...
              'Units', 'normalized', 'Position', [0.05 buttonY 0.15 0.04], ...
              'FontSize', 10, 'Callback', @(src,evt) addMarkerBeforeFirst(fig));
    
    uicontrol('Style', 'pushbutton', 'String', 'Add 74s After Last →', ...
              'Units', 'normalized', 'Position', [0.21 buttonY 0.15 0.04], ...
              'FontSize', 10, 'Callback', @(src,evt) addMarkerAfterLast(fig));
    
    uicontrol('Style', 'pushbutton', 'String', 'Clear All', ...
              'Units', 'normalized', 'Position', [0.37 buttonY 0.1 0.04], ...
              'FontSize', 10, 'ForegroundColor', [0.8 0 0], ...
              'Callback', @(src,evt) clearAllMarkers(fig));
    
    % Info text
    guiData.infoText = uicontrol('Style', 'text', ...
              'String', sprintf('Current markers: %d/%d', length(EEG.event), guiData.nTasks), ...
              'Units', 'normalized', 'Position', [0.05 buttonY-0.05 0.4 0.03], ...
              'HorizontalAlignment', 'left', 'FontSize', 10);
    
    % Done button
    uicontrol('Style', 'pushbutton', 'String', '✓ DONE (Save & Close)', ...
              'Units', 'normalized', 'Position', [0.75 buttonY 0.2 0.06], ...
              'FontSize', 12, 'FontWeight', 'bold', ...
              'BackgroundColor', [0.2 0.8 0.2], ...
              'Callback', @(src,evt) closeMarkerEditor(fig));
    
    % Click on plot to add marker
    set(guiData.ax, 'ButtonDownFcn', @(src,evt) addMarkerByClick(fig));
    
    % Store GUI data
    guidata(fig, guiData);
    
    % Initial plot update
    updateMarkerPlot(fig);
    
    % Wait for figure to close
    uiwait(fig);
    
    % Retrieve modified EEG
    if ishandle(fig)
        guiData = guidata(fig);
        EEG = guiData.EEG;
        close(fig);
    end
end

function newEvent = createMarkerEvent(EEG, latency)
%CREATEMARKEREVENT Create a new event structure matching existing format
    if ~isempty(EEG.event)
        newEvent = EEG.event(1);
        fields = fieldnames(newEvent);
        for i = 1:length(fields)
            switch fields{i}
                case 'type'
                    newEvent.type = 'marker';
                case 'latency'
                    newEvent.latency = latency;
                case 'duration'
                    newEvent.duration = 64;
                case 'urevent'
                    newEvent.urevent = [];
                case 'position'
                    newEvent.position = 1;
                otherwise
                    if isnumeric(newEvent.(fields{i}))
                        newEvent.(fields{i}) = 0;
                    else
                        newEvent.(fields{i}) = '';
                    end
            end
        end
    else
        newEvent.type = 'marker';
        newEvent.latency = latency;
        newEvent.duration = 64;
    end
end

function updateMarkerPlot(fig)
%UPDATEMARKERPLOT Update plot with current marker positions
    guiData = guidata(fig);
    EEG = guiData.EEG;
    
    % Delete old marker lines
    if ~isempty(guiData.markerLines)
        delete(guiData.markerLines(ishandle(guiData.markerLines)));
    end
    guiData.markerLines = [];
    
    % Draw new marker lines
    if ~isempty(EEG.event)
        for i = 1:length(EEG.event)
            latency_s = EEG.event(i).latency / EEG.srate;
            guiData.markerLines(i) = xline(guiData.ax, latency_s, 'r-', ...
                sprintf('M%d: %.1fs', i, latency_s), ...
                'LineWidth', 2, 'FontSize', 9, 'LabelOrientation', 'horizontal');
        end
    end
    
    % Update text fields and delete buttons
    for i = 1:guiData.nTasks
        if i <= length(EEG.event)
            latency_s = EEG.event(i).latency / EEG.srate;
            set(guiData.textFields(i), 'String', sprintf('%.2f', latency_s));
            set(guiData.deleteButtons(i), 'Enable', 'on');
        else
            set(guiData.textFields(i), 'String', '');
            set(guiData.deleteButtons(i), 'Enable', 'off');
        end
    end
    
    % Update info text
    set(guiData.infoText, 'String', ...
        sprintf('Current markers: %d/%d', length(EEG.event), guiData.nTasks));
    
    guidata(fig, guiData);
end

function addMarkerBeforeFirst(fig)
%ADDMARKERBEFOREFIRST Add marker 74s before the first marker
    guiData = guidata(fig);
    EEG = guiData.EEG;
    
    if isempty(EEG.event)
        warndlg('No markers exist! Click on plot to add first marker.', 'No Markers');
        return;
    end
    
    firstLatency = EEG.event(1).latency;
    newLatency = firstLatency - 74 * EEG.srate;
    
    if newLatency < 1
        warndlg('Cannot add marker before time 0!', 'Invalid Position');
        return;
    end
    
    newEvent = createMarkerEvent(EEG, newLatency);
    EEG.event = [newEvent, EEG.event];
    guiData.EEG = EEG;
    guidata(fig, guiData);
    updateMarkerPlot(fig);
end

function addMarkerAfterLast(fig)
%ADDMARKERAFTERLAST Add marker 74s after the last marker
    guiData = guidata(fig);
    EEG = guiData.EEG;
    
    if isempty(EEG.event)
        warndlg('No markers exist! Click on plot to add first marker.', 'No Markers');
        return;
    end
    
    lastLatency = EEG.event(end).latency;
    newLatency = lastLatency + 74 * EEG.srate;
    
    if newLatency > EEG.pnts
        warndlg('Cannot add marker beyond recording end!', 'Invalid Position');
        return;
    end
    
    newEvent = createMarkerEvent(EEG, newLatency);
    EEG.event(end+1) = newEvent;
    guiData.EEG = EEG;
    guidata(fig, guiData);
    updateMarkerPlot(fig);
end

function deleteMarkerAtIndex(fig, markerIdx)
%DELETEMARKERATINDEX Delete the marker at the specified index
    guiData = guidata(fig);
    EEG = guiData.EEG;
    
    if markerIdx > length(EEG.event)
        return;
    end
    
    EEG.event(markerIdx) = [];
    guiData.EEG = EEG;
    guidata(fig, guiData);
    updateMarkerPlot(fig);
end

function clearAllMarkers(fig)
%CLEARALLMARKERS Remove all markers after confirmation
    response = questdlg('Delete all markers?', 'Confirm', 'Yes', 'No', 'No');
    if strcmp(response, 'Yes')
        guiData = guidata(fig);
        
        if ~isempty(guiData.EEG.event)
            fields = fieldnames(guiData.EEG.event);
            emptyStruct = cell2struct(cell(length(fields), 0), fields, 1);
            guiData.EEG.event = emptyStruct;
        else
            guiData.EEG.event = struct('type', {}, 'latency', {}, 'duration', {});
        end
        
        guidata(fig, guiData);
        updateMarkerPlot(fig);
    end
end

function updateFromTextField(fig, markerIdx)
%UPDATEFROMTEXTFIELD Update marker position from text field edit
    guiData = guidata(fig);
    EEG = guiData.EEG;
    
    if markerIdx > length(EEG.event)
        return;
    end
    
    newValue_s = str2double(get(guiData.textFields(markerIdx), 'String'));
    
    if isnan(newValue_s)
        warndlg('Invalid number!', 'Error');
        updateMarkerPlot(fig);
        return;
    end
    
    newLatency = newValue_s * EEG.srate;
    
    if newLatency < 1 || newLatency > EEG.pnts
        warndlg('Marker position outside recording!', 'Error');
        updateMarkerPlot(fig);
        return;
    end
    
    EEG.event(markerIdx).latency = newLatency;
    guiData.EEG = EEG;
    guidata(fig, guiData);
    updateMarkerPlot(fig);
end

function addMarkerByClick(fig)
%ADDMARKERBYCLICK Add marker at clicked position on plot
    guiData = guidata(fig);
    EEG = guiData.EEG;
    
    point = get(guiData.ax, 'CurrentPoint');
    clickTime_s = point(1,1);
    
    if clickTime_s < 0 || clickTime_s > EEG.xmax
        return;
    end
    
    if length(EEG.event) >= guiData.nTasks
        warndlg(sprintf('Already have %d markers! Delete one first.', guiData.nTasks), ...
                'Too Many Markers');
        return;
    end
    
    newEvent = createMarkerEvent(EEG, clickTime_s * EEG.srate);
    
    % Insert in chronological order
    if isempty(EEG.event)
        EEG.event = newEvent;
    else
        insertIdx = length(EEG.event) + 1;
        for i = 1:length(EEG.event)
            if EEG.event(i).latency > newEvent.latency
                insertIdx = i;
                break;
            end
        end
        
        if insertIdx > length(EEG.event)
            EEG.event(end+1) = newEvent;
        else
            EEG.event = [EEG.event(1:insertIdx-1), newEvent, EEG.event(insertIdx:end)];
        end
    end
    
    guiData.EEG = EEG;
    guidata(fig, guiData);
    updateMarkerPlot(fig);
end

function closeMarkerEditor(fig)
%CLOSEMARKEREDITOR Finalize and close the marker editor
    guiData = guidata(fig);
    EEG = guiData.EEG;
    
    % Assign task names if we have the correct number of markers
    if length(EEG.event) == guiData.nTasks
        for i = 1:guiData.nTasks
            EEG.event(i).type = guiData.taskNames{i};
        end
    end
    
    guiData.EEG = EEG;
    guidata(fig, guiData);
    uiresume(fig);
end

%% ========================================================================
%  VISUALIZATION: RAW VS CLEANED COMPARISON
%  ========================================================================

function visualizeRawVsCleaned(EEG_raw, EEG_cleaned, channelLabels, ...
                               trim_start_idx, trim_end_idx, interpolated_channels, titleStr)
%VISUALIZERAWVSCLEANED Display raw (top) and cleaned (bottom) EEG for validation
%
%   Inputs:
%       EEG_raw               - Original raw EEG dataset (before trim)
%       EEG_cleaned           - Final processed EEG dataset
%       channelLabels         - Cell array of channel labels
%       trim_start_idx        - Sample index where trimming started
%       trim_end_idx          - Sample index where trimming ended
%       interpolated_channels - Logical vector or indices of interpolated channels
%       titleStr              - Title for the figure

    % Setup
    n_chans = min(EEG_raw.nbchan, length(channelLabels));
    times_raw = (0:EEG_raw.pnts-1) / EEG_raw.srate;
    
    % Calculate channel spacing for waterfall plot
    spacing = 4 * mean(std(EEG_raw.data(1:n_chans,:), 0, 2));
    if spacing == 0 || isnan(spacing), spacing = 50; end
    offsets = (0:n_chans-1)' * spacing;
    
    % Create figure
    fig = figure('Name', ['Validation: ' titleStr], 'Color', 'w', ...
                 'Position', [50, 50, 1600, 900]);
    
    %% TOP PANEL: Raw (unprocessed) data
    ax1 = subplot(2, 1, 1);
    hold on;
    
    % Plot raw data - all channels in blue
    raw_plot_data = EEG_raw.data(1:n_chans, :) + offsets;
    for ch = 1:n_chans
        plot(times_raw, raw_plot_data(ch,:), 'b', 'LineWidth', 0.5);
    end
    
    % Highlight trimmed-out regions in red
    if trim_start_idx > 1
        xregion(0, times_raw(trim_start_idx), 'FaceColor', 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end
    if trim_end_idx < EEG_raw.pnts
        xregion(times_raw(trim_end_idx), times_raw(end), 'FaceColor', 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end
    
    % Plot markers on raw data
    plotMarkersOnAxis(ax1, EEG_raw, times_raw, offsets(end) + spacing);
    
    % Formatting
    title(ax1, sprintf('RAW (Unprocessed) - %s', titleStr), 'FontSize', 12, 'FontWeight', 'bold');
    set(ax1, 'YTick', offsets, 'YTickLabel', channelLabels(1:n_chans), 'YDir', 'normal');
    xlabel(ax1, 'Time (s)');
    ylabel(ax1, 'Channels');
    ylim(ax1, [-spacing, offsets(end) + spacing*2]);
    xlim(ax1, [0, times_raw(end)]);
    grid(ax1, 'on');
    
    %% BOTTOM PANEL: Cleaned (processed) data
    ax2 = subplot(2, 1, 2);
    hold on;
    
    % Map cleaned data back to raw timeline
    % The cleaned data corresponds to [trim_start_idx : trim_end_idx] in the raw timeline
    data_cleaned_mapped = nan(n_chans, EEG_raw.pnts);
    
    % Check for ASR sample mask (indicates which samples were kept/modified)
    if isfield(EEG_cleaned.etc, 'clean_sample_mask')
        mask = EEG_cleaned.etc.clean_sample_mask;
        trimmed_len = trim_end_idx - trim_start_idx + 1;
        
        if length(mask) == trimmed_len
            kept_indices_relative = find(mask);
            if length(kept_indices_relative) == size(EEG_cleaned.data, 2)
                kept_indices_absolute = trim_start_idx + kept_indices_relative - 1;
                data_cleaned_mapped(:, kept_indices_absolute) = EEG_cleaned.data(1:n_chans, :);
            else
                % Fallback: place data starting at trim_start
                end_idx = min(trim_start_idx + size(EEG_cleaned.data, 2) - 1, EEG_raw.pnts);
                data_cleaned_mapped(:, trim_start_idx:end_idx) = ...
                    EEG_cleaned.data(1:n_chans, 1:(end_idx-trim_start_idx+1));
            end
        else
            % Mask length mismatch - use simple placement
            end_idx = min(trim_start_idx + size(EEG_cleaned.data, 2) - 1, EEG_raw.pnts);
            data_cleaned_mapped(:, trim_start_idx:end_idx) = ...
                EEG_cleaned.data(1:n_chans, 1:(end_idx-trim_start_idx+1));
        end
    else
        % No mask - assume continuous data in trimmed region
        end_idx = min(trim_start_idx + size(EEG_cleaned.data, 2) - 1, EEG_raw.pnts);
        data_cleaned_mapped(:, trim_start_idx:end_idx) = ...
            EEG_cleaned.data(1:n_chans, 1:(end_idx-trim_start_idx+1));
    end
    
    % Plot faint raw data as background reference
    raw_plot_faint = EEG_raw.data(1:n_chans, :) + offsets;
    for ch = 1:n_chans
        plot(times_raw, raw_plot_faint(ch,:), 'Color', [0.85 0.85 0.85], 'LineWidth', 0.3);
    end
    
    % Plot cleaned data
    clean_plot_data = data_cleaned_mapped + offsets;
    
    % Convert interpolated_channels to logical if needed
    if islogical(interpolated_channels)
        interp_mask = interpolated_channels(1:n_chans);
    else
        interp_mask = false(1, n_chans);
        interp_mask(interpolated_channels(interpolated_channels <= n_chans)) = true;
    end
    
    for ch = 1:n_chans
        if interp_mask(ch)
            % Interpolated channels in magenta
            plot(times_raw, clean_plot_data(ch,:), 'm', 'LineWidth', 0.8);
        else
            % Normal channels in blue
            plot(times_raw, clean_plot_data(ch,:), 'b', 'LineWidth', 0.6);
        end
    end
    
    % Mark trimmed/cut regions in red overlay
    if trim_start_idx > 1
        xregion(0, times_raw(trim_start_idx), 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
    if trim_end_idx < EEG_raw.pnts
        xregion(times_raw(trim_end_idx), times_raw(end), 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
    
    % Formatting
    title(ax2, 'CLEANED (Processed) - Blue = retained, Magenta = interpolated, Red shading = cut', ...
          'FontSize', 12, 'FontWeight', 'bold');
    set(ax2, 'YTick', offsets, 'YTickLabel', channelLabels(1:n_chans), 'YDir', 'normal');
    xlabel(ax2, 'Time (s)');
    ylabel(ax2, 'Channels');
    ylim(ax2, [-spacing, offsets(end) + spacing*2]);
    xlim(ax2, [0, times_raw(end)]);  % Same x-axis as raw for alignment
    grid(ax2, 'on');
    
    % Link axes for synchronized zooming/panning
    linkaxes([ax1, ax2], 'x');
    
    % Add legend
    hold(ax2, 'on');
    dummy_handles = gobjects(4, 1);
    dummy_handles(1) = plot(ax2, NaN, NaN, 'b-', 'LineWidth', 1.5);
    dummy_handles(2) = plot(ax2, NaN, NaN, 'm-', 'LineWidth', 1.5);
    dummy_handles(3) = plot(ax2, NaN, NaN, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 1.5);
    dummy_handles(4) = plot(ax2, NaN, NaN, 'rs', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
    
    legend(dummy_handles, {'Cleaned', 'Interpolated', 'Raw (reference)', 'Cut/Trimmed'}, ...
           'Location', 'southeast', 'FontSize', 9);
end

function plotMarkersOnAxis(ax, EEG, times, labelY)
%PLOTMARKERSONAXIS Plot event markers as vertical lines
    if isempty(EEG.event), return; end
    
    for e = 1:length(EEG.event)
        lat = EEG.event(e).latency;
        eventType = EEG.event(e).type;
        
        % Skip boundary and intermediate markers for cleaner display
        if contains(char(eventType), 'boundary'), continue; end
        if contains(char(eventType), '_stim'), continue; end
        
        if lat >= 1 && lat <= length(times)
            t_sec = times(round(lat));
            xline(ax, t_sec, 'k--', 'LineWidth', 1.5, 'Label', char(eventType), ...
                  'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', ...
                  'Interpreter', 'none', 'FontSize', 8);
        end
    end
end
