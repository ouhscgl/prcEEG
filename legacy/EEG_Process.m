function [EEGStats, demograph] = ...
    EEG_Process(load_path, eeglab_path, user_vars)
% EEG_Process - Main processing core for EEG pipeline.
% 
% Usage:
%   [EEGStats, demograph, events] = EEG_Process(load_path, eeglab_path, 
%                                               user_vars)
%
% Inputs:
%   load_path    - Path to EEG file(s) or folder system
%   eeglab_path  - Path to EEGLAB toolbox
%   user_vars    - Structure with customizable parameters for the pipeline
%
% Outputs:
%   EEGStats     - Structure containing processed EEG data and statistics
%   demograph    - Table with demographic information
%
% Description:
%   This function serves as the main processing core for EEG analysis.
%   It handles loading, preprocessing, and analyzing EEG data with options
%   for frequency band analysis, artifact rejection, ICA, and connectivity.

% Data loading and parameter validation ___________________________________
% -- Logfile setup (needs work)
global logfile
logfile = struct();

% -- adding EEGLab to path
try addpath(eeglab_path), eeglab('nogui')
catch e, error(['Unable to load EEGLab: ',e.message]), end

%-- Default analysis variables
defaults = struct();
defaults.reference_channels = {'CCP4h'}; 
defaults.bandpass_filter = {0.5 45};        
defaults.notch_filter = [50];                
defaults.run_ica = true;                     
defaults.epoch_window = {-0.2 0.8};               
defaults.baseline_window = [-200 0];           
defaults.artifact_threshold = 100;             
defaults.frequency_bands = struct(... 
    'delta',     [0.5  4], ...
    'theta',     [4    8], ...
    'alpha',     [8   13], ...
    'beta',      [13  30], ...
    'gamma',     [30  45], ...
    'broadband', [0.5 45]);
defaults.iclabel_thresholds = [
    NaN NaN;  %
    0.4   1;  %
    0.4   1;  %
    0.4   1;  % 
    0.4   1;  %
    0.4   1;  %
    NaN NaN]; %
defaults.channelLabels = {
    'AFp1',  'AFp2', 'AFF5h', 'AFF1h', 'AFF2h', 'AFF6h', 'FFC5h', ...
    'FFC3h', 'FFC4h', 'FFC6h', 'FCC5h', 'FCC3h','FCC4h', 'FCC6h', ...
    'CCP3h', 'CCP4h'};
defaults.connectivity_metric = 'rt_DCCA';   
defaults.connectivity_bands = {'alpha', 'theta'};  
defaults.connectivity_window = 1000;              
defaults.connectivity_step = 500;            
defaults.save_preprocessed = true;         
defaults.output_folder = '/Users/medicabg/Documents/Projects/012_CamillaPilotTCDS/002_ScogPreliminary/stat';
defaults.load_format = 'set';
defaults.folder_structure = {'group','subject', 'session'}; 

defaults.channelLocations = [eeglab_path filesep 'plugins/dipfit/',...
    'standard_BESA/standard-10-5-cap385.elp'];
defaults.resample_rate = 300;
defaults.conditionNames = {'rest1','rest2','s1','s2','s3','s4'};
defaults.conditionTimestamps = '';
defaults.conditionDurations = {72, 72, 72, 72, 72, 72};
defaults.epoch_type = 'block';
defaults.trimToEvents = true;
%-- Validating user variables, setting to default if variable not present
user_vars = validateAnalyticParameters(user_vars, defaults);

%% Preprocessing __________________________________________________________
% Internal working variable setup
% -- initialize statistics data carrier
EEGStats = struct('bands',        {},...
                  'power',        {},...
                  'erp',          {},...
                  'connectivity', {},...
                  'timefreq',     {},...
                  'subject',      {},...
                  'stimulus',     {},...
                  'chanLocs',     {});
% -- initialize demographics data carrier
demograph= struct('subject',      {},...
                  'age',          {},...
                  'gender',       {},...
                  'group',        {},...
                  'session',      {},...
                  'condition',    {});

% -- find all appropriate EEG files
% Parameters:
%   user_vars.folder_structure   - folder hierarchy
%   user_vars.load_format        - specify file type to be loaded
file_list = loadEEGPaths(load_path, user_vars);

% Main preprocessing loop
for i = 1:length(file_list)
cnt = ['[',num2str(i),'/',num2str(length(file_list)),'] '];
% -- Data loading
[EEG, logfile] = loadEEGData(file_list(i), logfile);
disp([cnt,'Loading file: ', file_list(i).filepath]);

% -- Channel labeling (and positioning)
if ~isempty(user_vars.channelLabels)
    channels = user_vars.channelLabels;
    [EEG.chanlocs(1:length(channels)).labels] = channels{:};
    EEG = pop_chanedit(EEG, 'lookup', user_vars.channelLocations);
    disp([cnt,'Labelling, assigning topography.']);
end

% -- Re-sampling frequency
if ~isempty(user_vars.resample_rate)
    sf = user_vars.resample_rate;
    if EEG.srate ~= sf, EEG = pop_resample(EEG, sf); end
    disp([cnt,'Re-sampling to: ',num2str(user_vars.resample_rate),'Hz']);
end

% -- Formatting conditions
% >> names
if ~isempty(user_vars.conditionNames)
    names = user_vars.conditionNames;
    ermsg = 'Number of names must match number of conditions';
    if length(names) ~= length(EEG.event), error(ermsg), end
    [EEG.event(1:length(names)).type] = names{:};
end
% >> timestamps
if ~isempty(user_vars.conditionTimestamps)
    latencies = user_vars.conditionTimestamps;
    ermsg = 'Number of timestamps must match number of conditions';
    if length(latencies) ~= length(EEG.event), error(ermsg), end
    [EEG.event(1:length(latencies)).latency] = latencies{:};
end
% >> durations
if ~isempty(user_vars.conditionDurations)
    durations = user_vars.conditionDurations;
    ermsg = 'Number of durations must match number of conditions';
    if length(durations) ~= length(EEG.event), error(ermsg), end
    [EEG.event(1:length(durations)).duration] = durations{:};
end

% -- Trim to events
if user_vars.trimToEvents && ~isempty(EEG.event)
    [st, ed] = deal(EEG.event(1).latency, EEG.event(end).latency);
    if isfield(EEG.event(end),'duration')
        if ~isempty(EEG.event(end).duration)
            temp = ed + EEG.srate * EEG.event(end).duration;
            ed = clip(temp, 0, size(EEG.data,2));
        end
    end
    EEG = pop_select( EEG, 'point', [round(st) round(ed)] );
    EEG = cbg_boundarybreak(EEG);
end

% -- Re-referencing electrodes
if ~isempty(user_vars.reference_channels)
    EEG = pop_reref(EEG, user_vars.reference_channels);
else
    EEG = pop_reref(EEG, [], 'interpchan', []);
end

% -- Frequency filtering
% >> band-pass
if ~isempty(user_vars.bandpass_filter)
    [hp, lp] = deal(user_vars.bandpass_filter{:});
    EEG = pop_eegfiltnew(EEG, 'locutoff', hp, 'hicutoff', lp);
end
% >> notch
if ~isempty(user_vars.notch_filter)
    for f = user_vars.notch_filter
        EEG = pop_eegfiltnew(EEG, f-1, f+1, [], 1);
    end
end
%%    
% -- Independent Component Analysis
if user_vars.run_ica
    % -- pre-ICA high-pass (1Hz)
    % Winkler I, Debener S, Müller KR, Tangermann M.  On the  influence 
    % of high-pass filtering on ICA-based artifact reduction in EEG-ERP
    % Annu  Int  Conf  IEEE   Eng  Med   Biol   Soc.   2015;2015:4101-5
    % doi:          10.1109/EMBC.2015.7319296.          PMID: 26737196.
    EEG_for_ica = pop_eegfiltnew(EEG, 1, []);
    
    % -- perform ICA
    %data_rank = rank(reshape(EEG_for_ica.data, EEG_for_ica.nbchan, []));
    % data_rank = EEG_for_ica.nbchan;
    EEG_for_ica = pop_runica(EEG_for_ica, 'icatype', 'runica', ...
                             'extended', 1);
    %EEG_for_ica = pop_runica(EEG_for_ica, 'icatype', 'runica', ...
    %                        'extended', 1, 'options', {'pca', data_rank});
    
    % -- copy ICA weights to original dataset
    EEG.icaweights  = EEG_for_ica.icaweights;
    EEG.icasphere   = EEG_for_ica.icasphere;
    EEG.icawinv     = EEG_for_ica.icawinv;
    EEG.icachansind = EEG_for_ica.icachansind;

    % -- automatic artifact removal via ICLabel
    EEG_for_ica = pop_iclabel(EEG_for_ica, 'default');
    EEG_for_ica = pop_icflag(EEG_for_ica, user_vars.iclabel_thresholds);
    EEG = pop_subcomp(EEG, find(EEG_for_ica.reject.gcompreject), 0);
end
%%
% -- Block / Epoch separation
% Parameters:
%   user_vars.epoch_type         - 'epoch', 'block' and '' for nothing
%   user_vars.epoch_window       - {start end} in seconds
%   user_vars.baseline_window    - baseline correction window in ms
%   user_vars.artifact_threshold - rejection threshold in µV
EEG = cbg_epochseparation(EEG, user_vars);
    
% -- Extract demographics if available

%% Processing
% After cbg_epochseparation:
if ~isscalar(EEG)
    
    for segIdx = 1:length(EEG)
        demograph = cbg_demographicextraction(EEG(segIdx), demograph, i);
        
        % Run processing modules
        EEGStats(end+1).subject = demograph(end).subject;
        EEGStats(end).stimulus = EEG(segIdx).condition;
        EEGStats(end).chanLocs = EEG(segIdx).chanlocs;
        [EEG(segIdx), EEGStats(end)] = cbg_PSD(EEG(segIdx), EEGStats(end), user_vars);
        [EEG(segIdx), EEGStats(end)] = cbg_ERP(EEG(segIdx), EEGStats(end));
        [EEG(segIdx), EEGStats(end)] = cbg_dFC(EEG(segIdx), EEGStats(end),user_vars);
        [EEG(segIdx), EEGStats(end)] = cbg_TFA(EEG(segIdx), EEGStats(end),user_vars);
        EEG(segIdx).etc.eegstats = EEGStats(end);
    end
else
    demograph = cbg_demographicextraction(EEG, demograph, i);
    EEGStats(end+1).subject = demograph(end).subject;
    EEGStats(end).stimulus = EEG.condition;
    EEGStats(end).chanLocs = EEG.chanlocs;
    [EEG, EEGStats(end)] = cbg_PSD(EEG, EEGStats(end), user_vars);
    [EEG, EEGStats(end)] = cbg_ERP(EEG, EEGStats(end));
    [EEG, EEGStats(end)] = cbg_dFC(EEG, EEGStats(end), user_vars);
    [EEG, EEGStats(end)] = cbg_TFA(EEG, EEGStats(end), user_vars);
    EEG.etc.eegstats = EEGStats(end);
end

% Save preprocessed data if requested
if user_vars.save_preprocessed
% -- setup output path
sdir = user_vars.output_folder; if ~isfolder(sdir); mkdir(sdir); end
% ..iterate through EEG segments    
for segIdx = 1:length(EEG)
    % -- dynamic working variables
    fn = ''; fv = {'group','subject','session','condition'};
    % -- set name, save file
    for o = 1:length(fv), fn = [fn EEG(segIdx).(fv{o}) '_']; end
    filename = fullfile(sdir, [fn(1:end-1) '.set']);
    pop_saveset(EEG(segIdx), 'filename', filename, 'savemode', 'onefile');
end
end
end
% Write demographic data
% demograph = struct2table(demograph);
% if user_vars.save_preprocessed, writetable(demograph, 'test.csv'), end
% % Display completion
% disp('EEG processing completed successfully!');
end

%% Processing functions ___________________________________________________
% Epoch / Block separation ________________________________________________
function EEG = cbg_epochseparation(EEG, user_vars)
global logfile
% -- only epoch if markers exist
logfile.(EEG.subject).epoch.operation = 'none';
logfile.(EEG.subject).epoch.output = 'No events in file.';
if isempty(EEG.event), return, end

% -- type of epoching
% >> no epoch
logfile.(EEG.subject).epoch.output = 'No epoching specified.';
if isempty(user_vars.epoch_type), return, end
% >> 'block' or 'epoch' (default: 'block')
logfile.(EEG.subject).epoch.operation = 'block';
do_epoch = false;
switch user_vars.epoch_type
    case 'epoch'
        do_epoch = true;
        logfile.(EEG.subject).epoch.operation = 'epoch';
    case 'block'
        do_epoch = false;
end

% -- Epoch separation
if do_epoch
% >> static working variable definition
event_types = unique({EEG.event.type});
[wl, wh] = deal(user_vars.epoch_window{:});

% >> epoch separation, baseline correction, artifact rejection
EEG = pop_epoch(EEG, event_types, [wl wh]);
EEG = pop_rmbase(EEG, user_vars.baseline_window);

% >> artifact rejection based on amplitude threshold
if isempty(user_vars.artifact_threshold), return, end
thresh = user_vars.artifact_threshold;
[EEG, rejected_epochs] = ...
    pop_eegthresh(EEG, 1, 1:EEG.nbchan, -thresh, thresh, wl, wh, 0, 1);
logfile.(EEG.subject).epoch.output = rejected_epochs;
EEG.condition = 'entire';

% -- Block separation
else
% >> static working variable definition
EEG.event(end+1).type = 'end marker';
EEG.event(end).position = 10;
EEG.event(end).latency = size(EEG.data,2);
EEG.event(end).urevent = length(EEG.event);
% .. iterate through events
for e = 1:length(EEG.event) - 1
% >> block duration
[st, ed] = deal(EEG.event(e).latency, EEG.event(e + 1).latency);
if isfield(EEG.event,'duration') && ~isempty(EEG.event(e).duration)
    ed = clip(st + EEG.srate * EEG.event(e).duration, 0, size(EEG.data,2));
end
temp = pop_select( EEG, 'point', [round(st) round(ed)] );
temp = cbg_boundarybreak(temp);
EEG_block(e) = temp;
EEG_block(e).condition = EEG.event(e).type;
end
logfile.(EEG.subject).epoch.output =['Created ',num2str(e),' blocks.'];
% >> assign return variables
EEG = EEG_block;
end
end
% _________________________________________________________________________
% Demographic extraction
function demograph = cbg_demographicextraction(EEG, demograph, enum)
    fields = {'subject', 'age', 'gender', 'group', 'session'};
    mask = isfield(EEG, fields);
    for f = 1:length(mask)
        if mask(f), demograph(enum).(fields{f}) = EEG.(fields{f}); end
    end
    
    % Handle event structures by flattening them
    if isfield(EEG, 'event') && ~isempty(EEG.event)
        % Get all field names from the event structure
        eventFields = fieldnames(EEG.event);
        
        % For each event field, create a cell array of values
        for i = 1:length(eventFields)
            fieldName = eventFields{i};
            fieldValues = {EEG.event.(fieldName)};
            
            % Convert numerical arrays to strings for consistent handling
            if isnumeric(fieldValues{1})
                fieldValues = cellfun(@num2str, fieldValues, ...
                    'UniformOutput', false);
            end
            
            % Join values with a delimiter for storage
            demograph(enum).(['event_' fieldName]) = ...
                strjoin(fieldValues, '|');
        end
        
        % Store the count of events
        demograph(enum).event_count = length(EEG.event);
    else
        % No events found
        demograph(enum).event_count = 0;
    end
end

function [EEG, EEGStats] = cbg_PSD(EEG, EEGStats, user_vars)
% -- spectral power (PSD) for each frequency band
% >> static working variable definition
band_names = fieldnames(user_vars.frequency_bands);
EEGStats.bands = struct();
EEGStats.power = struct();
% .. iterating through pre-defined bands
for b = 1:length(band_names)
    % >> dynamic working variable definition
    band_name = band_names{b};
    band_range = user_vars.frequency_bands.(band_name);
    
    % >> filter to the band of interest
    EEG_band = pop_eegfiltnew(EEG, band_range(1), band_range(2));
    EEGStats.bands.(band_name) = EEG_band;
    
    % >> calculate PSD
    [power, ~] = spectopo(EEG_band.data, 0, EEG_band.srate, 'plot', 'off');
    EEGStats.power.(band_name) = power;
end
end

function [EEG, EEGStats] = cbg_ERP(EEG, EEGStats)
% -- event related potential (ERP) for epoched data
if EEG.trials > 1
    % >> static working variable definition
    event_types = unique({EEG.event.type});
    EEGStats.erp = struct();
    % .. iterating through event types
    for e = 1:length(event_types)
        % >> dynamic working variable definition
        event_type = event_types{e};
        event_epochs = find(strcmp({EEG.epoch.eventtype}, event_type));
        
        % >> calculate ERP
        if ~isempty(event_epochs)
            erp_data = mean(EEG.data(:,:,event_epochs), 3);
            EEGStats.erp.(sanitizeFieldName(event_type)) = erp_data;
        end
    end
end
end

function [EEG, EEGStats] = cbg_dFC(EEG, EEGStats, user_vars)
% -- connectivity (dFC) on-demand
if ~isempty(user_vars.connectivity_bands) && ...
        strcmpi(user_vars.connectivity_metric, 'rt_DCCA')
    % >> static working variable definition
    EEGStats.connectivity = struct();
    
    % .. iterating through pre-defined bands
    for b = 1:length(user_vars.connectivity_bands)
        band_name = user_vars.connectivity_bands{b};
        
        if isfield(EEGStats.bands, band_name)
            data = EEGStats.bands.(band_name).data;
            
            % Calculate number of samples for window and step based on ms
            ssr = EEG.srate / 1000;
            window_samples = round(user_vars.connectivity_window * ssr);
            step_samples = round(user_vars.connectivity_step * ssr);
            scales = ...
                round([window_samples/4,window_samples/2, window_samples]);
            num_channels = EEG.nbchan;
            
            if EEG.trials == 1
                % Process all channels simultaneously
                [dccc_matrix, ~, ~] = rt_DCCA(data, scales, ...
                    'WindowSize', window_samples, ...
                    'StepSize', step_samples, ...
                    'BandFrequency', band_name, ...
                    'SamplingRate', EEG.srate, ...
                    'OutputMode', 'final');
                
                % Extract the DCCC values for the largest scale
                conn_matrix = dccc_matrix;
            else
                epoch_conn_matrices = zeros(num_channels, num_channels, ...
                    length(scales), EEG.trials);           
                for trial = 1:EEG.trials
                    trial_data = squeeze(data(:, :, trial));
                    
                    % Process all channels simultaneously for this trial
                    [dccc_matrix, ~, ~] = rt_DCCA(trial_data, scales, ...
                        'WindowSize', window_samples, ...
                        'StepSize', step_samples, ...
                        'BandFrequency', 'broadband', ...
                        'SamplingRate', EEG.srate, ...
                        'OutputMode', 'final');
                    
                    % Store the DCCC matrix for this trial
                    epoch_conn_matrices(:, :, :, trial) = dccc_matrix;
                end
                conn_matrix = mean(epoch_conn_matrices, 4);
            end
            EEGStats.connectivity.(band_name).scales = scales;
            EEGStats.connectivity.(band_name).dccc = conn_matrix;
        end
    end
end
end

function [EEG, EEGStats] = cbg_TFA(EEG, EEGStats, user_vars)
% -- time-frequency analysis (TFA) for epoched data
if EEG.trials > 1
% >> static working variable definition
event_types = unique({EEG.event.type});
EEGStats.timefreq = struct();
% .. iterate through event types
for e = 1:length(event_types)
% >> dynamic working variable definition
event_type = event_types{e};
event_epochs = find(strcmp({EEG.epoch.eventtype}, event_type));
    if ~isempty(event_epochs)
    % Calculate time-frequency transform for each channel
    tf_data = cell(EEG.nbchan, 1);
    for ch = 1:EEG.nbchan
    % Extract channel data for specified epochs
    channel_data = squeeze(EEG.data(ch, :, event_epochs));
    
    % Perform time-frequency analysis using EEGLAB's newtimef
    [wl, wh] = deal(user_vars.epoch_window{:});
    [tf_power, ~, ~, times, freqs] = ...
        newtimef(channel_data, EEG.pnts, [wl*1000, wh*1000], ...
        EEG.srate, 'cycles', [3 0.5], 'baseline', ...
        user_vars.baseline_window, 'plotitc', 'off', 'plotersp', 'off');
    tf_data{ch} = tf_power;
    end
    EEGStats.timefreq.(sanitizeFieldName(event_type)).power = tf_data;
    EEGStats.timefreq.(sanitizeFieldName(event_type)).times = times;
    EEGStats.timefreq.(sanitizeFieldName(event_type)).freqs = freqs;
    end
end
end
end

%% Auxiliary functions ____________________________________________________
% List all EEG files with metadata
function file_list = loadEEGPaths(load_path, user_vars)
    folder_structure = user_vars.folder_structure;
    file_ext = user_vars.load_format;
    file_list = struct('filepath', {}, 'metadata', {});
    
    % First check if the input is a file
    if isfile(load_path)
        [~, ~, ext] = fileparts(load_path);
        if contains(ext, file_ext, 'IgnoreCase', true)
            file_list(1).filepath = load_path;
            file_list(1).metadata = struct(); % Empty metadata structure
            return; % Return early if it's a file
        else
            warning('Input file does not match the specified format: %s', file_ext);
            return;
        end
    end
    
    % If we got here, load_path is a directory
    % Get all matching files recursively
    files = dir(fullfile(load_path, '**/*.*'));  % Recursive search
    files = files(~[files.isdir]);
    
    % Filter by extensions
    if ischar(file_ext)
        ext_mask = contains({files.name}, file_ext, 'IgnoreCase', true);
    elseif iscell(file_ext)
        ext_mask = false(size(files));
        for i = 1:length(file_ext)
            ext_mask = ext_mask | contains({files.name}, file_ext{i}, 'IgnoreCase', true);
        end
    else
        error('file_ext must be a string or cell array of strings');
    end
    files = files(ext_mask);
    
    % Process each file
    for i = 1:length(files)
        full_path = fullfile(files(i).folder, files(i).name);
        
        % Extract metadata from folder structure
        rel_path = strrep(files(i).folder, load_path, '');
        if startsWith(rel_path, filesep)
            rel_path = rel_path(2:end);
        end
        
        path_parts = strsplit(rel_path, filesep);
        metadata = struct();
        
        % Assign folder structure values
        for j = 1:min(length(path_parts), length(folder_structure))
            field_name = folder_structure{j};
            metadata.(field_name) = path_parts{j};
        end
        
        % Add to file list
        file_list(end+1) = struct(...
            'filepath', full_path, ...
            'metadata', metadata ...
        );
    end
end
% Load in EEG data
function [EEG, logfile] = loadEEGData(pathmeta, logfile)
    [folderPath, fileName, fileExt] = fileparts(pathmeta.filepath);
    if isempty(fileExt)
        error('No file extension found in path: %s', pathmeta.filepath);
    end
    
    try
        switch fileExt
            case '.set'
                EEG = pop_loadset('filename',[fileName fileExt],...
                                  'filepath',folderPath);
            case '.hdf5'
                EEG = pop_loadhdf5('filename', [fileName fileExt], ...
                                   'filepath', folderPath, ...
                                   'rejectchans', []);
            case '.edf'
                EEG = pop_biosig(pathmeta.filepath);
            otherwise
                error('Unsupported format: %s.', fileExt);
        end
        % -- adding metadata fields (assumption: id in folder/filename)
        if isfield(pathmeta, 'metadata')
            valid_fields = {'subject','group','condition','session'};
            meta_fields = fieldnames(pathmeta.metadata);
            for f = 1:length(meta_fields)
                if ~ismember(meta_fields{f}, valid_fields), continue, end
                EEG.(meta_fields{f}) = pathmeta.metadata.(meta_fields{f});
            end
        end
        temp_name = split(matlab.lang.makeValidName(fileName),'_');
        if isempty(EEG.group),     EEG.group     = 'group';        end
        if isempty(EEG.subject),   EEG.subject   = temp_name{1};   end
        if isempty(EEG.session),   EEG.session   = 'session';      end
        if isempty(EEG.condition), EEG.condition = 'condition';    end
        % -- log on success
        logfile.(EEG.subject).load.operation = fileExt;
        logfile.(EEG.subject).load.output = 'File loaded.';
    catch e
        % -- log on failure
        logfile.(fileName).load.operation = fileExt;
        logfile.(fileName).load.output = e.message;
    end
end

function params = validateAnalyticParameters(params, defaults)
% Validate and set default parameters if not provided
    if ~exist('params', 'var') || ~isstruct(params)
        params = struct();
    end
    
    field_names = fieldnames(defaults);
    for i = 1:length(field_names)
        field = field_names{i};
        if ~isfield(params, field) || isempty(params.(field))
            params.(field) = defaults.(field);
        end
    end
    
    % These fields should remain as cell arrays
    fields_requiring_cells = {
        'folder_structure',
        'connectivity_bands',
        'bandpass_filter',
        'epoch_window',
        'conditionDurations'
    };
    
    for i = 1:length(fields_requiring_cells)
        field = fields_requiring_cells{i};
        if isfield(params, field) && ~isempty(params.(field))
            % Check if the field is a numeric array and not already a cell
            if isnumeric(params.(field)) && ~iscell(params.(field))
                params.(field) = num2cell(params.(field));
                disp(['Converted numeric array to cell array for: ', field]);
            end
            % Ensure it's an appropriate cell array (not nested cells)
            if iscell(params.(field)) && ~isempty(params.(field)) && iscell(params.(field){1})
                % If we got a cell array containing a cell array (like {{1,2}}), 
                % take the inner cell
                params.(field) = params.(field){1};
            end
        end
    end
    
    % Handle the bandpass_filter special case
    if isfield(params, 'bandpass_filter')
        if iscell(params.bandpass_filter) && numel(params.bandpass_filter) > 2
            % If we got more than 2 elements, take just the first two
            params.bandpass_filter = {params.bandpass_filter{1}, params.bandpass_filter{2}};
        elseif isnumeric(params.bandpass_filter) && numel(params.bandpass_filter) == 2
            params.bandpass_filter = {params.bandpass_filter(1), params.bandpass_filter(2)};
        end
    end
    
    % Handle the epoch_window special case
    if isfield(params, 'epoch_window')
        if iscell(params.epoch_window) && numel(params.epoch_window) > 2
            % If we got more than 2 elements, take just the first two
            params.epoch_window = {params.epoch_window{1}, params.epoch_window{2}};
        elseif isnumeric(params.epoch_window) && numel(params.epoch_window) == 2
            params.epoch_window = {params.epoch_window(1), params.epoch_window(2)};
        end
    end
    
    % Handle the baseline_window special case
    if isfield(params, 'baseline_window')
        if iscell(params.baseline_window)
            % Convert from cell to numeric array since baseline_window is used as numeric
            if numel(params.baseline_window) == 2
                params.baseline_window = [params.baseline_window{1}, params.baseline_window{2}];
            elseif numel(params.baseline_window) == 1 && isnumeric(params.baseline_window{1}) && numel(params.baseline_window{1}) == 2
                % Handle case where cell contains an array: {[x y]}
                params.baseline_window = params.baseline_window{1};
            end
        end
        % Ensure it's a row vector
        if isnumeric(params.baseline_window) && numel(params.baseline_window) == 2
            params.baseline_window = reshape(params.baseline_window, 1, 2);
        end
    end
    
    % Handle the frequency_bands special case
    if isfield(params, 'frequency_bands')
        % frequency_bands is a struct with band ranges
        band_names = fieldnames(params.frequency_bands);
        for b = 1:length(band_names)
            band = band_names{b};
            band_value = params.frequency_bands.(band);
            
            % If the band range is a cell array, convert to numeric
            if iscell(band_value)
                if numel(band_value) == 2
                    % Cell array like {x, y}
                    params.frequency_bands.(band) = [band_value{1}, band_value{2}];
                elseif numel(band_value) == 1 && isnumeric(band_value{1}) && numel(band_value{1}) == 2
                    % Cell containing array: {[x y]}
                    params.frequency_bands.(band) = band_value{1};
                end
            end
            
            % Ensure it's a row vector
            if isnumeric(params.frequency_bands.(band)) && numel(params.frequency_bands.(band)) == 2
                params.frequency_bands.(band) = reshape(params.frequency_bands.(band), 1, 2);
            end
        end
    end
    
    % Additional checks for specific fields
    if isfield(params, 'folder_structure') && ~iscell(params.folder_structure)
        warning('folder_structure should be a cell array. Converting...');
        params.folder_structure = {params.folder_structure};
    end
    
    if isfield(params, 'connectivity_bands') && ~iscell(params.connectivity_bands)
        warning('connectivity_bands should be a cell array. Converting...');
        params.connectivity_bands = {params.connectivity_bands};
    end
end

function EEG = cbg_boundarybreak(EEG)
    non_boundary_idx = find(~strcmp({EEG.event.type}, 'boundary'));
    EEG.event = EEG.event(non_boundary_idx);
    EEG = eeg_checkset(EEG, 'eventconsistency');
end

function field_name = sanitizeFieldName(name)
% Convert a string to a valid MATLAB field name
    % Remove spaces and special characters
    field_name = regexprep(name, '[^a-zA-Z0-9]', '_');
    
    % Ensure it starts with a letter
    if ~isletter(field_name(1))
        field_name = ['x', field_name];
    end
end