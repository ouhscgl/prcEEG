
%% Dependencies
% > edfread()
% > eeglab()
% Written in MATLAB R2022a (2024).
% Author: ZK

%% User variables
% > root_dirx: Root directory, defaults to current directory if non given.
% Keep in mind, that in that case default save directories will be used.
% You can disable this behaviour by inputting a [space] as root.
% Example: root_dirx = 'C:\Users\user\test\' OR '' OR ' '

% > eegl_locs: Location of the eeglab toolbox, required dependency.
% Example: eeg_locs = 'C:\Users\user\test\eeglab2022.0'

% > file_load: Main sub-directory where the input files are located,
% required input. Defaults to searching current directory.
% Example: 'C:\Users\user\test\data\' OR [root_dirx, 'data',filesep];

% > fedf_save, fset_save, fmat_save: Save locations for the data. If folder
% doesn't exist, it will be created in the specified locations. If left
% empty, data in that format will not be saved.
% Example: 'C:\Users\user\test\edfs\' OR [root_dirx, 'edfs',filesep] OR ''

% > adds_appx: Appendix to add to the segmented data. Numbering is 2-digit.
% Example: '_segment_' = subject_01.edf >>> subject_01_segment_01.edf

% > mrkr_chan: Index of marker channel.
% Example: 21

% > samp_freq: Sampling frequency of the dataset.
% Example: 256

% > task_lens: Define in case of requiring specific lengths of data.
% Defaults to extracting between markers.
% Example: 72 OR ''

% > user_chan: Defined labels of channel data that will be extracted.
% Defaults to all.
% Example: {'AF3','F7','F3','FC5','T7','P7','O1','O2','P8','T8','FC6',...
%           'F4','F8','AF4'};

% USER VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root_dirx = '';
eegl_locs = 'C:\Users\user\test\eeglab2022.0';
edfr_locs = '';
file_load = [root_dirx, 'sample_data',filesep];
fedf_save = [root_dirx, 'edf_data',filesep];
fset_save = [root_dirx, 'set_data',filesep];
fmat_save = [root_dirx, 'mat_data',filesep];
adds_appx = '_segment_';
mrkr_chan = 21;
samp_freq = 256;
task_lens = 30;
user_chan = {'AF3','F7','F3','FC5','T7','P7','O1','O2','P8','T8','FC6',...
             'F4','F8','AF4'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Environment variables, user error handling
% Addition of files to workspace, initial user error handling and file
% management ______________________________________________________________
if isempty(root_dirx)
    root_dirx = [cd, filesep];
    file_load = [root_dirx, 'sel_data',filesep];
    fedf_save = [root_dirx, 'edf_data',filesep];
    fset_save = [root_dirx, 'set_data',filesep];
    fmat_save = [root_dirx, 'mat_data',filesep];
end
fileList = dir(fullfile(file_load,'**','*.edf'));
if isempty(fileList)
    ErrorLog(1, fullfile(file_load,'**','*.edf'));
    return
end
% Error handle sampling frequency
if isempty(samp_freq) || ~isa(samp_freq,'double')
    ErrorLog(7, 'samp_freq');
    return
end
% Error handle task length
if ~isempty(task_lens) && ~isa(task_lens,'double')
    ErrorLog(7, 'task_lens');
    return
end
% Error handle extract channels
if ~isempty(user_chan) && ~isa(user_chan,'cell')
    ErrorLog(7, 'user_chan');
    return
end
% Error handle file paths
if ~isa(root_dirx,'char') || ~isa(eegl_locs,'char') || ...
   ~isa(file_load,'char') || ~isa(fedf_save,'char') || ...
   ~isa(fset_save,'char') || ~isa(fmat_save,'char')
    ErrorLog(7, 'File paths.');
end %______________________________________________________________________

% Add paths of dependecies to workspace ___________________________________
if isempty(eegl_locs)
    ErrorLog(1, 'eegl_locs');
    return
end
addpath(eegl_locs)
if ~isempty(eegl_locs)
    addpath(edfr_locs)
end %______________________________________________________________________

% Add debug logs __________________________________________________________
global debg_logs, global load_logs, global save_logs
[debg_logs, load_logs, save_logs] = deal("","","");
%__________________________________________________________________________

% Initialize EEGLab _______________________________________________________
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% _________________________________________________________________________

%% Segmenting EEG data
% File list loop >>> Segment data
for i=1:length(fileList)

    % Load in data file, extract markers __________________________________
    fname_in  = [fileList(i).folder, filesep ,fileList(i).name];
    load_logs = [load_logs; string(fname_in)];
    try
        [hdr,rec] = edfread(fname_in);
    catch
         ErrorLog(5, fileList(i).name);
         continue
    end
    
    if mrkr_chan > size(rec,1)
        ErrorLog(6, fileList(i).name);
        continue
    end
    mrk = find(rec(mrkr_chan,:));
    % _____________________________________________________________________
    
    % Identify the number of markers contained in the dataset [0 or X] ____
    switch size(mrk,1)
        case 0
            ErrorLog(1, fileList(i).name);
            scd = cell(1,1);
            seq = zeros(1,3);
        otherwise
            scd = cell(length(mrk),1);
            seq = zeros(length(mrk)-1,3);
    end % _________________________________________________________________
    
    % Identify the mode of segmentation [lengthwise or between markers] ___
    ile = ~isempty(task_lens);
    % _____________________________________________________________________

    % Extract sequences ___________________________________________________
    switch length(mrk)
        case 0
            ErrorLog(2, fileList(i).name);
            seq    = [1, size(rec,2), 1+task_lens*samp_freq];
            scd{1} = rec(:,1:seq(2+ile));
        case 1
            ErrorLog(3, fileList(i).name);
            seq    = [mrk(1), size(rec,2), mrk(1)+task_lens*samp_freq];
            scd{1} = rec(:,mrk(1):seq(2+ile));
        otherwise
            for o=1:length(mrk)-1+ile
                seq(o,:) = [mrk(o), mrk(min(o+1,length(mrk))), ...
                            mrk(o)+task_lens*samp_freq];
                
                % ErrorHandle: Data len too low for fixed length segments
                if ile == 1 && size(rec,2) < seq(o,2+ile)
                    ErrorLog(4,[fileList(i).name,' segment: ',num2str(o)]);
                    seq(o,2+ile) = size(rec,2);
                end
                
                scd{o} = rec(:,seq(o,1):seq(o,2+ile));
            end
    end % _________________________________________________________________
    
    % Single file loop >>> Save data
    for o=1:length(scd)
        
    % Save .mat file if save path is given ________________________________
    if ~isempty(fmat_save)
        slc = append(fmat_save, ...
          strrep(fileList(i).folder, file_load, ''),filesep);
        sfn=[fileList(i).name(1:end-4),adds_appx,sprintf('%02d',o),'.mat'];
        if ~exist(slc,'dir')
            mkdir(slc)
        end
        segm_data = scd{o};
        
        if ~isempty(user_chan)
            idx = [];
            for p = 1:numel(hdr.label)
            if any(strcmp(hdr.label{p}, user_chan))
                idx = [idx, p];
            end
            end
            segm_data = segm_data(idx,:);
        end

        save([slc,sfn],'segm_data')
        save_logs = [save_logs; string([slc,sfn])];
    end % _________________________________________________________________
    
    % Run EEGLab if either .set or .edf data is required as output ________
    if ~isempty(fset_save) || ~isempty(fedf_save)
        EEG = pop_biosig(fname_in,'blockrange',...
              [floor(seq(o,1)/samp_freq) floor(seq(o,2+ile)/samp_freq)]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
        if ~isempty(user_chan)
            EEG = pop_select( EEG, 'channel',user_chan);
        end
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off');
        EEG = eeg_checkset( EEG );
        
    % Save .set file if save path is given ________________________________
    if ~isempty(fset_save)
        slc = append(fset_save, ...
              strrep(fileList(i).folder, file_load, ''),filesep);
        sfn=[fileList(i).name(1:end-4),adds_appx,sprintf('%02d',o) '.set'];
        if ~exist(slc,'dir')
            mkdir(slc)
        end
        EEG = pop_saveset( EEG, 'filename',sfn,'filepath',slc);
        save_logs = [save_logs; string([slc,sfn])];
    end % _________________________________________________________________
        
    % Save .edf file if save path is given ________________________________
    if ~isempty(fedf_save)           
        slc = append(fedf_save, ...
              strrep(fileList(i).folder, file_load, ''),filesep);
        sfn=[fileList(i).name(1:end-4),adds_appx,sprintf('%02d',o) '.edf'];
        if ~exist(slc,'dir')
            mkdir(slc)
        end
        pop_writeeeg(EEG,[slc,sfn],'TYPE','EDF')
        save_logs = [save_logs; string([slc,sfn])];
    end % _________________________________________________________________
    end % _________________________________________________________________
    end
end
% Finished
beep
disp('Files have finished processing.')
%% Error handling function - unclutterring main script
function ErrorLog(errno, instc)
    global debg_logs
    error_list = {'ERROR: No files exist in specified location that match criteria.',...
                  'ERROR: No markers exist in the specified file. Logged, continuing as fullfile.',...
                  'ERROR: Only a single marker exists in the specified file. Logged, continuing as fullfile.',...
                  'ERROR: Segment length smaller than expected. Logged, continuing with remainder of the file.',...
                  'ERROR: EDF file could not be loaded.',...
                  'ERROR: User specified marker channel index is out of bounds. Check data dimensions.',...
                  'ERROR: Please provide valid variable value. Refer to user instructions when needed.'
                 };
    error = [error_list{errno},' @ ',char(instc)];
    disp(error)
    debg_logs = [debg_logs;string(error)];
    return 
end