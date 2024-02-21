
%% Dependencies
% > edfread() - mathworks.com/matlabcentral/fileexchange/31900-edfread (@02/18/24)
% > eeglab()  - sccn.ucsd.edu/eeglab/downloadtoolbox.php (@02/18/24)
% Written in MATLAB R2022a (2024).
% Author: ZK

%% User variables
% > root_dirx: Root directory, defaults to current directory if non given.
% Keep in mind, that in that case default save directories will be used.
% You can disable this behaviour by inputting a [space] as root.
% Example: root_dirx = 'C:\Users\user\test\' OR '' OR ' '

% > eegl_locs: Location of the eeglab toolbox, required dependency.
% Example: eeg_locs = 'C:\Users\user\test\eeglab2022.0\'

% > edfr_locs: Location of the edfread script, required dependency.
% Example: eeg_locs = 'C:\Users\user\test\'

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

% > task_lens: Define in case of requiring specific lengths of data.
% Defaults to extracting between markers.
% Example: 72 OR ''

% > user_chan: Defined labels of channel data that will be extracted.
% Defaults to all.
% Example: {'AF3','F7','F3','FC5','T7','P7','O1','O2','P8','T8','FC6',...
%           'F4','F8','AF4'};

% USER VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root_dirx = 'C:\Projects\CID\';
eegl_locs = 'C:\Projects\_extensions\eeglab2023.1\';
edfr_locs = 'C:\Projects\_extensions\';
file_load = 'mrk_data';
fedf_save = 'seg\edf_data';
fset_save = 'seg\set_data';
fmat_save = 'seg\mat_data';
adds_appx = '_seg_';
mrkr_chan = '';
task_lens = '';
user_chan = {'AF3','F7','F3','FC5','T7','P7','O1','O2','P8','T8','FC6',...
             'F4','F8','AF4'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Environment variables, user error handling
% Addition of files to workspace, initial user error handling and file
% management ______________________________________________________________
if isempty(root_dirx)
    root_dirx = [cd, filesep];
end
file_load = [root_dirx, file_load, filesep];
fedf_save = [root_dirx, fedf_save, filesep];
fset_save = [root_dirx, fset_save, filesep];
fmat_save = [root_dirx, fmat_save, filesep];

fileList = dir(fullfile(file_load,'**','*.edf'));
if isempty(fileList)
    ErrorLog(1, fullfile(file_load,'**','*.edf'));
    return
end
% Error handle task length
if ~isempty(task_lens) && ~isa(task_lens,'char')
    ErrorLog(7, 'task_lens');
    return
end
task_lens = str2double(task_lens);
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
if ~isempty(edfr_locs)
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
    fname_in  = [fileList(i).folder, filesep, fileList(i).name];
    if isempty(mrkr_chan)
        temp_a = dir(fullfile([fileList(i).folder, filesep, '*.mat']));
        markr_in  = [fileList(i).folder, filesep, temp_a.name];
        clear temp_a
    else
        mrkr_chan = str2double(mrkr_chan);
    end
    load_logs = [load_logs; string(fname_in)];
    
    try
        [hdr,rec] = edfread(fname_in);
        sfr = hdr.samples(1);
    catch
         ErrorLog(5, fileList(i).name);
         continue
    end
    
    if mrkr_chan > size(rec,1)
        ErrorLog(6, fileList(i).name);
        continue
    end
    
    % Identify the mode of segmentation [lengthwise or between markers] ___
    ile = isempty(task_lens);
    % _____________________________________________________________________

    % Identify the number of markers contained in the dataset [0 or X] ____
    if ~isempty(mrkr_chan)
        mrk = find(rec(mrkr_chan,:));
        switch size(mrk,1)
            case 1
                if isempty(mrk)
                    ErrorLog(2, fileList(i).name);
                    scd = cell(1,1);
                    mrk      = [1,...
                                size(rec,2),...
                                (1+task_lens*sfr)*ile];
                    scd{1} = rec(:,mrk(1):mrk(2+ile));
                else
                    ErrorLog(3, fileList(i).name);
                    scd = cell(1,1);
                    mrk      = [mrk(1),...
                                size(rec,2),...
                                (mrk(1)+task_lens*sfr)*ile];
                    scd{1} = rec(:,mrk(1):mrk(2+ile));
                end
            otherwise
                scd = cell(size(mrk,1),1);
                for p=1:size(mrk,1)
                mrk(p,:) = [mrk(p),...
                            mrk(min(p+1,length(mrk))),...
                            (mrk(p)+task_lens*sfr)*ile];
                scd{p} = rec(:,mrk(1):mrk(2+ile));
                end
        end 
    else
        load(markr_in)
        mrk = mrkr_data;
        scd = cell(size(mrk,1),1);
        for p=1:size(mrk,1)
        scd{p} = rec(:,mrk(1):mrk(2+ile));
        end
    end
    % _____________________________________________________________________
    
    % Single file loop >>> Save data
    for o=1:length(scd)
        
    % Save .mat file if save path is given ________________________________
    if ~isempty(fmat_save)
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
        
        SaveData(fileList(i),file_load,fmat_save,...
            struct('segm_data',segm_data),'.mat',...
            [adds_appx,sprintf('%02d',o)])
    end % _________________________________________________________________
    
    % Run EEGLab if either .set or .edf data is required as output ________
    if ~isempty(fset_save) || ~isempty(fedf_save)
        EEG = pop_biosig(fname_in,'blockrange',...
              [floor(mrk(o,1)/sfr) floor(mrk(o,2+ile)/sfr)]);
        if ~isempty(user_chan)
            EEG = pop_select( EEG, 'channel',user_chan);
        end
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
        EEG = eeg_checkset( EEG );
        
    % Save .set file if save path is given ________________________________
    if ~isempty(fset_save)
        SaveData(fileList(i),file_load,fset_save,...
            EEG,'.set',...
            [adds_appx,sprintf('%02d',o)])
    end % _________________________________________________________________
        
    % Save .edf file if save path is given ________________________________
    if ~isempty(fedf_save)           
        SaveData(fileList(i),file_load,fedf_save,...
            EEG,'.edf',...
            [adds_appx,sprintf('%02d',o)])
    end % _________________________________________________________________
    end % _________________________________________________________________
    end
end
% Finished ________________________________________________________________
beep
disp('Files have finished processing.')
save('user_logs','debg_logs','load_logs','save_logs')
% _________________________________________________________________________
%% Save function - unclutterring main script
function SaveData(fileList,path_main,path_curr,data,type,adds_appx)
    global save_logs         
    slc = append(path_curr, ...
                 strrep(fileList.folder, path_main, ''),filesep);
    sfn=[fileList.name(1:end-4),adds_appx,type];
    if ~exist(slc,'dir')
        mkdir(slc);
    end
    if      type == ".mat"
        save(fullfile(slc,sfn),'-struct','data');
    elseif  type == ".edf"
        pop_writeeeg(data,[slc,sfn],'TYPE','EDF');
    elseif  type == ".set"
        EEG = pop_saveset( data, 'filename',sfn,'filepath',slc);
    end
    save_logs(end+1,1) = string([slc,sfn]);
    disp(['Saved ',slc,sfn,'.'])
end
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