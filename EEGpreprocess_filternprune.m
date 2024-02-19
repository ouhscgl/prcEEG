
%% Dependencies
% > eeglab()  - sccn.ucsd.edu/eeglab/downloadtoolbox.php (@02/18/24)
% > (optional) MARA() - github.com/irenne/MARA           (@02/18/24)
% > (optional) ADJUST() - nitrc.org/projects/adjust/     (@02/18/24)
% > (optional) EEGHEADER_FAUX.mat - github.com/ouhscgl/prcEEG 
% Written in MATLAB R2022a (2024).
% Author: ZK

%% User variables
% > root_dirx: Root directory, defaults to current directory if non given.
% Keep in mind, that in that case default save directories will be used.
% Example: root_dirx = 'C:\Users\user\test\' OR ''

% > eegl_locs: Location of the eeglab toolbox, required dependency.
% Example: eeg_locs = 'C:\Users\user\test\eeglab2022.0\'

% > file_load: Main sub-directory where the input files are located,
% required input. Defaults to searching current directory.
% Example: 'C:\Users\user\test\data\' OR [root_dirx, 'data',filesep];

% > fedf_save, fset_save, fmat_save: Save locations for the data. If folder
% doesn't exist, it will be created in the specified locations. If left
% empty, data in that format will not be saved.
% > comp_save, vlla_save: Save locations for the independent components and
% the unpruned post-ICA data.
% Example: 'C:\Users\user\test\edfs\' OR [root_dirx, 'edfs',filesep] OR ''

% > artx_rejt: Type of component based artifact rejection to utilize. 
% NOTE:ADJUST has not been implemented due to inefficiencies in the source.
% Example: ['ICLABEL','MARA','ADJUST'] OR ['ICLABELMARAADJUST'] OR etc.

% > save_type, load_type: Specify the file format to be loaded / saved in.
% Save accepts all 3 types, load is single key. Keys are case sensitive.
% Example: ['EDF','SET','MAT'] OR ['EDF','MAT'] OR etc.

% > save_vlla, save_comp: Determines whether un-pruned post-ICA data and
% component data will be saved.
% Example: 1 OR 0

% > lowp_filt, higp_filt: Low-pass and high-pass filter values to use with
% EEGLab's built in defailt band-pass FIR filter. Setting values to 0
% disables pass filter.
% Example: lowp_filt = 0.5 ; higp_filt = 45

% > user_chan: Defined labels of channel data that will be extracted.
% Defaults to all.
% Example: {'AF3','F7','F3','FC5','T7','P7','O1','O2','P8','T8','FC6',...
%           'F4','F8','AF4'};

% USER VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root_dirx = '';
eegl_locs = 'G:\Projects\OUHSCgl\_exe\eeglab2023.1\';
chan_locs = 'plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc';
file_load = 'set_data';

comp_save = 'set\pruned_cmp_data';
vlla_save = 'set\unpruneica_data';
fedf_save = 'set\pruned_edf_data';
fset_save = 'set\pruned_set_data';
fmat_save = 'set\pruned_mat_data';

artx_rejt = ['ICLABEL','MARA','ADJUST'];
save_type = ['EDF','SET','MAT'];
load_type = ['SET'];
save_vlla = 1;
save_comp = 1;

lowp_filt = 0.5;
higp_filt = 45;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Environment variables, user error handling
% Addition of files to workspace, initial user error handling and file
% Add paths of dependecies to workspace ___________________________________
if isempty(eegl_locs)
    ErrorLog(1, 'eegl_locs');
    return
end
addpath(eegl_locs)
% _________________________________________________________________________

% Add file access paths to workspace ______________________________________
if isempty(root_dirx)
    root_dirx = [cd, filesep];
end 
file_load = [root_dirx, file_load,filesep];
comp_save = [root_dirx, comp_save,filesep];
vlla_save = [root_dirx, vlla_save,filesep];
fedf_save = [root_dirx, fedf_save,filesep];
fset_save = [root_dirx, fset_save,filesep];
fmat_save = [root_dirx, fmat_save,filesep];
% _________________________________________________________________________

% File loading, error handing _____________________________________________
temp_ptv = ['*.',lower(load_type)];
fileList = dir(fullfile(file_load,'**',temp_ptv));
if isempty(fileList)
    ErrorLog(1, fullfile(file_load,'**',temp_ptv));
    return
end
clear temp_ptv
% _________________________________________________________________________

% Add debug logs __________________________________________________________
global debg_logs, global load_logs, global save_logs
[debg_logs, load_logs, save_logs] = deal("","","");
%__________________________________________________________________________

% Initialize EEGLab _______________________________________________________
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% _________________________________________________________________________

% Identify the type of component based artifact rejection _________________
rem = [contains(artx_rejt,'ICLABEL'),...
       contains(artx_rejt,'MARA'),...
       contains(artx_rejt,'ADJUST')];
% _________________________________________________________________________

% Identify the format(s) of saving data ___________________________________
    sdt = [contains(save_type,'EDF'),...
           contains(save_type,'SET'),...
           contains(save_type,'MAT')];
 % ________________________________________________________________________

% Error handle helper variable for injecting .mat into EEGLab _____________
if rem(2) && load_type == "MAT"
    try
    load(faux_eegs)
    catch
    ErrorLog(1, 'faux_eegs');
    return
    end
end
% _________________________________________________________________________

%% Artifact rejection of EEG data
for i=1:length(fileList)
    fname_in  = [fileList(i).folder, filesep ,fileList(i).name];
    load_logs = [load_logs; string(fname_in)];
    
    % Import data to EEGLab workspace _____________________________________
    switch load_type
        case "EDF"
            EEG = pop_biosig(fname_in);
        case "MAT"
            try
            EEG       = load(faux_eegs);
            EEG.data  = load(fname_in);
            EEG.times = 1:1000/EEG.rate:size(EEG.data,2)*1000/EEG.rate;
            catch
            ErrorLog(1, 'mat file');   
            end
        case "SET"
            EEG = pop_loadset('filename',fileList(i).name,...
                  'filepath',[fileList(i).folder, filesep]);
        otherwise
            ErrorLog(1, 'save_type'); 
            return
    end % _________________________________________________________________

    % Setup EEGLab workspace, specify channel map _________________________
    % Advisable channel montage.
    % Ref.: EEGLab user manual regarding filtering with ICA
    EEG=pop_chanedit(EEG,'lookup',[eegl_locs, filesep,chan_locs]);
    EEG = eeg_checkset( EEG );
    % _____________________________________________________________________

    % STEP #1: Frequency based artifact rejection / band-pass FIR filter /_
    if higp_filt ~= 0
    EEG = pop_eegfiltnew(EEG, 'locutoff', lowp_filt,...                    
        'hicutoff', higp_filt, 'plotfreqz',0);
    % Advisable before ICA, as convenience. Doesn't affect calculations.
    % Ref.: EEGLab user manual regarding re-referencing
    EEG = pop_reref(EEG,[]);                                               
    EEG = eeg_checkset( EEG );
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt',...
        'on','pca',EEG.nbchan);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off',...
    'setname','Post ICA');
    EEG = eeg_checkset( EEG );
    end
    % _____________________________________________________________________

    % STEP #2: Component based artifact rejection _________________________
    % #####################################################################
    % Method1: ICLabel
    if rem(1)
    EEG = pop_iclabel(EEG, 'default');
    cl = EEG.etc.ic_classification.ICLabel.classifications;
    EEG = eeg_checkset( EEG );
    
    % Save [dataset] after ICA and labeling for future reference __________
    if save_vlla
    SaveData(fileList(i),file_load,vlla_save,EEG,'.set','_postica')
    end % _________________________________________________________________
    
    % Save [component] : [labels],[dataset] for future reference __________
    if save_comp
    labl_comp = zeros(EEG.nbchan,1);
    for chn = 1:EEG.nbchan
        [~,tmp] = max(cl(chn,:));
        labl_comp(chn,1) = tmp;
    end
    % MSE of operation is ~5.2e-05
    data_comp = (EEG.icaweights*EEG.icasphere)*EEG.data;
    SaveData(fileList(i),file_load,comp_save,...
        struct('data_comp',data_comp,'labl_comp',labl_comp),...
        '.mat','_iclabel')
    end
    % _____________________________________________________________________
    % NEEDS REVISION FROM INDEX TO PURE
    % Establish rejected component list, remove components ________________
    labl_rica = [];
    for chn = 1:EEG.nbchan
        [~,tmp] = max(cl(chn,:));
        % Select components to keep
        if tmp ~= 1 && tmp ~= 7
        labl_rica = [labl_rica chn];
        end
    end
    if length(labl_rica) == EEG.nbchan
        ErrorLog(9,fname_in)
        continue
    end
    EEG = pop_subcomp( EEG, labl_rica, 0);
    EEG = eeg_checkset( EEG );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off',...
                              'setname','After pruned with ICLabel');
    temp_lbl = zeros(EEG.nbchan,1);
    temp_lbl(labl_rica) = 1;
    labl_rica = temp_lbl;
    clear temp_lbl
    % _____________________________________________________________________

    % Save processed set __________________________________________________
    % Save .mat [iclabel] : [labels], [dataset] for further processing
    if sdt(3)
        data_rica = EEG.data;
        SaveData(fileList(i),file_load,fmat_save,...
            struct('data_rica',data_rica,'labl_rica',labl_rica),...
            '.mat','_iclabel')
    end
    
    % Save .set [iclabel]: [dataset] for future processing ________________
    if sdt(2)
        SaveData(fileList(i),file_load,fset_save,EEG,'.set','_iclabel')
    end % _________________________________________________________________

    % Save .edf [iclabel]: [dataset] for future processing ________________
    if sdt(1)
        SaveData(fileList(i),file_load,fedf_save,EEG,'.edf','_iclabel')
    end % _________________________________________________________________
    
    % USERLOG _____________________________________________________________
    save_logs(end,2)                     = num2str(sum(labl_rica));
    if ~isempty(labl_rica)
    save_logs(end,3:2+EEG.nbchan) = labl_rica;
    end
    disp(['ICLabel: Removed ',num2str(sum(labl_rica)),' components.'])
    % _____________________________________________________________________
    end
    % #####################################################################
    % Method2: MARA
    if rem(2)
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,...
        'retrieve',1,'study',0);
    [ALLEEG EEG CURRENTSET] = processMARA ( ALLEEG,EEG,CURRENTSET );
    labl_mara = EEG.reject.gcompreject;
    EEG = pop_subcomp( EEG, [], 0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off',...
                              'setname','After pruned with MARA');
    % Save processed set __________________________________________________
    % Save .mat [mara] : [labels], [dataset] for further processing
    if sdt(3)
        data_mara = EEG.data;
        SaveData(fileList(i),file_load,fmat_save,...
            struct('data_mara',data_mara,'labl_mara',labl_mara),...
            '.mat','_icamara')
    end
    
    % Save .set [iclabel]: [dataset] for future processing ________________
    if sdt(2)
        SaveData(fileList(i),file_load,fset_save,EEG,'.set','_icamara')
    end % _________________________________________________________________

    % Save .edf [iclabel]: [dataset] for future processing ________________
    if sdt(1)
        SaveData(fileList(i),file_load,fedf_save,EEG,'.edf','_icamara')
    end % _________________________________________________________________
    
    % USERLOG _____________________________________________________________
    save_logs(end,2)                     = num2str(sum(labl_mara));
    if ~isempty(labl_mara)
    save_logs(end,3:2+EEG.nbchan) = labl_mara;
    end
    disp(['ICLabel: Removed ',num2str(sum(labl_mara)),' components.'])
    % _____________________________________________________________________
    end
    % #####################################################################
    % Method3: ADJUST
    if rem(3)
        % Can't be automated, doesn't work.
        ErrorLog(8,' Please select a different artifact rejection method.')
    end
    % #####################################################################
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
                  'ERROR: Please provide valid variable value. Refer to user instructions when needed.',...
                  'WARNING: ADJUST has not been implemented.',...
                  'ERROR: All components rejected. Skipping analysis.'
                 };
    error = [error_list{errno},' @ ',char(instc)];
    disp(error)
    debg_logs = [debg_logs;string(error)];
    return 
end