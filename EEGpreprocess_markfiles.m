%% Dependencies
% Written in MATLAB R2022a (2024).
% Author: ZK

%% User variables
% USER VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root_dirx = 'G:\Projects\OUHSCgl\prcEEG\';
file_load = 'sample_data';
fedf_save = 'mrk_data';
samp_freq = 128;
mark_chan = 23;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Environment variables, user error handling
% Addition of files to workspace, initial user error handling and file
% management ______________________________________________________________
if isempty(root_dirx)
    root_dirx = [cd, filesep];
end
file_load = [root_dirx, file_load,filesep];
fedf_save = [root_dirx, fedf_save,filesep];

fileList = dir(fullfile(file_load,'**','*.edf'));
global debg_log
debg_log = "";

for i=1:length(fileList)
    % Identify movement data files and exclude them _______________________
    if contains(fileList(i).name,{'.md.','.mdedf'})
        continue
    end % _________________________________________________________________
   
    % Load in data file, catch markers ____________________________________
    fname_in  = [fileList(i).folder, filesep ,fileList(i).name];
    temp_a = dir(fullfile(fileList(i).folder,'*SubjectTimeLog.mat'));
    fmrkr_in  = [fileList(i).folder, filesep ,temp_a.name];
    clear temp_a
    disp(['IMPRT: ',fileList(i).name,' and started extraction.'])
    try
        [hdr,rec] = edfread(fname_in);
    catch
        err = "ERROR: EDF not readable.";
        debg_log(end+1) = err;
        disp(err)
    end

    fname_ex = fullfile(...
               append(fedf_save,...
               strrep(fileList(i).folder,file_load,''),...
               filesep));
    if ~exist(fname_ex,'dir')
        mkdir(fname_ex);
    end

    srq = [1,3,25];
    mrk = find(rec(mark_chan,:));
    val = rec(mark_chan,mrk);
    
    switch length(mrk)
        case 0
            err = string(['[ ',fileList(i).name,' ] ERROR: Includes',...
                num2str(length(mrk)),' markers.']);
            debg_log(end+1) = err;
            disp(err)
            continue
        case 1
            try
                load(fmrkr_in)
            catch
                err = string(['[ ',fileList(i).name,' ] ERROR: Includes',...
                    ' no marker file.']);
                debg_log(end+1) = err;
                disp(err)
            end
            [x,y] = find(SubjectTimeLog == "Instructions");
            mrk(1:4) = (str2double(SubjectTimeLog(x,1))*samp_freq + mrk(1));
            val(1:4) = 25;
        case 4
            err = string(['[ ',fileList(i).name,' ] ERROR: Includes',...
                num2str(length(mrk)),' markers.']);
            debg_log(end+1) = err;
            disp(err)
            continue
        case 5
            mrk      = mrk([1,3,5]);
            val      = val([1,3,5]);
            try
                load(fmrkr_in)
            catch
                err = string(['[ ',fileList(i).name,' ] ERROR: Includes',...
                    ' no marker file.']);
                debg_log(end+1) = err;
                disp(err)
            end
            [x,y] = find(SubjectTimeLog == "Instructions");
            mrk(3:6) = str2double(SubjectTimeLog(x,1))*samp_freq + mrk(3);
            val(3:6) = 25;
        otherwise
            err = string(['[ ',fileList(i).name,' ] ERROR: Includes',...
                num2str(length(mrk)),' markers.']);
            debg_log(end+1) = err;
            disp(err)
            continue
    end
    
    % Cycle through markers _______________________________________________
    mrkr_data = [];
    for p=1:length(mrk)
        switch val(p)
            case 1
                mrkr_data(p,:) = [mrk(p) mrk(p)+16*samp_freq mrk(p)+16*samp_freq];
            case 3
                mrkr_data(p,:) = [mrk(p) mrk(p)+16*samp_freq mrk(p)+16*samp_freq];
            case 25
                mrkr_data(p,:) = [mrk(p) mrk(p)+73*samp_freq mrk(p)+73*samp_freq];
            otherwise
                err = string(['[ ',fileList(i).name,' ] ERROR: Includes ',...
                    'unexpected marker value of ',num2str(length(mrk))]);
                debg_log(end+1) = err;
                disp(err)
        end
    end % _________________________________________________________________

    % Save data, marker ___________________________________________________
    fname_in = fullfile(fileList(i).folder, fileList(i).name);
    [~, lfn, ~] = fileparts(fileList(i).folder);
    fname_tm = ['subjectid_',lfn,'_mrk'];
    fname_ex = fullfile(...
                    append(fedf_save,...
                           strrep(fileList(i).folder,file_load,''),...
                           filesep));
    if ~exist(fname_ex,'dir')
        mkdir(fname_ex);
    end
    fname_ex = [fname_ex,fname_tm];
    if ~isfile([fname_ex,'.edf'])
        copyfile(fname_in, [fname_ex,'.edf']);
        save([fname_ex,'.mat'],"mrkr_data");
    else
        copyfile(fname_in, [fname_ex,num2str(i),'.edf'])
        save([fname_ex,num2str(i),'.mat'],"mrkr_data");
    end
    disp(['SAVED: ',fname_ex,'.edf with markers:'])
    disp(mrkr_data/samp_freq)
    % _____________________________________________________________________
end