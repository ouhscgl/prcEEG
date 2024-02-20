%% Dependencies
% Written in MATLAB R2022a (2024).
% Author: ZK

%% User variables
% USER VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root_dirx = '';
file_load = 'raw_data';
fedf_save = 'mrk_data';
samp_freq = 128;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Environment variables, user error handling
% Addition of files to workspace, initial user error handling and file
% management ______________________________________________________________
if isempty(root_dirx)
    root_dirx = [cd, filesep];
    file_load = [root_dirx, file_load,filesep];
    fedf_save = [root_dirx, fedf_save,filesep];
end

fileList = dir(fullfile(file_load,'**','*.edf'));
mrkrList = dir(fullfile(file_load,'**','*SubjectTrialLog.mat'));

for i=1:length(fileList)
    % Identify movement data files and exclude them _______________________
    if contains(fileList(i).name,{'.md.'})
        continue
    end % _________________________________________________________________
   
    % Load in data file, catch markers ____________________________________
    fname_in  = [fileList(i).folder, filesep ,fileList(i).name];
    fmrkr_in  = [mrkrList(i).folder, filesep ,mrkrList(i).name];
    disp(['IMPRT: ',fileList(i).name,' and started extraction.'])
    try
        [hdr,rec] = edfread(fname_in);
    catch
        disp('ERROR: EDF not readable.')
    end   
    mrk = find(rec(26,:) == ['1','3','25']);
    val = rec(26,mrk);
    
    switch length(mrk)
        case 0
            disp(['[ ',fileList(i).name,' ] ERROR: Includes',...
                num2str(length(mrk)),' markers.'])
            continue
        case 1
            load(fmrkr_in)
            [x,y] = find(SubjectTimeLog == "Instructions");
            mrk(2:5) = str2double(SubjectTimeLog(1,x));
        case 4
            disp(['[ ',fileList(i).name,' ] ERROR: Includes',...
                num2str(length(mrk)),' markers.'])
            continue
        case 5
            mrk      = mrk([1,3,5]);
            val      = val([1,3,5]);
            load(fmrkr_in)
            [x,y] = find(SubjectTimeLog == "Instructions");
            mrk(4:7) = str2double(SubjectTimeLog(1,x))*samp_freq;
            val(4:7) = 25;
        otherwise
            disp(['[ ',fileList(i).name,' ] ERROR: Includes',...
                num2str(length(mrk)),' markers.'])
            continue
    end
    
    % Cycle through markers _______________________________________________
    for p=1:length(mrk)
        switch val
            case 1
                mrkr_data(p) = [mrk(p) mrk(p)+15*samp_freq];
            case 3
                mrkr_data(p) = [mrk(p) mrk(p)+15*samp_freq];
            case 25
                mrkr_data(p) = [mrk(p) mrk(p)+72*samp_freq];
            otherwise
                disp(['[ ',fileList(i).name,' ] ERROR: Includes ',...
                    'unexpected marker value of ',num2str(length(mrk))])
        end
    end % _________________________________________________________________

    % Save data, marker ___________________________________________________
    fname_in = fullfile(fileList(i).folder, fileList(i).name);
    [~, lfn, ~] = fileparts(fileList(i).folder);
    fname_tm = ['subjectid_',lfn,'_mrk'];
    fname_ex = fullfile(...
                    append(fedf_save,...
                           strrep(fileList(i).folder,file_load,''),...
                           filesep,fname_tm));
    if ~isfile([fname_ex,'.edf'])
        movefile(fname_in, [fname_ex,'.edf']);
        save([fname_ex,'.mat'],"mrkr_data");
    else
        movefile(fname_in, [fname_ex,num2str(i),'.edf'])
        save([fname_ex,num2str(i),'.mat'],"mrkr_data");
    end
    disp(['SAVED: ',fname_ex,'.edf with markers:'])
    disp(mrkr_data/samp_freq)
    % _____________________________________________________________________
end