

settings.eegLabels = {'AF3', 'F7', 'F3', 'FC5', 'T7', 'P7', 'O1', ...
                    'O2', 'P8', 'T8', 'FC6', 'F4', 'F8', 'AF4'};
settings.bandPass = [0.5 45];
settings.hampelWindow = 11;
settings.overlap = 4.5;

func_EEG_spectral_analysis(fn, fs, nchan, ncomp, preprocess, varargin)

root_dirx = 'G:\Projects\OUHSCgl\prcEEG\';
edfr_locs = 'G:\Projects\OUHSCgl\_exe\';
eegl_locs = 'G:\Projects\OUHSCgl\_exe\eeglab2023.1\';

file_load = 'fpn\mat_data';
fcsv_save = 'aly\spectralanalysis_csv_data';

% Addition of files to workspace, initial user error handling and file
[debg_logs, save_logs, load_logs] = deal("","","");
% Add paths of dependecies to workspace ___________________________________
if isempty(eegl_locs)
    debg_logs(end+1,1) = mdc_errorlog(1, 'eegl_locs');
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


fileList = dir(fullfile(file_load,'**','.mat'));

fields = fieldnames(kwargs);
for n=1:numel(fields)
	kwargs.(fields{n}) = settings.(fields{n});
end

for i=1:length(fileList)
    fname_in  = [fileList(i).folder, filesep, fileList(i).name];
    data = load(fname_in);

end