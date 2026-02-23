%% Enhanced EEG Analysis Extensions
% This script provides additional analysis modules to extend the main EEG processing pipeline.
% It includes:
% 1. Statistical comparison between EEG and fNIRS results
% 2. Advanced visualization for multimodal connectivity
% 3. Source localization for EEG data
% 4. Machine learning classification approaches
% 5. Quality control metrics and reporting

%% 1. Statistical Correlation Analysis Between EEG and fNIRS Metrics
function runEEGfNIRSCorrelation(eeg_extn_path, nir_extn_path, tasks)
    % This function analyzes correlations between EEG and fNIRS metrics
    % across subjects to identify relationships between the modalities
    %
    % Usage:
    %   runEEGfNIRSCorrelation(eeg_extn_path, nir_extn_path, {'nbk', 'ftp'})
    
    disp('=== Starting EEG-fNIRS Statistical Correlation Analysis ===');
    
    if ~exist('tasks', 'var')
        tasks = {'nbk', 'ftp'};
    end
    
    % Create output directory
    output_dir = fullfile(eeg_extn_path, 'correlation_analysis');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Analyze each task
    for t = 1:length(tasks)
        task = tasks{t};
        disp(['Analyzing task: ' task]);
        
        % Load EEG metrics (power in key frequency bands)
        eeg_power_files = dir(fullfile(eeg_extn_path, task, '*_power_*.csv'));
        
        % Load fNIRS metrics (hemodynamic responses)
        fnirs_files = dir(fullfile(nir_extn_path, task, '*.csv'));
        fnirs_files = fnirs_files(~contains({fnirs_files.name}, 'group'));
        
        % Extract subject IDs
        eeg_subjects = cellfun(@extractSubjectID, {eeg_power_files.name}, 'UniformOutput', false);
        fnirs_subjects = cellfun(@extractSubjectID, {fnirs_files.name}, 'UniformOutput', false);
        
        % Find common subjects
        common_subjects = intersect(eeg_subjects, fnirs_subjects);
        
        if isempty(common_subjects)
            warning('No common subjects found between EEG and fNIRS data for task %s', task);
            continue;
        end
        
        disp(['Found ' num2str(length(common_subjects)) ' common subjects for analysis']);
        
        % Initialize correlation matrices
        eeg_bands = {'alpha', 'theta', 'beta', 'delta', 'gamma'};
        fnirs_types = {'hbo', 'hbr'};
        
        correlation_results = struct();
        
        % For each EEG band, correlate with fNIRS measures
        for band_idx = 1:length(eeg_bands)
            band = eeg_bands{band_idx};
            
            for fnirs_idx = 1:length(fnirs_types)
                fnirs_type = fnirs_types{fnirs_idx};
                
                % Initialize arrays for correlation
                eeg_values = [];
                fnirs_values = [];
                
                % Collect values for each subject
                for s = 1:length(common_subjects)
                    subject = common_subjects{s};
                    
                    % Find EEG file for this subject and band
                    eeg_idx = find(contains({eeg_power_files.name}, subject) & ...
                                  contains({eeg_power_files.name}, band));
                    
                    % Find fNIRS file for this subject
                    fnirs_idx = find(contains({fnirs_files.name}, subject));
                    
                    if ~isempty(eeg_idx) && ~isempty(fnirs_idx)
                        % Load EEG data
                        eeg_data = readtable(fullfile(eeg_power_files(eeg_idx(1)).folder, ...
                                                     eeg_power_files(eeg_idx(1)).name));
                        
                        % Load fNIRS data
                        fnirs_data = readtable(fullfile(fnirs_files(fnirs_idx(1)).folder, ...
                                                       fnirs_files(fnirs_idx(1)).name));
                        
                        % Extract mean EEG power for this band
                        mean_eeg_power = mean(eeg_data.Var2);
                        
                        % Extract mean fNIRS response for requested type
                        fnirs_cols = fnirs_data.Properties.VariableNames;
                        fnirs_type_cols = fnirs_cols(contains(fnirs_cols, fnirs_type));
                        
                        if ~isempty(fnirs_type_cols)
                            mean_fnirs_value = mean(fnirs_data{:, fnirs_type_cols}, 'omitnan');
                            
                            % Add to arrays
                            eeg_values = [eeg_values; mean_eeg_power];
                            fnirs_values = [fnirs_values; mean_fnirs_value];
                        end
                    end
                end
                
                % Calculate correlation if we have enough data points
                if length(eeg_values) >= 3
                    [r, p] = corrcoef(eeg_values, fnirs_values);
                    correlation_results.([band '_' fnirs_type]).r = r(1,2);
                    correlation_results.([band '_' fnirs_type]).p = p(1,2);
                    correlation_results.([band '_' fnirs_type]).n = length(eeg_values);
                    
                    % Create scatter plot
                    figure;
                    scatter(eeg_values, fnirs_values, 50, 'filled');
                    hold on;
                    
                    % Add fit line
                    p = polyfit(eeg_values, fnirs_values, 1);
                    x_fit = linspace(min(eeg_values), max(eeg_values), 100);
                    y_fit = polyval(p, x_fit);
                    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
                    
                    % Add labels
                    xlabel(['EEG ' band ' band power']);
                    ylabel(['fNIRS ' fnirs_type ' response']);
                    title(['Correlation between EEG ' band ' and fNIRS ' fnirs_type ...
                           ' (r = ' num2str(r(1,2), '%.2f') ...
                           ', p = ' num2str(p(1,2), '%.3f') ')']);
                    grid on;
                    
                    % Save figure
                    saveas(gcf, fullfile(output_dir, [task '_corr_' band '_' fnirs_type '.png']));
                    saveas(gcf, fullfile(output_dir, [task '_corr_' band '_' fnirs_type '.fig']));
                    close;
                end
            end
        end
        
        % Save correlation results
        save(fullfile(output_dir, [task '_correlation_results.mat']), 'correlation_results');
        
        % Create summary table
        summary_table = table();
        bands_all = {};
        fnirs_all = {};
        r_values = [];
        p_values = [];
        n_values = [];
        
        fields = fieldnames(correlation_results);
        for f = 1:length(fields)
            parts = strsplit(fields{f}, '_');
            bands_all{end+1} = parts{1};
            fnirs_all{end+1} = parts{2};
            r_values(end+1) = correlation_results.(fields{f}).r;
            p_values(end+1) = correlation_results.(fields{f}).p;
            n_values(end+1) = correlation_results.(fields{f}).n;
        end
        
        summary_table.EEG_Band = bands_all';
        summary_table.NIRS_Type = fnirs_all';
        summary_table.Correlation = r_values';
        summary_table.P_Value = p_values';
        summary_table.N_Subjects = n_values';
        
        % Save summary table
        writetable(summary_table, fullfile(output_dir, [task '_correlation_summary.csv']));
        
        disp(['Completed correlation analysis for task: ' task]);
    end
    
    disp('=== EEG-fNIRS Correlation Analysis Completed ===');
end

%% 2. Advanced Multimodal Connectivity Visualization
function createAdvancedConnectivityVis(eeg_extn_path, nir_extn_path, tasks)
    % This function creates advanced visualizations for multimodal connectivity
    % including:
    % - Network graphs
    % - Circular plots
    % - Brain maps
    %
    % Usage:
    %   createAdvancedConnectivityVis(eeg_extn_path, nir_extn_path, {'nbk', 'ftp'})
    
    disp('=== Starting Advanced Connectivity Visualization ===');
    
    if ~exist('tasks', 'var')
        tasks = {'nbk', 'ftp'};
    end
    
    % Process each task
    for t = 1:length(tasks)
        task = tasks{t};
        disp(['Creating visualizations for task: ' task]);
        
        % Output directory
        output_dir = fullfile(eeg_extn_path, 'multimodal', task, 'advanced_vis');
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        
        % Find connectivity matrices
        conn_files = dir(fullfile(eeg_extn_path, 'multimodal', task, '*_multimodal_conn.csv'));
        
        if isempty(conn_files)
            warning('No connectivity files found for task %s', task);
            continue;
        end
        
        disp(['Found ' num2str(length(conn_files)) ' connectivity files']);
        
        % Process each subject
        for f = 1:length(conn_files)
            % Extract subject ID
            [~, file_name, ~] = fileparts(conn_files(f).name);
            parts = strsplit(file_name, '_');
            subject_id = parts{1};
            
            disp(['Processing subject: ' subject_id]);
            
            % Load connectivity matrix
            conn_matrix = csvread(fullfile(conn_files(f).folder, conn_files(f).name));
            
            % 1. Create network graph
            try
                createNetworkGraph(conn_matrix, subject_id, output_dir);
            catch e
                warning('Error creating network graph: %s', e.message);
            end
            
            % 2. Create circular plot
            try
                createCircularPlot(conn_matrix, subject_id, output_dir);
            catch e
                warning('Error creating circular plot: %s', e.message);
            end
            
            % 3. Create brain map
            try
                createBrainMap(conn_matrix, subject_id, output_dir);
            catch e
                warning('Error creating brain map: %s', e.message);
            end
        end
        
        % Create group average visualization
        try
            % Load all matrices
            all_matrices = zeros(size(conn_matrix, 1), size(conn_matrix, 2), length(conn_files));
            
            for f = 1:length(conn_files)
                all_matrices(:,:,f) = csvread(fullfile(conn_files(f).folder, conn_files(f).name));
            end
            
            % Calculate average
            avg_matrix = mean(all_matrices, 3);
            
            % Create visualizations for group average
            createNetworkGraph(avg_matrix, 'group_average', output_dir);
            createCircularPlot(avg_matrix, 'group_average', output_dir);
            createBrainMap(avg_matrix, 'group_average', output_dir);
            
            % Save average matrix
            csvwrite(fullfile(output_dir, 'group_average_conn.csv'), avg_matrix);
        catch e
            warning('Error creating group average visualizations: %s', e.message);
        end
    end
    
    disp('=== Advanced Connectivity Visualization Completed ===');
end

function createNetworkGraph(conn_matrix, subject_id, output_dir)
    % Create network graph visualization
    
    % Threshold the matrix to show only strong connections
    threshold = mean(conn_matrix(:)) + 0.5 * std(conn_matrix(:));
    thresholded_matrix = conn_matrix;
    thresholded_matrix(thresholded_matrix < threshold) = 0;
    
    % Create network graph
    figure('Position', [100, 100, 800, 600]);
    
    % Calculate number of nodes
    n_eeg = size(conn_matrix, 1);
    n_nirs = size(conn_matrix, 2);
    
    % Create node positions
    theta = linspace(0, 2*pi, n_eeg + n_nirs + 1);
    theta = theta(1:end-1);
    
    % EEG nodes on left, NIRS nodes on right
    x = cos(theta);
    y = sin(theta);
    
    % Plot nodes
    hold on;
    % EEG nodes (blue)
    scatter(x(1:n_eeg), y(1:n_eeg), 300, 'b', 'filled');
    % NIRS nodes (red)
    scatter(x(n_eeg+1:end), y(n_eeg+1:end), 300, 'r', 'filled');
    
    % Add node labels
    for i = 1:n_eeg
        text(x(i)*1.1, y(i)*1.1, ['EEG' num2str(i)], 'FontSize', 12, 'HorizontalAlignment', 'center');
    end
    
    for i = 1:n_nirs
        text(x(n_eeg+i)*1.1, y(n_eeg+i)*1.1, ['NIRS' num2str(i)], 'FontSize', 12, 'HorizontalAlignment', 'center');
    end
    
    % Plot connections
    max_conn = max(thresholded_matrix(:));
    for i = 1:n_eeg
        for j = 1:n_nirs
            if thresholded_matrix(i,j) > 0
                % Scale line width by connection strength
                line_width = 0.5 + 4 * thresholded_matrix(i,j) / max_conn;
                
                % Scale color by connection strength (0=blue, 1=red)
                color_val = thresholded_matrix(i,j) / max_conn;
                conn_color = [color_val, 0, 1-color_val];
                
                % Draw line
                line([x(i), x(n_eeg+j)], [y(i), y(n_eeg+j)], ...
                     'LineWidth', line_width, 'Color', conn_color);
            end
        end
    end
    
    % Add title and remove axes
    title(['EEG-fNIRS Connectivity Network: ' strrep(subject_id, '_', '\_')], 'FontSize', 16);
    axis equal off;
    
    % Save the figure
    saveas(gcf, fullfile(output_dir, [subject_id '_network.png']));
    saveas(gcf, fullfile(output_dir, [subject_id '_network.fig']));
    close;
end

function createCircularPlot(conn_matrix, subject_id, output_dir)
    % Create circular plot visualization
    
    % Calculate number of nodes
    n_eeg = size(conn_matrix, 1);
    n_nirs = size(conn_matrix, 2);
    
    % Create figure
    figure('Position', [100, 100, 800, 800]);
    
    % Calculate positions on the circle
    theta = linspace(0, 2*pi, n_eeg + n_nirs + 1);
    theta = theta(1:end-1);
    radius = 5;
    
    % Node positions
    x = radius * cos(theta);
    y = radius * sin(theta);
    
    % Setup plot area
    axis equal;
    axis([-radius*1.2, radius*1.2, -radius*1.2, radius*1.2]);
    hold on;
    
    % Draw circular outline
    plot(radius * cos(linspace(0, 2*pi, 100)), radius * sin(linspace(0, 2*pi, 100)), 'k-', 'LineWidth', 1);
    
    % Draw EEG node points (first half)
    for i = 1:n_eeg
        plot(x(i), y(i), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
        text(x(i)*1.1, y(i)*1.1, ['EEG' num2str(i)], 'FontSize', 10);
    end
    
    % Draw NIRS node points (second half)
    for i = 1:n_nirs
        plot(x(n_eeg+i), y(n_eeg+i), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        text(x(n_eeg+i)*1.1, y(n_eeg+i)*1.1, ['NIRS' num2str(i)], 'FontSize', 10);
    end
    
    % Draw connections
    % Normalize connection strengths
    max_conn = max(conn_matrix(:));
    min_conn = min(conn_matrix(:));
    norm_matrix = (conn_matrix - min_conn) / (max_conn - min_conn);
    
    % Threshold to only show moderately strong connections
    threshold = 0.3;
    
    % Draw connection curves
    for i = 1:n_eeg
        for j = 1:n_nirs
            if norm_matrix(i,j) > threshold
                % Calculate control points for Bezier curve
                % The control point is inside the circle to create a curve
                ctrl_factor = 0.5; % How much it curves inward
                ctrl_point = ctrl_factor * [x(i) + x(n_eeg+j), y(i) + y(n_eeg+j)];
                
                % Create Bezier curve points
                t = linspace(0, 1, 50);
                curve_x = (1-t).^2 .* x(i) + 2*(1-t).*t .* ctrl_point(1) + t.^2 .* x(n_eeg+j);
                curve_y = (1-t).^2 .* y(i) + 2*(1-t).*t .* ctrl_point(2) + t.^2 .* y(n_eeg+j);
                
                % Line width based on connection strength
                line_width = 0.5 + 3 * norm_matrix(i,j);
                
                % Color based on connection strength (blue to red)
                conn_color = [norm_matrix(i,j), 0, 1-norm_matrix(i,j)];
                
                % Plot the curve
                plot(curve_x, curve_y, 'Color', conn_color, 'LineWidth', line_width);
            end
        end
    end
    
    % Add title
    title(['Circular Connectivity Plot: ' strrep(subject_id, '_', '\_')], 'FontSize', 16);
    axis off;
    
    % Save the figure
    saveas(gcf, fullfile(output_dir, [subject_id '_circular.png']));
    saveas(gcf, fullfile(output_dir, [subject_id '_circular.fig']));
    close;
end

function createBrainMap(conn_matrix, subject_id, output_dir)
    % Create a simplified brain map visualization
    
    figure('Position', [100, 100, 1000, 500]);
    
    % Left side: Top view of brain
    subplot(1, 2, 1);
    
    % Draw simple brain outline
    t = linspace(0, 2*pi, 100);
    x = cos(t);
    y = sin(t);
    
    % Brain outline (slightly egg-shaped)
    x_brain = 10 * (x * 0.8 + 0.2 * (x .* (x > 0)));
    y_brain = 10 * y;
    
    % Plot brain outline
    plot(x_brain, y_brain, 'k-', 'LineWidth', 2);
    hold on;
    
    % Plot midline
    plot([0, 0], [-10, 10], 'k--', 'LineWidth', 1);
    
    % Calculate number of nodes
    n_eeg = size(conn_matrix, 1);
    n_nirs = size(conn_matrix, 2);
    
    % Place EEG channels at appropriate positions (simplified 10-20 system)
    eeg_pos_x = [];
    eeg_pos_y = [];
    
    % Simple distribution of EEG channels
    angle_step = pi / (n_eeg + 1);
    for i = 1:n_eeg
        angle = pi/2 - i * angle_step;
        eeg_pos_x(i) = 7 * cos(angle);
        eeg_pos_y(i) = 7 * sin(angle);
        
        % Plot EEG channel
        plot(eeg_pos_x(i), eeg_pos_y(i), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
        text(eeg_pos_x(i)+0.5, eeg_pos_y(i)+0.5, ['E' num2str(i)], 'FontSize', 10);
    end
    
    % Place NIRS optodes
    nirs_pos_x = [];
    nirs_pos_y = [];
    
    % Simple distribution of NIRS channels
    for i = 1:n_nirs
        angle = pi/2 - i * angle_step;
        nirs_pos_x(i) = 5 * cos(angle);
        nirs_pos_y(i) = 5 * sin(angle);
        
        % Plot NIRS channel
        plot(nirs_pos_x(i), nirs_pos_y(i), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        text(nirs_pos_x(i)+0.5, nirs_pos_y(i)+0.5, ['N' num2str(i)], 'FontSize', 10);
    end
    
    % Draw connections
    max_conn = max(conn_matrix(:));
    for i = 1:n_eeg
        for j = 1:n_nirs
            if conn_matrix(i,j) > 0.5 * max_conn
                % Calculate control points for Bezier curve
                ctrl_factor = 0.2;
                ctrl_x = (eeg_pos_x(i) + nirs_pos_x(j))/2;
                ctrl_y = (eeg_pos_y(i) + nirs_pos_y(j))/2 - ctrl_factor * 10;
                
                % Create Bezier curve points
                t = linspace(0, 1, 50);
                curve_x = (1-t).^2 .* eeg_pos_x(i) + 2*(1-t).*t .* ctrl_x + t.^2 .* nirs_pos_x(j);
                curve_y = (1-t).^2 .* eeg_pos_y(i) + 2*(1-t).*t .* ctrl_y + t.^2 .* nirs_pos_y(j);
                
                % Line width based on connection strength
                line_width = 0.5 + 2 * conn_matrix(i,j) / max_conn;
                
                % Color based on connection strength (blue to red)
                conn_color = [conn_matrix(i,j)/max_conn, 0, 1-conn_matrix(i,j)/max_conn];
                
                % Plot the curve
                plot(curve_x, curve_y, 'Color', conn_color, 'LineWidth', line_width);
            end
        end
    end
    
    % Title and axis settings
    title(['Brain Map (Top View): ' strrep(subject_id, '_', '\_')], 'FontSize', 14);
    axis equal off;
    
    % Right side: Connection matrix as heatmap
    subplot(1, 2, 2);
    imagesc(conn_matrix);
    colormap(jet);
    colorbar;
    
    % Add labels
    xlabel('NIRS Channels');
    ylabel('EEG Channels');
    title('Connectivity Matrix', 'FontSize', 14);
    
    % Add channel numbers as ticks
    xticks(1:n_nirs);
    yticks(1:n_eeg);
    xticklabels(arrayfun(@num2str, 1:n_nirs, 'UniformOutput', false));
    yticklabels(arrayfun(@num2str, 1:n_eeg, 'UniformOutput', false));
    
    % Main title
    suptitle(['EEG-fNIRS Connectivity: ' strrep(subject_id, '_', '\_')]);
    
    % Save the figure
    saveas(gcf, fullfile(output_dir, [subject_id '_brainmap.png']));
    saveas(gcf, fullfile(output_dir, [subject_id '_brainmap.fig']));
    close;
end

%% 3. Source Localization for EEG Data
function runEEGSourceLocalization(eeg_extn_path, task, event_type)
    % This function performs basic source localization for EEG data
    % Uses EEGLAB's DIPFIT plugin if available
    %
    % Usage:
    %   runEEGSourceLocalization(eeg_extn_path, 'nbk', 'n2a')
    
    disp('=== Starting EEG Source Localization ===');
    
    % Check for DIPFIT plugin
    try
        which('pop_dipfit_settings');
    catch
        warning('DIPFIT plugin not found. Please install it in EEGLAB.');
        return;
    end
    
    % Create output directory
    output_dir = fullfile(eeg_extn_path, task, 'source_localization');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Find preprocessed EEG files
    eeg_files = dir(fullfile(eeg_extn_path, task, 'preprocessed', '*.set'));
    
    if isempty(eeg_files)
        warning('No preprocessed EEG files found for task %s', task);
        return;
    end
    
    disp(['Found ' num2str(length(eeg_files)) ' EEG files for source localization']);
    
    % Process each file
    for f = 1:length(eeg_files)
        file_path = fullfile(eeg_files(f).folder, eeg_files(f).name);
        [~, file_base, ~] = fileparts(eeg_files(f).name);
        
        disp(['Processing ' file_base '...']);
        
        try
            % Load dataset
            EEG = pop_loadset(file_path);
            
            % Check for channel locations
            if ~isfield(EEG, 'chanlocs') || isempty(EEG.chanlocs) || ...
               ~isfield(EEG.chanlocs, 'X') || isempty([EEG.chanlocs.X])
                warning('Channel locations not available for %s. Skipping...', file_base);
                continue;
            end
            
            % Find epochs of the requested event type
            if ~isempty(event_type) && EEG.trials > 1
                event_epochs = find(strcmp({EEG.epoch.eventtype}, event_type));
                
                if isempty(event_epochs)
                    warning('Event type %s not found in %s. Skipping...', event_type, file_base);
                    continue;
                end
                
                % Create average ERP for this event type
                erp_data = mean(EEG.data(:,:,event_epochs), 3);
            else
                % Use all data
                if EEG.trials > 1
                    erp_data = mean(EEG.data, 3);
                else
                    % For continuous data, just use a segment
                    segment_length = min(1000, EEG.pnts);
                    start_point = round(EEG.pnts/2 - segment_length/2);
                    erp_data = EEG.data(:, start_point:(start_point+segment_length-1));
                end
            end
            
            % Create a new dataset with just the ERP
            ERP = EEG;
            ERP.data = erp_data;
            ERP.trials = 1;
            ERP.pnts = size(erp_data, 2);
            
            % Setup the head model
            ERP = pop_dipfit_settings(ERP, ...
                'hdmfile', 'standard_BEM/standard_vol.mat', ...
                'coordformat', 'MNI', ...
                'mrifile', 'standard_BEM/standard_mri.mat', ...
                'chanfile', 'standard_BEM/elec/standard_1005.elc', ...
                'chansel', 1:ERP.nbchan);
            
            % Find component maps
            ERP = pop_multifit(ERP, 1:min(10, ERP.nbchan), 'threshold', 100, 'dipplot', 'off');
            
            % Plot dipoles
            pop_dipplot(ERP, 1:min(5, ERP.nbchan), 'mri', 'standard_BEM/standard_mri.mat', ...
                'normlen', 'on', 'num', 'on');
            
            % Save figure
            saveas(gcf, fullfile(output_dir, [file_base '_dipoles.png']));
            saveas(gcf, fullfile(output_dir, [file_base '_dipoles.fig']));
            close;
            
            % Save localization results
            save(fullfile(output_dir, [file_base '_dipfit.mat']), 'ERP');
            
            disp(['Completed source localization for ' file_base]);
        catch e
            warning('Error during source localization for %s: %s', file_base, e.message);
        end
    end
    
    disp('=== Source Localization Completed ===');
end

%% 4. Machine Learning Classification
function runMLClassification(eeg_extn_path, task)
    % This function performs machine learning classification on EEG data
    % to distinguish between different event types
    %
    % Usage:
    %   runMLClassification(eeg_extn_path, 'nbk')
    
    disp('=== Starting Machine Learning Classification ===');
    
    % Create output directory
    output_dir = fullfile(eeg_extn_path, task, 'ml_classification');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Find preprocessed EEG files
    eeg_files = dir(fullfile(eeg_extn_path, task, 'preprocessed', '*.set'));
    
    if isempty(eeg_files)
        warning('No preprocessed EEG files found for task %s', task);
        return;
    end
    
    disp(['Found ' num2str(length(eeg_files)) ' EEG files for classification']);
    
    % Determine event types based on task
    if strcmp(task, 'nbk')
        event_types = {'rEC', 'rEO', 'n0a', 'n1a', 'n0b', 'n2a'};
        % For binary classification, choose two event types
        binary_classes = {'n0a', 'n2a'}; % e.g., 0-back vs 2-back
    elseif strcmp(task, 'ftp')
        event_types = {'L01', 'R01', 'L02', 'R02', 'L03', 'R03'};
        % For binary classification, choose left vs right
        binary_classes = {'L01', 'R01'}; % e.g., left vs right finger tapping
    else
        warning('Unknown task type %s', task);
        return;
    end
    
    % Initialize arrays for features and labels
    all_features = [];
    all_labels = [];
    subject_ids = {};
    
    % Extract features for each file
    for f = 1:length(eeg_files)
        file_path = fullfile(eeg_files(f).folder, eeg_files(f).name);
        [~, file_base, ~] = fileparts(eeg_files(f).name);
        
        disp(['Extracting features from ' file_base '...']);
        
        try
            % Load dataset
            EEG = pop_loadset(file_path);
            
            % Check if the data is epoched
            if EEG.trials <= 1
                warning('Data in %s is not epoched. Skipping...', file_base);
                continue;
            end
            
            % Extract subject ID
            subject_id = extractSubjectID(file_base);
            
            % Find epochs of the requested binary classes
            class1_epochs = find(strcmp({EEG.epoch.eventtype}, binary_classes{1}));
            class2_epochs = find(strcmp({EEG.epoch.eventtype}, binary_classes{2}));
            
            if isempty(class1_epochs) || isempty(class2_epochs)
                warning('One or both classes not found in %s. Skipping...', file_base);
                continue;
            end
            
            % Extract features for each epoch
            for epoch_idx = [class1_epochs; class2_epochs]'
                % Get epoch data
                epoch_data = EEG.data(:,:,epoch_idx);
                
                % Extract features
                % 1. Spectral features
                % Calculate power in different frequency bands
                features = [];
                
                freq_bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
                band_ranges = {[0.5 4], [4 8], [8 13], [13 30], [30 45]};
                
                for b = 1:length(freq_bands)
                    band = band_ranges{b};
                    % Filter data to the band
                    filtered_data = eegfilt(epoch_data, EEG.srate, band(1), band(2), 0, []);
                    
                    % Calculate band power for each channel
                    band_power = mean(filtered_data.^2, 2);
                    
                    % Add to features
                    features = [features; band_power];
                end
                
                % 2. Time domain features
                % Calculate mean, std, min, max for each channel
                channel_mean = mean(epoch_data, 2);
                channel_std = std(epoch_data, 0, 2);
                channel_min = min(epoch_data, [], 2);
                channel_max = max(epoch_data, [], 2);
                
                features = [features; channel_mean; channel_std; channel_min; channel_max];
                
                % 3. Connectivity features (simplification)
                % Use a subset of channels for simplicity
                channels_subset = min(5, EEG.nbchan);
                if channels_subset > 1
                    conn_matrix = zeros(channels_subset, channels_subset);
                    for ch1 = 1:channels_subset
                        for ch2 = (ch1+1):channels_subset
                            conn_matrix(ch1, ch2) = corr(epoch_data(ch1,:)', epoch_data(ch2,:)');
                            conn_matrix(ch2, ch1) = conn_matrix(ch1, ch2);
                        end
                    end
                    
                    % Add upper triangle of connectivity matrix to features
                    upper_triu = conn_matrix(triu(true(size(conn_matrix)), 1));
                    features = [features; upper_triu];
                end
                
                % Add to the overall feature matrix
                all_features = [all_features, features];
                
                % Add label (1 for class1, 2 for class2)
                if ismember(epoch_idx, class1_epochs)
                    all_labels = [all_labels, 1];
                else
                    all_labels = [all_labels, 2];
                end
                
                % Add subject ID
                subject_ids{end+1} = subject_id;
            end
        catch e
            warning('Error extracting features from %s: %s', file_base, e.message);
        end
    end
    
    % Check if we have enough data
    if isempty(all_features) || size(all_features, 2) < 10
        warning('Not enough data for classification');
        return;
    end
    
    disp(['Extracted features from ' num2str(size(all_features, 2)) ' epochs']);
    
    % Normalize features
    all_features = zscore(all_features, [], 2);
    
    % Perform classification
    try
        % Convert to expected format
        X = all_features';
        y = all_labels';
        subject_ids = subject_ids';
        
        % Get unique subjects
        unique_subjects = unique(subject_ids);
        
        % Results storage
        all_accuracies = [];
        all_cm = zeros(2, 2);
        
        % Cross-validation by subject (leave-one-subject-out)
        for s = 1:length(unique_subjects)
            test_subject = unique_subjects{s};
            
            % Split data
            test_idx = strcmp(subject_ids, test_subject);
            train_idx = ~test_idx;
            
            X_train = X(train_idx, :);
            y_train = y(train_idx);
            X_test = X(test_idx, :);
            y_test = y(test_idx);
            
            if isempty(X_train) || isempty(X_test)
                continue;
            end
            
            % Train SVM classifier
            svm_model = fitcsvm(X_train, y_train, 'KernelFunction', 'rbf', ...
                'Standardize', true, 'ClassNames', [1, 2]);
            
            % Make predictions
            y_pred = predict(svm_model, X_test);
            
            % Calculate accuracy
            accuracy = sum(y_pred == y_test) / length(y_test);
            all_accuracies = [all_accuracies; accuracy];
            
            % Calculate confusion matrix
            cm = confusionmat(y_test, y_pred);
            
            % Ensure cm is 2x2
            if size(cm, 1) == 1
                if cm(1,1) == 1
                    cm = [cm 0; 0 0];
                else
                    cm = [0 0; 0 cm(1,1)];
                end
            end
            
            all_cm = all_cm + cm;
            
            disp(['Subject ' test_subject ': Classification accuracy = ' num2str(accuracy*100, '%.1f') '%']);
        end
        
        % Calculate overall performance
        mean_accuracy = mean(all_accuracies);
        std_accuracy = std(all_accuracies);
        
        % Calculate sensitivity, specificity, F1 score
        TP = all_cm(1,1);
        FP = all_cm(2,1);
        FN = all_cm(1,2);
        TN = all_cm(2,2);
        
        sensitivity = TP / (TP + FN);
        specificity = TN / (TN + FP);
        precision = TP / (TP + FP);
        recall = sensitivity;
        f1_score = 2 * precision * recall / (precision + recall);
        
        % Save results as text file
        fid = fopen(fullfile(output_dir, [task '_classification_results.txt']), 'w');
        fprintf(fid, 'EEG Classification Results\n');
        fprintf(fid, '-------------------------\n');
        fprintf(fid, 'Task: %s\n', task);
        fprintf(fid, 'Classes: %s vs %s\n', binary_classes{1}, binary_classes{2});
        fprintf(fid, 'Number of subjects: %d\n', length(unique_subjects));
        fprintf(fid, 'Number of epochs: %d\n', size(X, 1));
        fprintf(fid, 'Number of features: %d\n', size(X, 2));
        fprintf(fid, '\nPerformance Metrics:\n');
        fprintf(fid, '  Mean Accuracy: %.2f%% (±%.2f%%)\n', mean_accuracy*100, std_accuracy*100);
        fprintf(fid, '  Sensitivity: %.2f%%\n', sensitivity*100);
        fprintf(fid, '  Specificity: %.2f%%\n', specificity*100);
        fprintf(fid, '  Precision: %.2f%%\n', precision*100);
        fprintf(fid, '  Recall: %.2f%%\n', recall*100);
        fprintf(fid, '  F1 Score: %.2f\n', f1_score);
        fprintf(fid, '\nConfusion Matrix:\n');
        fprintf(fid, '  %d\t%d\n', all_cm(1,1), all_cm(1,2));
        fprintf(fid, '  %d\t%d\n', all_cm(2,1), all_cm(2,2));
        fprintf(fid, '\nPer-Subject Accuracies:\n');
        for s = 1:length(unique_subjects)
            if s <= length(all_accuracies)
                fprintf(fid, '  %s: %.2f%%\n', unique_subjects{s}, all_accuracies(s)*100);
            end
        end
        fclose(fid);
        
        % Create visualization
        figure('Position', [100, 100, 800, 600]);
        
        % Confusion matrix
        subplot(2, 2, 1);
        imagesc(all_cm);
        colormap('hot');
        colorbar;
        title('Confusion Matrix');
        xlabel('Predicted Class');
        ylabel('True Class');
        xticks([1 2]);
        yticks([1 2]);
        xticklabels(binary_classes);
        yticklabels(binary_classes);
        
        % Add text labels to confusion matrix
        for i = 1:2
            for j = 1:2
                text(j, i, num2str(all_cm(i,j)), 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'FontWeight', 'bold');
            end
        end
        
        % Bar chart of performance metrics
        subplot(2, 2, 2);
        metrics = [mean_accuracy, sensitivity, specificity, precision, f1_score];
        bar(metrics * 100);
        ylim([0 100]);
        title('Performance Metrics');
        xticklabels({'Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'F1'});
        ylabel('Percentage (%)');
        
        % Individual subject accuracies
        subplot(2, 2, [3 4]);
        bar(all_accuracies * 100);
        hold on;
        plot([0, length(all_accuracies)+1], [mean_accuracy*100, mean_accuracy*100], 'r--', 'LineWidth', 2);
        xlim([0, length(all_accuracies)+1]);
        ylim([0, 100]);
        title('Individual Subject Accuracies');
        xlabel('Subject');
        ylabel('Accuracy (%)');
        xticks(1:length(all_accuracies));
        xticklabels(unique_subjects);
        
        % Add main title
        suptitle(['EEG Classification: ' task ' - ' binary_classes{1} ' vs ' binary_classes{2}]);
        
        % Save figure
        saveas(gcf, fullfile(output_dir, [task '_classification_results.png']));
        saveas(gcf, fullfile(output_dir, [task '_classification_results.fig']));
        close;
        
        disp(['Classification results saved to ' fullfile(output_dir, [task '_classification_results.txt'])]);
    catch e
        warning('Error during classification: %s', e.message);
    end
    
    disp('=== Machine Learning Classification Completed ===');
end

%% 5. Quality Control and Reporting
function generateQualityReport(eeg_extn_path, nir_extn_path, tasks)
    % This function generates a comprehensive quality report for the entire dataset
    % Includes data quality metrics, processing statistics, and result summaries
    %
    % Usage:
    %   generateQualityReport(eeg_extn_path, nir_extn_path, {'nbk', 'ftp'})
    
    disp('=== Starting Quality Control and Reporting ===');
    
    if ~exist('tasks', 'var')
        tasks = {'nbk', 'ftp'};
    end
    
    % Create output directory
    output_dir = fullfile(eeg_extn_path, 'quality_reports');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Process each task
    for t = 1:length(tasks)
        task = tasks{t};
        disp(['Generating report for task: ' task]);
        
        % Create task output directory
        task_dir = fullfile(output_dir, task);
        if ~exist(task_dir, 'dir')
            mkdir(task_dir);
        end
        
        % Find all EEG files
        eeg_files = dir(fullfile(eeg_extn_path, task, 'preprocessed', '*.set'));
        
        if isempty(eeg_files)
            warning('No preprocessed EEG files found for task %s', task);
            continue;
        end
        
        % Find all fNIRS files
        fnirs_files = dir(fullfile(nir_extn_path, task, '*.csv'));
        fnirs_files = fnirs_files(~contains({fnirs_files.name}, 'group'));
        
        % Generate Report
        % Create main HTML file
        html_file = fullfile(task_dir, [task '_quality_report.html']);
        fid = fopen(html_file, 'w');
        
        % HTML Header
        fprintf(fid, '<!DOCTYPE html>\n');
        fprintf(fid, '<html>\n<head>\n');
        fprintf(fid, '<title>EEG-fNIRS Quality Report: %s</title>\n', task);
        fprintf(fid, '<style>\n');
        fprintf(fid, 'body { font-family: Arial, sans-serif; margin: 40px; }\n');
        fprintf(fid, 'h1, h2, h3 { color: #2c3e50; }\n');
        fprintf(fid, 'table { border-collapse: collapse; width: 100%%; margin-bottom: 20px; }\n');
        fprintf(fid, 'th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n');
        fprintf(fid, 'th { background-color: #f2f2f2; }\n');
        fprintf(fid, 'tr:nth-child(even) { background-color: #f9f9f9; }\n');
        fprintf(fid, '.good { color: green; }\n');
        fprintf(fid, '.warning { color: orange; }\n');
        fprintf(fid, '.error { color: red; }\n');
        fprintf(fid, '.summary { background-color: #eef; padding: 15px; border-radius: 5px; margin-bottom: 20px; }\n');
        fprintf(fid, '.section { margin-bottom: 30px; }\n');
        fprintf(fid, '.chart-container { display: flex; justify-content: space-around; margin-bottom: 30px; }\n');
        fprintf(fid, '.chart { width: 45%%; }\n');
        fprintf(fid, '</style>\n');
        fprintf(fid, '</head>\n<body>\n');
        
        % Report Header
        fprintf(fid, '<h1>EEG-fNIRS Quality Report: %s</h1>\n', upper(task));
        fprintf(fid, '<p>Generated on %s</p>\n', datestr(now));
        
        % Summary section
        fprintf(fid, '<div class="summary section">\n');
        fprintf(fid, '<h2>Summary</h2>\n');
        fprintf(fid, '<p>Total EEG datasets: %d</p>\n', length(eeg_files));
        fprintf(fid, '<p>Total fNIRS datasets: %d</p>\n', length(fnirs_files));
        
        % Extract subject IDs
        eeg_subjects = cellfun(@extractSubjectID, {eeg_files.name}, 'UniformOutput', false);
        fnirs_subjects = cellfun(@extractSubjectID, {fnirs_files.name}, 'UniformOutput', false);
        common_subjects = intersect(eeg_subjects, fnirs_subjects);
        
        fprintf(fid, '<p>Subjects with both EEG and fNIRS: %d</p>\n', length(common_subjects));
        fprintf(fid, '</div>\n');
        
        % Process EEG datasets
        fprintf(fid, '<div class="section">\n');
        fprintf(fid, '<h2>EEG Data Quality</h2>\n');
        fprintf(fid, '<table>\n');
        fprintf(fid, '<tr><th>Subject</th><th>Channels</th><th>Sampling Rate</th><th>Duration (s)</th><th>Epochs</th><th>Rejected (%%)</th><th>SNR (dB)</th><th>Status</th></tr>\n');
        
        eeg_quality = struct('subject', {}, 'channels', {}, 'srate', {}, 'duration', {}, ...
                          'epochs', {}, 'rejected_pct', {}, 'snr', {}, 'status', {});
        
        for f = 1:length(eeg_files)
            try
                % Load EEG data
                EEG = pop_loadset(fullfile(eeg_files(f).folder, eeg_files(f).name));
                
                % Calculate quality metrics
                subject = extractSubjectID(eeg_files(f).name);
                channels = EEG.nbchan;
                srate = EEG.srate;
                
                if EEG.trials > 1
                    duration = EEG.pnts * EEG.trials / EEG.srate;
                    epochs = EEG.trials;
                    
                    % Estimate reject rate (based on metadata if available)
                    if isfield(EEG, 'etc') && isfield(EEG.etc, 'clean_sample_mask')
                        rejected_pct = 100 * (1 - mean(EEG.etc.clean_sample_mask));
                    else
                        rejected_pct = NaN;
                    end
                else
                    duration = EEG.pnts / EEG.srate;
                    epochs = 0;
                    rejected_pct = NaN;
                end
                
                % Calculate SNR (signal-to-noise ratio)
                % Simple approach: ratio of signal power to baseline noise power
                if EEG.trials > 1
                    % For epoched data, use average ERP vs single-trial variance
                    mean_signal = mean(EEG.data, 3);
                    noise = EEG.data - repmat(mean_signal, [1, 1, EEG.trials]);
                    
                    signal_power = mean(mean_signal.^2, 'all');
                    noise_power = mean(noise.^2, 'all');
                    
                    if noise_power > 0
                        snr = 10 * log10(signal_power / noise_power);
                    else
                        snr = Inf;
                    end
                else
                    % For continuous data, estimate from spectrum
                    [spectrum, freqs] = spectopo(EEG.data, 0, EEG.srate, 'plot', 'off');
                    
                    % Find signal vs noise frequency bands
                    signal_idx = find(freqs >= 1 & freqs <= 30);
                    noise_idx = find(freqs > 30 & freqs <= 50);
                    
                    if ~isempty(signal_idx) && ~isempty(noise_idx)
                        signal_power = mean(10.^(spectrum(:,signal_idx)/10), 'all');
                        noise_power = mean(10.^(spectrum(:,noise_idx)/10), 'all');
                        
                        if noise_power > 0
                            snr = 10 * log10(signal_power / noise_power);
                        else
                            snr = Inf;
                        end
                    else
                        snr = NaN;
                    end
                end
                
                % Determine status
                if snr > 10
                    status = 'good';
                    status_label = '<span class="good">Good</span>';
                elseif snr > 5
                    status = 'acceptable';
                    status_label = '<span class="warning">Acceptable</span>';
                else
                    status = 'poor';
                    status_label = '<span class="error">Poor</span>';
                end
                
                % Write to table
                fprintf(fid, '<tr><td>%s</td><td>%d</td><td>%.1f</td><td>%.1f</td><td>%d</td><td>%.1f</td><td>%.1f</td><td>%s</td></tr>\n', ...
                      subject, channels, srate, duration, epochs, rejected_pct, snr, status_label);
                
                % Save to structure
                eeg_quality(end+1).subject = subject;
                eeg_quality(end).channels = channels;
                eeg_quality(end).srate = srate;
                eeg_quality(end).duration = duration;
                eeg_quality(end).epochs = epochs;
                eeg_quality(end).rejected_pct = rejected_pct;
                eeg_quality(end).snr = snr;
                eeg_quality(end).status = status;
                
            catch e
                warning('Error processing EEG file %s: %s', eeg_files(f).name, e.message);
            end
        end
        
        fprintf(fid, '</table>\n');
        fprintf(fid, '</div>\n');
        
        % Create EEG quality charts
        try
            % SNR distribution
            snr_values = [eeg_quality.snr];
            snr_subjects = {eeg_quality.subject};
            
            figure('Position', [100, 100, 800, 400]);
            
            % Bar chart of SNR by subject
            bar(snr_values);
            hold on;
            
            % Add reference lines
            plot([0, length(snr_values)+1], [10, 10], 'g--', 'LineWidth', 2);
            plot([0, length(snr_values)+1], [5, 5], 'y--', 'LineWidth', 2);
            
            xlabel('Subject');
            ylabel('SNR (dB)');
            title('EEG Signal-to-Noise Ratio by Subject');
            xticks(1:length(snr_values));
            xticklabels(snr_subjects);
            xtickangle(45);
            
            % Add legend
            legend({'SNR', 'Good (>10 dB)', 'Acceptable (>5 dB)'}, 'Location', 'best');
            
            % Save figure
            snr_chart_file = [task '_eeg_snr_chart.png'];
            saveas(gcf, fullfile(task_dir, snr_chart_file));
            close;
            
            % Add chart to report
            fprintf(fid, '<div class="chart-container">\n');
            fprintf(fid, '<div class="chart">\n');
            fprintf(fid, '<h3>EEG Signal Quality</h3>\n');
            fprintf(fid, '<img src="%s" alt="EEG SNR Chart" style="width:100%%">\n', snr_chart_file);
            fprintf(fid, '</div>\n');
        catch e
            warning('Error creating EEG quality charts: %s', e.message);
        end
        
        % Process fNIRS datasets
        fprintf(fid, '<div class="section">\n');
        fprintf(fid, '<h2>fNIRS Data Quality</h2>\n');
        
        % Check if we can access a sample dataset to extract quality metrics
        try
            % Find one NIRS dataset to analyze
            sample_dirs = dir(fullfile(nir_extn_path, '..', '..', 'temp', task, '**', '*.wl1'));
            
            if ~isempty(sample_dirs)
                % Load first dataset
                nirs_folder = sample_dirs(1).folder;
                data_raw = nirs.io.loadNIRx(nirs_folder);
                
                % Extract some info
                fprintf(fid, '<p>fNIRS Data Information:</p>\n');
                fprintf(fid, '<ul>\n');
                fprintf(fid, '<li>Number of sources: %d</li>\n', data_raw.probe.nsrc);
                fprintf(fid, '<li>Number of detectors: %d</li>\n', data_raw.probe.ndet);
                fprintf(fid, '<li>Sampling rate: %.2f Hz</li>\n', data_raw.Fs);
                fprintf(fid, '<li>Number of channels: %d</li>\n', size(data_raw.data, 2));
                fprintf(fid, '</ul>\n');
                
                % Add some quality metrics about fNIRS data
                quality_metrics = struct();
                
                % Check SNR across subjects
                fnirs_quality = [];
                
                for s = 1:length(fnirs_files)
                    try
                        % Read CSV data
                        fnirs_data = readtable(fullfile(fnirs_files(s).folder, fnirs_files(s).name));
                        
                        % Extract subject ID
                        subject = extractSubjectID(fnirs_files(s).name);
                        
                        % Find columns with hbo data
                        fnirs_cols = fnirs_data.Properties.VariableNames;
                        hbo_cols = find(contains(fnirs_cols, 'hbo'));
                        hbr_cols = find(contains(fnirs_cols, 'hbr'));
                        
                        % Calculate quality metrics if possible
                        if ~isempty(hbo_cols) && ~isempty(hbr_cols)
                            % Calculate coefficient of variation (CV) as quality indicator
                            hbo_values = table2array(fnirs_data(:,hbo_cols));
                            hbr_values = table2array(fnirs_data(:,hbr_cols));
                            
                            if ~isempty(hbo_values) && ~isempty(hbr_values)
                                hbo_cv = std(hbo_values, 'all', 'omitnan') / mean(hbo_values, 'all', 'omitnan');
                                hbr_cv = std(hbr_values, 'all', 'omitnan') / mean(hbr_values, 'all', 'omitnan');
                                
                                % Store quality info
                                fnirs_quality(end+1).subject = subject;
                                fnirs_quality(end).hbo_cv = hbo_cv;
                                fnirs_quality(end).hbr_cv = hbr_cv;
                                
                                % Determine status
                                if hbo_cv < 0.5 && hbr_cv < 0.5
                                    fnirs_quality(end).status = 'good';
                                elseif hbo_cv < 1 && hbr_cv < 1
                                    fnirs_quality(end).status = 'acceptable';
                                else
                                    fnirs_quality(end).status = 'poor';
                                end
                            end
                        end
                    catch e
                        warning('Error processing fNIRS file %s: %s', fnirs_files(s).name, e.message);
                    end
                end
                
                % Create fNIRS quality table
                fprintf(fid, '<table>\n');
                fprintf(fid, '<tr><th>Subject</th><th>HbO CV</th><th>HbR CV</th><th>Status</th></tr>\n');
                
                for q = 1:length(fnirs_quality)
                    % Determine status class
                    if strcmp(fnirs_quality(q).status, 'good')
                        status_label = '<span class="good">Good</span>';
                    elseif strcmp(fnirs_quality(q).status, 'acceptable')
                        status_label = '<span class="warning">Acceptable</span>';
                    else
                        status_label = '<span class="error">Poor</span>';
                    end
                    
                    fprintf(fid, '<tr><td>%s</td><td>%.2f</td><td>%.2f</td><td>%s</td></tr>\n', ...
                        fnirs_quality(q).subject, fnirs_quality(q).hbo_cv, fnirs_quality(q).hbr_cv, status_label);
                end
                
                fprintf(fid, '</table>\n');
                
                % Create NIRS quality chart
                try
                    % Extract values
                    hbo_cv = [fnirs_quality.hbo_cv];
                    hbr_cv = [fnirs_quality.hbr_cv];
                    subjects = {fnirs_quality.subject};
                    
                    figure('Position', [100, 100, 800, 400]);
                    
                    % Bar chart of CV by subject
                    bar([hbo_cv', hbr_cv']);
                    hold on;
                    
                    % Add reference lines
                    plot([0, length(subjects)+1], [0.5, 0.5], 'g--', 'LineWidth', 2);
                    plot([0, length(subjects)+1], [1, 1], 'y--', 'LineWidth', 2);
                    
                    xlabel('Subject');
                    ylabel('Coefficient of Variation');
                    title('fNIRS Signal Quality by Subject');
                    xticks(1:length(subjects));
                    xticklabels(subjects);
                    xtickangle(45);
                    legend({'HbO', 'HbR', 'Good (<0.5)', 'Acceptable (<1)'}, 'Location', 'best');
                    
                    % Save figure
                    fnirs_chart_file = [task '_fnirs_quality_chart.png'];
                    saveas(gcf, fullfile(task_dir, fnirs_chart_file));
                    close;
                    
                    % Add chart to report
                    fprintf(fid, '<div class="chart">\n');
                    fprintf(fid, '<h3>fNIRS Signal Quality</h3>\n');
                    fprintf(fid, '<img src="%s" alt="fNIRS Quality Chart" style="width:100%%">\n', fnirs_chart_file);
                    fprintf(fid, '</div>\n');
                    fprintf(fid, '</div>\n');
                catch e
                    warning('Error creating fNIRS quality chart: %s', e.message);
                end
            end
        catch e
            warning('Error processing fNIRS data for quality report: %s', e.message);
            fprintf(fid, '<p>Could not access raw fNIRS data for detailed quality metrics.</p>\n');
        end
        
        fprintf(fid, '</div>\n');
        
        % Multimodal integration summary
        fprintf(fid, '<div class="section">\n');
        fprintf(fid, '<h2>Multimodal Integration</h2>\n');
        
        % Check for multimodal connectivity files
        multimodal_files = dir(fullfile(eeg_extn_path, 'multimodal', task, '*_multimodal_conn.csv'));
        
        if ~isempty(multimodal_files)
            fprintf(fid, '<p>Successfully generated multimodal connectivity for %d subjects.</p>\n', length(multimodal_files));
            
            % Add information about each subject's connectivity
            fprintf(fid, '<table>\n');
            fprintf(fid, '<tr><th>Subject</th><th>Mean Connectivity</th><th>Connectivity Range</th></tr>\n');
            
            all_conn_means = [];
            
            for m = 1:length(multimodal_files)
                try
                    % Load connectivity matrix
                    [~, file_base, ~] = fileparts(multimodal_files(m).name);
                    parts = strsplit(file_base, '_');
                    subject = parts{1};
                    
                    conn_matrix = csvread(fullfile(multimodal_files(m).folder, multimodal_files(m).name));
                    
                    % Calculate metrics
                    conn_mean = mean(conn_matrix(:));
                    conn_min = min(conn_matrix(:));
                    conn_max = max(conn_matrix(:));
                    
                    % Determine status based on connectivity strength
                    if abs(conn_mean) > 0.4
                        status_class = 'good';
                    elseif abs(conn_mean) > 0.2
                        status_class = 'warning';
                    else
                        status_class = 'error';
                    end
                    
                    fprintf(fid, '<tr><td>%s</td><td class="%s">%.3f</td><td>[%.3f, %.3f]</td></tr>\n', ...
                        subject, status_class, conn_mean, conn_min, conn_max);
                    
                    all_conn_means(end+1) = conn_mean;
                catch e
                    warning('Error processing multimodal file %s: %s', multimodal_files(m).name, e.message);
                end
            end
            
            fprintf(fid, '</table>\n');
            
            % Plot histogram of connectivity
            if ~isempty(all_conn_means)
                figure;
                histogram(all_conn_means, 10, 'FaceColor', 'b', 'EdgeColor', 'k');
                xlabel('Mean Connectivity Strength');
                ylabel('Number of Subjects');
                title('Distribution of EEG-fNIRS Connectivity');
                grid on;
                
                % Save figure
                conn_hist_file = [task '_connectivity_hist.png'];
                saveas(gcf, fullfile(task_dir, conn_hist_file));
                close;
                
                fprintf(fid, '<div style="text-align:center">\n');
                fprintf(fid, '<h3>Distribution of EEG-fNIRS Connectivity</h3>\n');
                fprintf(fid, '<img src="%s" alt="Connectivity Histogram" style="width:60%%">\n', conn_hist_file);
                fprintf(fid, '</div>\n');
            end
        else
            fprintf(fid, '<p>No multimodal connectivity files found.</p>\n');
        end
        
        fprintf(fid, '</div>\n');
        
        % Summary and recommendations
        fprintf(fid, '<div class="section">\n');
        fprintf(fid, '<h2>Summary and Recommendations</h2>\n');
        
        % Count quality issues
        eeg_good = sum(strcmp({eeg_quality.status}, 'good'));
        eeg_acceptable = sum(strcmp({eeg_quality.status}, 'acceptable'));
        eeg_poor = sum(strcmp({eeg_quality.status}, 'poor'));
        
        fprintf(fid, '<h3>EEG Data Quality:</h3>\n');
        fprintf(fid, '<ul>\n');
        fprintf(fid, '<li>Good: %d subjects (%.1f%%)</li>\n', eeg_good, 100*eeg_good/length(eeg_quality));
        fprintf(fid, '<li>Acceptable: %d subjects (%.1f%%)</li>\n', eeg_acceptable, 100*eeg_acceptable/length(eeg_quality));
        fprintf(fid, '<li>Poor: %d subjects (%.1f%%)</li>\n', eeg_poor, 100*eeg_poor/length(eeg_quality));
        fprintf(fid, '</ul>\n');
        
        % Add recommendations
        fprintf(fid, '<h3>Recommendations:</h3>\n');
        fprintf(fid, '<ul>\n');
        
        if eeg_poor > 0
            fprintf(fid, '<li>Consider excluding subjects with poor EEG quality from analysis.</li>\n');
        end
        
        if length(common_subjects) < length(eeg_subjects) || length(common_subjects) < length(fnirs_subjects)
            fprintf(fid, '<li>Some subjects have only EEG or only fNIRS data. Consider collecting both modalities for all subjects in future studies.</li>\n');
        end
        
        fprintf(fid, '<li>For optimal multimodal analysis, ensure recording durations and sampling rates are compatible between modalities.</li>\n');
        fprintf(fid, '</ul>\n');
        
        fprintf(fid, '</div>\n');
        
        % Close HTML file
        fprintf(fid, '</body>\n</html>');
        fclose(fid);
        
        disp(['Generated quality report for ' task ': ' html_file]);
    end
    
    disp('=== Quality Control and Reporting Completed ===');
end

%% Helper function to extract subject ID from filename
function subject_id = extractSubjectID(filename)
    % Extract subject ID from filename
    % This is a simple implementation - might need customization
    parts = strsplit(filename, '_');
    if length(parts) > 0
        subject_id = parts{1};
    else
        subject_id = filename;
    end
end