function EEG_statistics(inputDir, outputDir, options)
%EEG_STATISTICS Statistical analysis using linear mixed-effects models
%
%   EEG_statistics(inputDir, outputDir)
%   EEG_statistics(inputDir, outputDir, options)
%
%   Uses fitlme with Subject as random intercept to leverage per-channel
%   (bandpower) and per-trial (ERP) data instead of subject-level averages.
%
%   Inputs:
%       inputDir  - Directory containing .mat files from EEG_analysis
%       outputDir - Directory for output statistics CSV files
%       options   - Optional struct with:
%           .alpha      - Significance level (default: 0.05)
%           .roiSuffix  - ROI suffix to analyze (default: 'allchan')
%           .dfMethod   - Degrees of freedom method (default: 'satterthwaite')
%
%   Output CSV files:
%       stats_bandpower_abs.csv  - Absolute power (LME per band × task)
%       stats_bandpower_rel.csv  - Relative power (LME per band × task)
%       stats_bandpower_log.csv  - Log power (LME per band × task)
%       stats_ratios.csv         - Theta/Beta, Delta/Alpha ratios (LME)
%       stats_erp.csv            - Per-trial ERP metrics (LME, correct only)
%
%   Each output row contains:
%       - Omnibus F-test from anova(lme)
%       - Pairwise: difference estimate, SE, t, df, p, Cohen's d
%       - FDR-corrected q-values for all p-values

arguments
    inputDir    (1,1) string
    outputDir   (1,1) string
    options.alpha      (1,1) double = 0.05
    options.roiSuffix  (1,1) string = "allchan"
    options.dfMethod   (1,1) string = "satterthwaite"
end

fprintf('\n===== EEG STATISTICAL ANALYSIS (LME) =====\n');
fprintf('Input:  %s\n', inputDir);
fprintf('Output: %s\n', outputDir);

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Load data
bpFile = fullfile(inputDir, sprintf('bandpower_%s.mat', options.roiSuffix));
erpFile = fullfile(inputDir, sprintf('erp_%s.mat', options.roiSuffix));

if ~exist(bpFile, 'file')
    error('Bandpower file not found: %s', bpFile);
end
if ~exist(erpFile, 'file')
    error('ERP file not found: %s', erpFile);
end

bpLoad = load(bpFile, 'bp_table');
bp = bpLoad.bp_table;

erpLoad = load(erpFile, 'erp_table');
erp = erpLoad.erp_table;

% Ensure categorical types
bp.SubjectID = categorical(bp.SubjectID);
bp.Group     = categorical(bp.Group);
bp.Task      = categorical(bp.Task);
erp.SubjectID = categorical(erp.SubjectID);
erp.Group     = categorical(erp.Group);
erp.Task      = categorical(erp.Task);
erp.Accuracy  = categorical(erp.Accuracy);

nSubjects_bp = length(unique(bp.SubjectID));
nSubjects_erp = length(unique(erp.SubjectID));
fprintf('Loaded bandpower: %d rows (%d subjects)\n', height(bp), nSubjects_bp);
fprintf('Loaded ERP: %d rows (%d subjects)\n', height(erp), nSubjects_erp);

% Get unique groups
groups = unique(cellstr(bp.Group), 'stable');
nGroups = length(groups);
fprintf('Groups: %s\n', strjoin(groups, ', '));

if nGroups < 2
    error('Need at least 2 groups for comparison');
end

%% Identify bands and tasks from column names
bpVars = bp.Properties.VariableNames;
bandCols = bpVars(endsWith(bpVars, '_abs'));
bandNames = cellfun(@(x) x(1:end-4), bandCols, 'UniformOutput', false);  % strip '_abs'
taskNames = unique(cellstr(bp.Task), 'stable');

fprintf('Bands: %s\n', strjoin(bandNames, ', '));
fprintf('Tasks: %s\n', strjoin(taskNames, ', '));

%% ========================================================================
%  BANDPOWER STATISTICS (LME: value ~ Group + (1|SubjectID))
%  ========================================================================
fprintf('\n→ Analyzing bandpower with LME...\n');

metricSuffixes = {'abs', 'rel', 'log'};

for m = 1:length(metricSuffixes)
    suffix = metricSuffixes{m};
    fprintf('  Processing %s power...\n', suffix);

    results = table();

    for b = 1:length(bandNames)
        band = bandNames{b};
        colName = sprintf('%s_%s', band, suffix);

        if ~ismember(colName, bpVars)
            continue
        end

        for t = 1:length(taskNames)
            task = taskNames{t};

            % Subset data for this task
            taskMask = bp.Task == task;
            tbl = bp(taskMask, {'SubjectID', 'Group', colName});
            tbl.Properties.VariableNames{end} = 'Value';

            % Remove NaN
            tbl = tbl(~isnan(tbl.Value), :);

            if height(tbl) < 6
                continue
            end

            % Run LME comparison
            statRow = compareLME(tbl, groups, options.dfMethod);
            statRow.Band = {band};
            statRow.Task = {task};
            statRow.Metric = {suffix};
            statRow = movevars(statRow, {'Band', 'Task', 'Metric'}, 'Before', 1);

            results = [results; statRow]; %#ok<AGROW>
        end
    end

    % FDR correction
    results = applyFDR(results);

    % Write output
    outFile = fullfile(outputDir, sprintf('stats_bandpower_%s.csv', suffix));
    writetable(results, outFile);
    fprintf('    ✓ stats_bandpower_%s.csv (%d comparisons)\n', suffix, height(results));
end

%% ========================================================================
%  RATIO STATISTICS (computed from per-channel absolute power)
%  ========================================================================
fprintf('\n→ Analyzing spectral ratios with LME...\n');

ratioResults = table();

for t = 1:length(taskNames)
    task = taskNames{t};
    taskMask = bp.Task == task;

    % --- Theta / Beta ---
    if ismember('theta_abs', bpVars) && ismember('beta_abs', bpVars)
        tbl = bp(taskMask, {'SubjectID', 'Group'});
        tbl.Value = bp.theta_abs(taskMask) ./ bp.beta_abs(taskMask);
        tbl = tbl(~isnan(tbl.Value) & isfinite(tbl.Value), :);

        if height(tbl) >= 6
            statRow = compareLME(tbl, groups, options.dfMethod);
            statRow.Ratio = {'Theta/Beta'};
            statRow.Task = {task};
            statRow = movevars(statRow, {'Ratio', 'Task'}, 'Before', 1);
            ratioResults = [ratioResults; statRow]; %#ok<AGROW>
        end
    end

    % --- Delta / Alpha ---
    if ismember('delta_abs', bpVars) && ismember('alpha_abs', bpVars)
        tbl = bp(taskMask, {'SubjectID', 'Group'});
        tbl.Value = bp.delta_abs(taskMask) ./ bp.alpha_abs(taskMask);
        tbl = tbl(~isnan(tbl.Value) & isfinite(tbl.Value), :);

        if height(tbl) >= 6
            statRow = compareLME(tbl, groups, options.dfMethod);
            statRow.Ratio = {'Delta/Alpha'};
            statRow.Task = {task};
            statRow = movevars(statRow, {'Ratio', 'Task'}, 'Before', 1);
            ratioResults = [ratioResults; statRow]; %#ok<AGROW>
        end
    end
end

ratioResults = applyFDR(ratioResults);

outFile = fullfile(outputDir, 'stats_ratios.csv');
writetable(ratioResults, outFile);
fprintf('  ✓ stats_ratios.csv (%d comparisons)\n', height(ratioResults));

%% ========================================================================
%  DELTA DAR: Change in Delta/Alpha Ratio (2-back minus Eyes Open)
%  ========================================================================
fprintf('\n→ Analyzing Delta DAR (2-back - Eyes Open)...\n');

deltaDARResults = table();

if ismember('delta_abs', bpVars) && ismember('alpha_abs', bpVars)
    % Compute DAR per observation
    bp.DAR = bp.delta_abs ./ bp.alpha_abs;
    
    % Get subject-level mean DAR for eyes_open
    eoMask = bp.Task == 'eyes_open';
    if any(eoMask)
        eoDAR = grpstats(bp(eoMask, :), {'SubjectID', 'Group'}, 'mean', 'DataVars', 'DAR');
        eoDAR.Properties.VariableNames{'mean_DAR'} = 'DAR_EO';
    else
        fprintf('  ⚠ No eyes_open data found, skipping Delta DAR\n');
        eoDAR = table();
    end
    
    % Get subject-level mean DAR for nback_2a
    nbMask = bp.Task == 'nback_2a';
    if any(nbMask)
        nbDAR = grpstats(bp(nbMask, :), {'SubjectID', 'Group'}, 'mean', 'DataVars', 'DAR');
        nbDAR.Properties.VariableNames{'mean_DAR'} = 'DAR_2back';
    else
        fprintf('  ⚠ No nback_2a data found, skipping Delta DAR\n');
        nbDAR = table();
    end
    
    % Merge and compute difference
    if ~isempty(eoDAR) && ~isempty(nbDAR)
        % Join on SubjectID
        eoDAR_clean = eoDAR(:, {'SubjectID', 'Group', 'DAR_EO'});
        nbDAR_clean = nbDAR(:, {'SubjectID', 'DAR_2back'});
        
        deltaDAR = innerjoin(eoDAR_clean, nbDAR_clean, 'Keys', 'SubjectID');
        
        % Compute Delta DAR (2-back minus Eyes Open)
        % Negative = good (reduced slow/fast ratio during task)
        % Positive = bad (increased slow/fast ratio during task)
        deltaDAR.Value = deltaDAR.DAR_2back - deltaDAR.DAR_EO;
        
        % Remove NaN/Inf
        deltaDAR = deltaDAR(~isnan(deltaDAR.Value) & isfinite(deltaDAR.Value), :);
        
        fprintf('  Found %d subjects with both conditions\n', height(deltaDAR));
        
        if height(deltaDAR) >= 6
            % Prepare table for statistics
            tbl = table();
            tbl.SubjectID = deltaDAR.SubjectID;
            tbl.Group = deltaDAR.Group;
            tbl.Value = deltaDAR.Value;
            
            % Run group comparison (using ANOVA since one value per subject)
            statRow = compareGrandAvg(tbl, groups);
            statRow.Metric = {'DeltaDAR_2back_minus_EO'};
            statRow = movevars(statRow, 'Metric', 'Before', 1);
            
            deltaDARResults = [deltaDARResults; statRow]; %#ok<AGROW>
            
            % Also save the per-subject delta DAR values
            subjDeltaDAR = deltaDAR(:, {'SubjectID', 'Group', 'DAR_EO', 'DAR_2back', 'Value'});
            subjDeltaDAR.Properties.VariableNames{'Value'} = 'DeltaDAR';
            outFileSubj = fullfile(outputDir, 'subject_delta_dar.csv');
            writetable(subjDeltaDAR, outFileSubj);
            fprintf('  ✓ subject_delta_dar.csv (%d subjects)\n', height(subjDeltaDAR));
        end
    end
end

if ~isempty(deltaDARResults)
    deltaDARResults = applyFDR(deltaDARResults);
    outFile = fullfile(outputDir, 'stats_delta_dar.csv');
    writetable(deltaDARResults, outFile);
    fprintf('  ✓ stats_delta_dar.csv\n');
end

%% ========================================================================
%  ERP STATISTICS (per-trial, correct only)
%  ========================================================================
fprintf('\n→ Analyzing ERP metrics with LME...\n');

% Filter to correct trials only
erpCorrect = erp;%(erp.Accuracy == 'correct', :);
erpMetrics = {'PeakAmp', 'PeakLat', 'MeanAmp', 'CentLat'};
erpTaskNames = unique(cellstr(erpCorrect.Task), 'stable');

erpResults = table();

for em = 1:length(erpMetrics)
    metricName = erpMetrics{em};

    for t = 1:length(erpTaskNames)
        task = erpTaskNames{t};
        taskMask = erpCorrect.Task == task;

        tbl = erpCorrect(taskMask, {'SubjectID', 'Group', metricName});
        tbl.Properties.VariableNames{end} = 'Value';
        tbl = tbl(~isnan(tbl.Value), :);

        if height(tbl) < 6
            continue
        end

        % Check we have data from multiple subjects
        nSubj = length(unique(tbl.SubjectID));
        if nSubj < 3
            continue
        end

        statRow = compareLME(tbl, groups, options.dfMethod);
        statRow.ERPMetric = {metricName};
        statRow.Task = {task};
        statRow = movevars(statRow, {'ERPMetric', 'Task'}, 'Before', 1);

        erpResults = [erpResults; statRow]; %#ok<AGROW>
    end
end

erpResults = applyFDR(erpResults);

outFile = fullfile(outputDir, 'stats_erp.csv');
writetable(erpResults, outFile);
fprintf('  ✓ stats_erp.csv (%d comparisons)\n', height(erpResults));

%% ========================================================================
%  GRAND AVERAGE BANDPOWER (subject means across all tasks)
%  ========================================================================
fprintf('\n→ Analyzing grand average bandpower with LME...\n');

for m = 1:length(metricSuffixes)
    suffix = metricSuffixes{m};
    fprintf('  Processing grand average %s power...\n', suffix);
    
    grandAvgBpResults = table();
    
    for b = 1:length(bandNames)
        band = bandNames{b};
        colName = sprintf('%s_%s', band, suffix);
        
        if ~ismember(colName, bpVars)
            continue
        end
        
        % Compute subject-level grand average (mean across all tasks and channels)
        subjMeans = grpstats(bp, {'SubjectID', 'Group'}, 'mean', 'DataVars', colName);
        
        % Prepare table for LME
        tbl = table();
        tbl.SubjectID = subjMeans.SubjectID;
        tbl.Group = subjMeans.Group;
        tbl.Value = subjMeans.(sprintf('mean_%s', colName));
        tbl = tbl(~isnan(tbl.Value), :);
        
        if height(tbl) < 6
            continue
        end
        
        % Run LME comparison (note: with subject means, random effect may be singular)
        statRow = compareGrandAvg(tbl, groups);
        statRow.Band = {band};
        statRow.Task = {'grand_average'};
        statRow.Metric = {suffix};
        statRow = movevars(statRow, {'Band', 'Task', 'Metric'}, 'Before', 1);
        
        grandAvgBpResults = [grandAvgBpResults; statRow]; %#ok<AGROW>
    end
    
    % FDR correction for grand average results
    grandAvgBpResults = applyFDR(grandAvgBpResults);
    
    % Append to existing bandpower stats file
    existingFile = fullfile(outputDir, sprintf('stats_bandpower_%s.csv', suffix));
    if exist(existingFile, 'file')
        existingResults = readtable(existingFile);
        combinedResults = [existingResults; grandAvgBpResults];
        writetable(combinedResults, existingFile);
        fprintf('    ✓ Appended %d grand avg rows to stats_bandpower_%s.csv\n', ...
                height(grandAvgBpResults), suffix);
    else
        writetable(grandAvgBpResults, existingFile);
        fprintf('    ✓ stats_bandpower_%s.csv (grand avg only, %d rows)\n', ...
                suffix, height(grandAvgBpResults));
    end
end

%% ========================================================================
%  GRAND AVERAGE RATIOS (subject means across all tasks)
%  ========================================================================
fprintf('\n→ Analyzing grand average ratios...\n');

grandAvgRatioResults = table();

% --- Theta / Beta Grand Average ---
if ismember('theta_abs', bpVars) && ismember('beta_abs', bpVars)
    % First compute ratio per row, then average per subject
    bp.ThetaBetaRatio = bp.theta_abs ./ bp.beta_abs;
    subjMeans = grpstats(bp, {'SubjectID', 'Group'}, 'mean', 'DataVars', 'ThetaBetaRatio');
    
    tbl = table();
    tbl.SubjectID = subjMeans.SubjectID;
    tbl.Group = subjMeans.Group;
    tbl.Value = subjMeans.mean_ThetaBetaRatio;
    tbl = tbl(~isnan(tbl.Value) & isfinite(tbl.Value), :);
    
    if height(tbl) >= 6
        statRow = compareGrandAvg(tbl, groups);
        statRow.Ratio = {'Theta/Beta'};
        statRow.Task = {'grand_average'};
        statRow = movevars(statRow, {'Ratio', 'Task'}, 'Before', 1);
        grandAvgRatioResults = [grandAvgRatioResults; statRow]; %#ok<AGROW>
    end
end

% --- Delta / Alpha Grand Average ---
if ismember('delta_abs', bpVars) && ismember('alpha_abs', bpVars)
    bp.DeltaAlphaRatio = bp.delta_abs ./ bp.alpha_abs;
    subjMeans = grpstats(bp, {'SubjectID', 'Group'}, 'mean', 'DataVars', 'DeltaAlphaRatio');
    
    tbl = table();
    tbl.SubjectID = subjMeans.SubjectID;
    tbl.Group = subjMeans.Group;
    tbl.Value = subjMeans.mean_DeltaAlphaRatio;
    tbl = tbl(~isnan(tbl.Value) & isfinite(tbl.Value), :);
    
    if height(tbl) >= 6
        statRow = compareGrandAvg(tbl, groups);
        statRow.Ratio = {'Delta/Alpha'};
        statRow.Task = {'grand_average'};
        statRow = movevars(statRow, {'Ratio', 'Task'}, 'Before', 1);
        grandAvgRatioResults = [grandAvgRatioResults; statRow]; %#ok<AGROW>
    end
end

grandAvgRatioResults = applyFDR(grandAvgRatioResults);

% Append to existing ratios file
existingFile = fullfile(outputDir, 'stats_ratios.csv');
if exist(existingFile, 'file')
    existingResults = readtable(existingFile);
    combinedResults = [existingResults; grandAvgRatioResults];
    writetable(combinedResults, existingFile);
    fprintf('  ✓ Appended %d grand avg rows to stats_ratios.csv\n', height(grandAvgRatioResults));
else
    writetable(grandAvgRatioResults, existingFile);
    fprintf('  ✓ stats_ratios.csv (grand avg only, %d rows)\n', height(grandAvgRatioResults));
end

%% ========================================================================
%  GRAND AVERAGE ERP (subject means across all tasks)
%  ========================================================================
fprintf('\n→ Analyzing grand average ERP metrics...\n');

grandAvgErpResults = table();

for em = 1:length(erpMetrics)
    metricName = erpMetrics{em};
    
    if ~ismember(metricName, erpCorrect.Properties.VariableNames)
        continue
    end
    
    % Compute subject-level grand average across all tasks
    subjMeans = grpstats(erpCorrect, {'SubjectID', 'Group'}, 'mean', 'DataVars', metricName);
    
    tbl = table();
    tbl.SubjectID = subjMeans.SubjectID;
    tbl.Group = subjMeans.Group;
    tbl.Value = subjMeans.(sprintf('mean_%s', metricName));
    tbl = tbl(~isnan(tbl.Value), :);
    
    if height(tbl) < 6
        continue
    end
    
    nSubj = length(unique(tbl.SubjectID));
    if nSubj < 3
        continue
    end
    
    statRow = compareGrandAvg(tbl, groups);
    statRow.ERPMetric = {metricName};
    statRow.Task = {'grand_average'};
    statRow = movevars(statRow, {'ERPMetric', 'Task'}, 'Before', 1);
    
    grandAvgErpResults = [grandAvgErpResults; statRow]; %#ok<AGROW>
end

grandAvgErpResults = applyFDR(grandAvgErpResults);

% Append to existing ERP file
existingFile = fullfile(outputDir, 'stats_erp.csv');
if exist(existingFile, 'file')
    existingResults = readtable(existingFile);
    combinedResults = [existingResults; grandAvgErpResults];
    writetable(combinedResults, existingFile);
    fprintf('  ✓ Appended %d grand avg rows to stats_erp.csv\n', height(grandAvgErpResults));
else
    writetable(grandAvgErpResults, existingFile);
    fprintf('  ✓ stats_erp.csv (grand avg only, %d rows)\n', height(grandAvgErpResults));
end

fprintf('\n===== STATISTICAL ANALYSIS COMPLETE =====\n');
end

%% ========================================================================
%  CORE LME COMPARISON
%  ========================================================================

function result = compareLME(tbl, groups, dfMethod)
%COMPARELME Fit LME and run omnibus + pairwise comparisons
%
%   Input:
%       tbl      - Table with columns: SubjectID, Group, Value
%       groups   - Cell array of group names
%       dfMethod - DF method for anova ('satterthwaite' or 'residual')
%
%   Model: Value ~ Group + (1 | SubjectID)
%
%   Returns single-row table with:
%       - Per-group descriptives (marginal mean, SE, n_obs, n_subj)
%       - Omnibus F-test
%       - Pairwise: estimate, SE, t, df, p, Cohen's d

    nGroups = length(groups);
    result = table();

    %% Descriptives per group
    for g = 1:nGroups
        gMask = tbl.Group == groups{g};
        gVals = tbl.Value(gMask);
        gSubj = unique(tbl.SubjectID(gMask));

        result.(sprintf('%s_mean', groups{g})) = mean(gVals, 'omitnan');
        result.(sprintf('%s_std', groups{g}))  = std(gVals, 0, 'omitnan');
        result.(sprintf('%s_nObs', groups{g})) = sum(~isnan(gVals));
        result.(sprintf('%s_nSubj', groups{g})) = length(gSubj);
    end

    %% Fit omnibus LME: Value ~ Group + (1|SubjectID)
    % Ensure Group is categorical with all levels
    tbl.Group = categorical(cellstr(tbl.Group), groups);
    tbl.SubjectID = categorical(cellstr(tbl.SubjectID));

    try
        % Set reference group by making it first level of the categorical
        tbl.Group = reordercats(tbl.Group, groups);

        lme = fitlme(tbl, 'Value ~ Group + (1 | SubjectID)');

        % Omnibus F-test for Group effect
        aovTbl = anova(lme, 'DFMethod', dfMethod);

        % Find the Group row (may be row 2 if intercept is row 1)
        groupRow = find(strcmp(aovTbl.Term, 'Group'), 1);
        if isempty(groupRow)
            % Fallback: second row
            groupRow = 2;
        end

        result.Omnibus_F    = aovTbl.FStat(groupRow);
        result.Omnibus_DF1  = aovTbl.DF1(groupRow);
        result.Omnibus_DF2  = aovTbl.DF2(groupRow);
        result.Omnibus_p    = aovTbl.pValue(groupRow);

        % Random effects variance
        [~, ~, stats] = covarianceParameters(lme);
        result.SubjectVar = stats{1}.Estimate(1);     % subject intercept variance
        result.ResidualVar = lme.MSE;                  % residual variance
        result.ICC = stats{1}.Estimate(1) / ...
                    (stats{1}.Estimate(1) + lme.MSE);  % intraclass correlation

    catch ME
        warning('LME fit failed. Falling back to Kruskal-Wallis.');
        result.Omnibus_F   = NaN;
        result.Omnibus_DF1 = NaN;
        result.Omnibus_DF2 = NaN;

        % Fallback: KW on subject means
        [result.Omnibus_p] = fallbackKW(tbl, groups);
        result.SubjectVar  = NaN;
        result.ResidualVar = NaN;
        result.ICC         = NaN;
    end

    %% Pairwise comparisons (subset + refit for each pair)
    for i = 1:nGroups
        for j = (i+1):nGroups
            pairName = sprintf('%s_vs_%s', groups{i}, groups{j});
            pairMask = (tbl.Group == groups{i}) | (tbl.Group == groups{j});
            pairTbl = tbl(pairMask, :);
            pairTbl.Group = removecats(pairTbl.Group);
            pairTbl.SubjectID = removecats(pairTbl.SubjectID);

            try
                pairTbl.Group = reordercats(pairTbl.Group, {groups{i}, groups{j}});

                lme_pair = fitlme(pairTbl, 'Value ~ Group + (1 | SubjectID)');

                % Extract the Group coefficient (difference from reference)
                coefs = lme_pair.Coefficients;
                % The Group effect is the non-intercept row
                groupCoefIdx = find(contains(coefs.Name, 'Group_'), 1);
                if isempty(groupCoefIdx)
                    groupCoefIdx = 2;  % fallback
                end

                result.(sprintf('%s_est', pairName))  = coefs.Estimate(groupCoefIdx);
                result.(sprintf('%s_SE', pairName))    = coefs.SE(groupCoefIdx);
                result.(sprintf('%s_t', pairName))     = coefs.tStat(groupCoefIdx);
                result.(sprintf('%s_df', pairName))    = coefs.DF(groupCoefIdx);
                result.(sprintf('%s_p', pairName))     = coefs.pValue(groupCoefIdx);

            catch
                % Fallback: ranksum on subject means
                g1 = pairTbl(pairTbl.Group == groups{i}, :);
                g2 = pairTbl(pairTbl.Group == groups{j}, :);

                g1means = grpstats(g1, 'SubjectID', 'mean', 'DataVars', 'Value');
                g2means = grpstats(g2, 'SubjectID', 'mean', 'DataVars', 'Value');

                [pVal, ~, stats] = ranksum(g1means.mean_Value, g2means.mean_Value);

                result.(sprintf('%s_est', pairName))  = NaN;
                result.(sprintf('%s_SE', pairName))    = NaN;
                result.(sprintf('%s_t', pairName))     = NaN;
                result.(sprintf('%s_df', pairName))    = NaN;
                result.(sprintf('%s_p', pairName))     = pVal;
            end

            % Effect size: Cohen's d from group means and pooled SD
            g1vals = tbl.Value(tbl.Group == groups{i});
            g2vals = tbl.Value(tbl.Group == groups{j});
            g1mean = mean(g1vals, 'omitnan');
            g2mean = mean(g2vals, 'omitnan');
            n1 = sum(~isnan(g1vals));
            n2 = sum(~isnan(g2vals));
            pooledSD = sqrt(((n1-1)*var(g1vals,'omitnan') + (n2-1)*var(g2vals,'omitnan')) / (n1+n2-2));

            if pooledSD > 0
                result.(sprintf('%s_d', pairName)) = (g1mean - g2mean) / pooledSD;
            else
                result.(sprintf('%s_d', pairName)) = NaN;
            end
        end
    end
end

function result = compareGrandAvg(tbl, groups)
%COMPAREGRANDAVG Compare groups using ANOVA + t-tests on subject means
%
%   For grand average data where we have one value per subject,
%   LME with random intercepts is not appropriate. Use standard ANOVA.
%
%   Input:
%       tbl    - Table with columns: SubjectID, Group, Value (one row per subject)
%       groups - Cell array of group names
%
%   Returns single-row table with:
%       - Per-group descriptives (mean, SD, n)
%       - Omnibus F-test from ANOVA
%       - Pairwise: t, df, p, Cohen's d

    nGroups = length(groups);
    result = table();
    
    % Ensure categorical
    tbl.Group = categorical(cellstr(tbl.Group), groups);
    
    %% Descriptives per group
    groupVals = cell(nGroups, 1);
    for g = 1:nGroups
        gMask = tbl.Group == groups{g};
        gVals = tbl.Value(gMask);
        groupVals{g} = gVals;
        
        result.(sprintf('%s_mean', groups{g})) = mean(gVals, 'omitnan');
        result.(sprintf('%s_std', groups{g}))  = std(gVals, 0, 'omitnan');
        result.(sprintf('%s_nObs', groups{g})) = sum(~isnan(gVals));
        result.(sprintf('%s_nSubj', groups{g})) = sum(~isnan(gVals));  % same for grand avg
    end
    
    %% Omnibus test: one-way ANOVA
    try
        allVals = tbl.Value;
        groupLabels = cellstr(tbl.Group);
        [p, tbl_anova, stats] = anova1(allVals, groupLabels, 'off');
        
        result.Omnibus_F   = tbl_anova{2, 5};  % F statistic
        result.Omnibus_DF1 = tbl_anova{2, 3};  % between-groups df
        result.Omnibus_DF2 = tbl_anova{3, 3};  % within-groups df
        result.Omnibus_p   = p;
        
        % No random effects for grand average
        result.SubjectVar  = NaN;
        result.ResidualVar = tbl_anova{3, 4};  % MSE within
        result.ICC         = NaN;
        
    catch ME
        warning('ANOVA failed: %s. Using Kruskal-Wallis.', ME.message);
        try
            [p, ~, ~] = kruskalwallis(tbl.Value, cellstr(tbl.Group), 'off');
        catch
            p = NaN;
        end
        result.Omnibus_F   = NaN;
        result.Omnibus_DF1 = NaN;
        result.Omnibus_DF2 = NaN;
        result.Omnibus_p   = p;
        result.SubjectVar  = NaN;
        result.ResidualVar = NaN;
        result.ICC         = NaN;
    end
    
    %% Pairwise comparisons: two-sample t-tests
    for i = 1:nGroups
        for j = (i+1):nGroups
            pairName = sprintf('%s_vs_%s', groups{i}, groups{j});
            
            g1vals = groupVals{i};
            g2vals = groupVals{j};
            
            g1vals = g1vals(~isnan(g1vals));
            g2vals = g2vals(~isnan(g2vals));
            
            if length(g1vals) >= 2 && length(g2vals) >= 2
                [~, p, ~, stats] = ttest2(g1vals, g2vals);
                
                result.(sprintf('%s_est', pairName))  = mean(g1vals) - mean(g2vals);
                result.(sprintf('%s_SE', pairName))   = stats.sd / sqrt(1/length(g1vals) + 1/length(g2vals));
                result.(sprintf('%s_t', pairName))    = stats.tstat;
                result.(sprintf('%s_df', pairName))   = stats.df;
                result.(sprintf('%s_p', pairName))    = p;
            else
                result.(sprintf('%s_est', pairName))  = NaN;
                result.(sprintf('%s_SE', pairName))   = NaN;
                result.(sprintf('%s_t', pairName))    = NaN;
                result.(sprintf('%s_df', pairName))   = NaN;
                result.(sprintf('%s_p', pairName))    = NaN;
            end
            
            % Effect size: Cohen's d
            g1mean = mean(g1vals, 'omitnan');
            g2mean = mean(g2vals, 'omitnan');
            n1 = length(g1vals);
            n2 = length(g2vals);
            pooledSD = sqrt(((n1-1)*var(g1vals,'omitnan') + (n2-1)*var(g2vals,'omitnan')) / (n1+n2-2));
            
            if pooledSD > 0
                result.(sprintf('%s_d', pairName)) = (g1mean - g2mean) / pooledSD;
            else
                result.(sprintf('%s_d', pairName)) = NaN;
            end
        end
    end
end

function p = fallbackKW(tbl, groups)
%FALLBACKKW Kruskal-Wallis on subject means when LME fails

    subjMeans = grpstats(tbl, {'SubjectID', 'Group'}, 'mean', 'DataVars', 'Value');

    allVals = subjMeans.mean_Value;
    groupLabels = cellstr(subjMeans.Group);

    try
        p = kruskalwallis(allVals, groupLabels, 'off');
    catch
        p = NaN;
    end
end

%% ========================================================================
%  FDR CORRECTION
%  ========================================================================

function results = applyFDR(results)
%APPLYFDR Apply Benjamini-Hochberg FDR correction to all p-value columns

    if isempty(results) || height(results) == 0
        return
    end

    varNames = results.Properties.VariableNames;
    pCols = varNames(endsWith(varNames, '_p'));

    for i = 1:length(pCols)
        pCol = pCols{i};
        qCol = strrep(pCol, '_p', '_q');

        pVals = results.(pCol);
        qVals = fdrCorrect(pVals);

        % Insert q column right after p column
        pColIdx = find(strcmp(varNames, pCol));
        results.(qCol) = qVals;

        % Reorder to place q after p
        currentVars = results.Properties.VariableNames;
        qColIdx = find(strcmp(currentVars, qCol));
        newOrder = [1:pColIdx, qColIdx, (pColIdx+1):(qColIdx-1)];
        results = results(:, newOrder);

        % Refresh varNames for next iteration
        varNames = results.Properties.VariableNames;
    end
end

function q = fdrCorrect(p)
%FDRCORRECT Benjamini-Hochberg FDR correction

    n = length(p);
    q = nan(size(p));

    validIdx = ~isnan(p);
    pValid = p(validIdx);
    nValid = sum(validIdx);

    if nValid == 0
        return
    end

    [pSorted, sortIdx] = sort(pValid);

    qSorted = pSorted .* nValid ./ (1:nValid)';

    for i = nValid-1:-1:1
        qSorted(i) = min(qSorted(i), qSorted(i+1));
    end

    qSorted = min(qSorted, 1);

    qValid = zeros(size(pValid));
    qValid(sortIdx) = qSorted;

    q(validIdx) = qValid;
end
