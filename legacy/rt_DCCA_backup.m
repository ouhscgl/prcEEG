function [DCCC, DCCA, scaling_exponent] = rt_DCCA(data, scales, varargin)
% RT_DCCA - Real-time Detrended Cross-Correlation Analysis for EEG signals
%
% This function performs real-time detrended cross-correlation analysis
% on EEG signals, based on the algorithm by Kaposzta et al. (2022).
%
% Usage:
%   [DCCC, DCCA, scaling_exponent] = rt_DCCA(data, scales, 'Param1', Value1, ...)
%
% Inputs:
%   data    - EEG data matrix [channels x timepoints] 
%   scales  - Vector of scales for analysis (e.g., [8 16 32 64 128])
%
% Optional Parameters:
%   'WindowSize'     - Size of the sliding window (default: max(scales)*4)
%   'StepSize'       - Step size for analysis updates (default: min(scales))
%   'BandFrequency'  - Specific band for filtering before analysis:
%                       'delta' (1-4 Hz), 'theta' (4-8 Hz), 'alpha' (8-13 Hz),
%                       'beta' (13-30 Hz), 'gamma' (30-45 Hz), 'broadband' (default)
%   'SamplingRate'   - EEG sampling rate in Hz (default: 250)
%   'ArtifactMask'   - Binary mask of artifacts (1=artifact, 0=clean) [1 x timepoints]
%   'OutputMode'     - 'full' for all time points, 'final' for last time point only (default)
%   'ChannelPairs'   - Cell array of channel pairs to analyze, or mask matrix [ch x ch]
%
% Outputs:
%   DCCC            - Detrended cross-correlation coefficient [ch x ch x scales x timepoints]
%   DCCA            - Detrended cross-correlation analysis values [ch x ch x scales x timepoints]
%   scaling_exponent - Estimated scaling exponents [ch x ch x timepoints]
%
% Reference:
%   Kaposzta et al. (2022). Real-Time Algorithm for Detrended Cross-
%   Correlation Analysis of Long-Range Coupled Processes.
%   Frontiers in Physiology. 13. 817268. doi: 10.3389/fphys.2022.817268.

% Parse input parameters
p = inputParser;
addRequired(p, 'data', @isnumeric);
addRequired(p, 'scales', @isnumeric);
addParameter(p, 'WindowSize', [], @isnumeric);
addParameter(p, 'StepSize', [], @isnumeric);
addParameter(p, 'BandFrequency', 'broadband', @ischar);
addParameter(p, 'SamplingRate', 250, @isnumeric);
addParameter(p, 'ArtifactMask', [], @isnumeric);
addParameter(p, 'OutputMode', 'final', @ischar);
addParameter(p, 'ChannelPairs', [], @(x) isempty(x) || isnumeric(x) || iscell(x));
parse(p, data, scales, varargin{:});

% Extract parameters
band = p.Results.BandFrequency;
smax = max(scales);
smin = min(scales);
win_size = p.Results.WindowSize;
step_size = p.Results.StepSize;
fs = p.Results.SamplingRate;
artifact_mask = p.Results.ArtifactMask;
output_mode = p.Results.OutputMode;
channel_pairs = p.Results.ChannelPairs;

% Set defaults if not specified
if isempty(win_size)
    win_size = smax * 4;
end
if isempty(step_size)
    step_size = smin;
end
if isempty(artifact_mask)
    artifact_mask = zeros(1, size(data, 2));
end

% Apply band-pass filtering if specified
if ~strcmpi(band, 'broadband')
    data = bandpass_filter(data, band, fs);
end

% Ensure data is [channels x timepoints]
if size(data, 1) > size(data, 2)
    data = data';
end

% Initialize variables
[nch, N] = size(data);
ns = length(scales);
n_windows = floor((N - win_size) / step_size) + 1;

% Set up channel pair mask
if ~isempty(channel_pairs)
    if iscell(channel_pairs)
        % Convert cell of channel pair indices to a mask
        channel_mask = zeros(nch, nch);
        for i = 1:length(channel_pairs)
            channel_mask(channel_pairs{i}(1), channel_pairs{i}(2)) = 1;
            channel_mask(channel_pairs{i}(2), channel_pairs{i}(1)) = 1;  % Make it symmetric
        end
    else
        % Use provided mask directly
        channel_mask = channel_pairs;
    end
    
    % Ensure diagonal is included
    channel_mask = channel_mask | eye(nch);
else
    % Analyze all channel pairs
    channel_mask = ones(nch, nch);
end

% Initialize output variables based on output mode
if strcmpi(output_mode, 'full')
    DCCC = zeros(nch, nch, ns, n_windows);
    DCCA = zeros(nch, nch, ns, n_windows);
    scaling_exponent = zeros(nch, nch, n_windows);
else % 'final' mode - only keep the last result
    DCCC = zeros(nch, nch, ns);
    DCCA = zeros(nch, nch, ns);
    scaling_exponent = zeros(nch, nch);
end

% Initialize helper variables
Xi = zeros(nch, ns);
sx = zeros(nch, ns);
sxi = zeros(nch, ns);
sx2 = zeros(nch, nch, ns);
ix = zeros(1, ns);

% Initialize storage for F^2_DCCA values
F2_DCCA_storage = cell(1, ns);
for i = 1:ns
    F2_DCCA_storage{i} = zeros(nch, nch, floor(win_size/scales(i)));
end

window_idx = 1;
for t = 1:N
    % Skip artifacts by replacing with zeros
    if t <= length(artifact_mask) && artifact_mask(t)
        data_t = zeros(nch, 1);
    else
        data_t = data(:, t);
    end
    
    % Handle NaN values
    data_t(isnan(data_t)) = 0;
    
    % Update cumulative sum 
    Xi = Xi + repmat(data_t, [1,ns]);
    
    % Update helper variables
    ix = ix + 1;
    sx = sx + Xi;
    sxi = sxi + Xi .* ix;
    for s = 1:ns
        sx2(:,:,s) = sx2(:,:,s) + Xi(:,s) * Xi(:,s)';
    end
    
    % Process each scale
    for s = 1:ns
        if mod(t, scales(s)) == 0
            % Calculate trend coefficients
            m = -6 * (-2*sxi(:,s) + (scales(s)+1).*sx(:,s)) ...
                ./ (scales(s) * (scales(s)^2 - 1));
            b = sx(:,s) ./ scales(s) - m .* (scales(s)+1) ./ 2;
            
            % Calculate f^2_DCCA for current window
            f2_DCCA = calculate_f2_DCCA(sx(:,s), sxi(:,s), sx2(:,:,s), m, b, scales(s), channel_mask);
            
            % Store in circular buffer
            buffer_idx = mod(t/scales(s)-1, floor(win_size/scales(s)))+1;
            F2_DCCA_storage{s}(:,:,buffer_idx) = f2_DCCA;
            
            % Reset helper variables for this scale
            Xi(:,s) = zeros(nch, 1);
            sx(:,s) = zeros(nch, 1);
            sxi(:,s) = zeros(nch, 1);
            sx2(:,:,s) = zeros(nch, nch);
            ix(s) = 0;
        end
    end
    
    % Compute metrics when we have enough data and at each step_size
    if t >= win_size && mod(t - win_size, step_size) == 0
        % Average F2_DCCA for each scale
        F2_DCCA = zeros(nch, nch, ns);
        for s = 1:ns
            F2_DCCA(:,:,s) = mean(F2_DCCA_storage{s}, 3);
        end
        
        % Calculate DCCC
        cur_DCCC = calculate_DCCC(F2_DCCA, nch, ns, channel_mask);
        
        % Calculate scaling exponent
        alpha = calculate_scaling_exponent(F2_DCCA, scales, nch, channel_mask);
        
        % Store results
        if strcmpi(output_mode, 'full')
            DCCA(:,:,:,window_idx) = F2_DCCA;
            DCCC(:,:,:,window_idx) = cur_DCCC;
            scaling_exponent(:,:,window_idx) = alpha;
        else % 'final' mode - overwrite previous results
            DCCA = F2_DCCA;
            DCCC = cur_DCCC;
            scaling_exponent = alpha;
        end
        
        % Update window index
        window_idx = window_idx + 1;
        
        % Adjust for the next window by removing oldest points
        if t + 1 <= N && step_size > 0
            % Remove step_size points of data from cumulative sum
            for s = 1:ns
                old_data = data(:, (t-win_size+1):(t-win_size+step_size));
                old_data(isnan(old_data)) = 0;
                Xi(:,s) = Xi(:,s) - sum(old_data, 2);
            end
        end
    end
end

% Helper functions
function filtered_data = bandpass_filter(data, band, fs)
    % Define frequency bands
    bands = struct('delta', [1 4], 'theta', [4 8], 'alpha', [8 13], ...
                   'beta', [13 30], 'gamma', [30 45]);
    
    % Skip filtering for broadband
    if strcmpi(band, 'broadband')
        filtered_data = data;
        return;
    end
    
    % Apply bandpass filter
    f_band = bands.(lower(band));
    [b, a] = butter(4, f_band/(fs/2), 'bandpass');
    filtered_data = filtfilt(b, a, double(data'))';
end

function f2_DCCA = calculate_f2_DCCA(sx, sxi, sx2, m, b, s, mask)
    % Only calculate for needed channel pairs to save computation
    % Calculate detrended covariance according to Equation 16 in the paper
    term1 = (1/s) * sx2;
    term2 = -(1/s) * (m * sxi' + sxi * m');
    term3 = -(1/s) * (b * sx' + sx * b');
    term4 = ((s+1)*(2*s+1)/6) * (m * m');
    term5 = ((s+1)/2) * (m * b' + b * m');
    term6 = b * b';
    
    % Combine all terms only for selected channel pairs
    f2_DCCA = zeros(size(term1));
    
    % Apply mask to only compute necessary pairs
    pairs = find(mask);
    for idx = 1:length(pairs)
        [i, j] = ind2sub(size(mask), pairs(idx));
        f2_DCCA(i,j) = term1(i,j) + term2(i,j) + term3(i,j) + term4(i,j) + term5(i,j) + term6(i,j);
    end
    
    % Ensure numerical stability
    f2_DCCA = max(f2_DCCA, eps);  % Ensure positive values
end

function DCCC = calculate_DCCC(F2_DCCA, nch, ns, mask)
    DCCC = zeros(nch, nch, ns);
    for i = 1:ns
        % Extract DFA values from diagonal
        diag_DCCA = diag(F2_DCCA(:,:,i));
        
        % Create denominator matrix for normalization
        DFA_matrix = sqrt(diag_DCCA) * sqrt(diag_DCCA)';
        
        % Calculate DCCC only for needed pairs
        temp_DCCC = zeros(nch, nch);
        pairs = find(mask);
        for idx = 1:length(pairs)
            [ch1, ch2] = ind2sub([nch, nch], pairs(idx));
            
            % Skip if denominator is zero
            if DFA_matrix(ch1, ch2) <= eps
                temp_DCCC(ch1, ch2) = 0;
                continue;
            end
            
            value = F2_DCCA(ch1, ch2, i) / DFA_matrix(ch1, ch2);
            
            % Ensure within -1 to 1 range
            temp_DCCC(ch1, ch2) = min(max(value, -1), 1);
        end
        
        % Set diagonal to 1 (autocorrelation)
        temp_DCCC(logical(eye(nch))) = 1;
        
        DCCC(:,:,i) = temp_DCCC;
    end
end

function alpha = calculate_scaling_exponent(F2_DCCA, scales, nch, mask)
    alpha = zeros(nch, nch);
    log_s = log(scales(:));
    
    % Only calculate for needed channel pairs
    pairs = find(mask);
    for idx = 1:length(pairs)
        [i, j] = ind2sub([nch, nch], pairs(idx));
        
        % Extract F2 values for this channel pair across scales
        log_F2 = log(squeeze(F2_DCCA(i,j,:)));
        
        % Check for invalid values
        if any(~isfinite(log_F2)) || any(log_F2 <= 0)
            alpha(i,j) = NaN;
            continue;
        end
        
        % Linear regression to find scaling exponent
        X = [ones(length(log_s), 1), log_s];
        b = X \ log_F2;
        alpha(i,j) = b(2); % Slope is the scaling exponent
    end
end
end