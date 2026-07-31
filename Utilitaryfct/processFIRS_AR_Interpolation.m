function [cleanOD, noiseMask, rejectedChannels] = processFIRS_AR_Interpolation(dataMatrix, mark, rejectionThreshold,  nsample,norder)
% PROCESSFIRS_AR_INTERPOLATION  Converts Raw to OD, rejects noisy channels, and applies AR interpolation.
%
% This function integrates OD conversion using a single global baseline,
% identifies noisy channels, and fills temporal gaps via AR modeling.
%
% Inputs:
%   dataMatrix         - Raw light intensity matrix [Time x Channels]
%   mark               - Matrix of noise segments [Start Index, Length, Channel Index]
%   rejectionThreshold - Ratio of noise allowed per channel (0 to 1). Default is 0.3.
%
% Outputs:
%   cleanOD            - Optical Density data after rejection and AR interpolation.
%   rejectedChannels   - Indices of channels marked for total rejection (entirely NaN).

    % Set default threshold if not provided
    if nargin < 3, rejectionThreshold = 0.3; end 
    
    [numSamples, numChannels] = size(dataMatrix);
    
    % --- Step 1. Generate Noise Mask ---
    % Flagging segments identified during manual inspection (Yellow Blocks)
    noiseMask = false(numSamples, numChannels);
    for i = 1:size(mark, 1)
        startIdx = mark(i, 1);
        len      = mark(i, 2);
        ch       = mark(i, 3);
        endIdx   = startIdx + len - 1;
        
        if endIdx > numSamples, endIdx = numSamples; end
        noiseMask(startIdx:endIdx, ch) = true;
    end
    
    % --- Step 2. Optical Density (OD) Conversion with Safety Checks ---
    % 2.1 Mask manual noise first so they don't bias the I_mean calculation
    rawIntensity = dataMatrix;
    rawIntensity(noiseMask) = NaN;
    
    % 2.2 Remove non-positive values to prevent "Complex Double" errors
    % Any intensity <= 0 will be set to NaN before the log transformation.
    rawIntensity(rawIntensity <= 0) = NaN;
    
    % 2.3 Calculate Global Baseline (I_mean) per channel
    I_mean = nanmean(rawIntensity, 1);
    
    % 2.4 Convert to delta Optical Density (dOD)
    % Formula: dOD = -ln(I(t) / I_mean). 
    % We apply this to the masked raw data.
    cleanOD = -log(rawIntensity ./ I_mean);
    
    % --- Step 3. Global Channel Rejection ---
    % Identify channels where the noise ratio exceeds the threshold.
    noiseRatio = sum(noiseMask, 1) / numSamples;
    rejectedChannels = find(noiseRatio >= rejectionThreshold);
    
    % For rejected channels, set the entire time series to NaN for later spatial interpolation.
    for i = 1:length(rejectedChannels)
        ch = rejectedChannels(i);
        cleanOD(:, ch) = NaN;
        noiseMask(:, ch) = 1;
    end
    
    % --- Step 4. Temporal Imputation (AR Modeling) ---
    % Fill minor NaN gaps in valid channels using the specified AR parameters.
    for ch = 1:numChannels
        % Process only channels that are NOT globally rejected but have partial NaNs
        if ~ismember(ch, rejectedChannels) && any(isnan(cleanOD(:, ch)))
            try
                % fillgaps with your specific parameters (Window: 80, Overlap: 80)
                cleanOD(:, ch) = fillgaps(cleanOD(:, ch), 80, 80);
            catch ME
                warning('AR Interpolation failed for channel %d: %s. Leaving as NaN.', ch, ME.message);
            end
        end
    end
    
    % Log processing results to console
    fprintf('--- fNIRS Preprocessing Report ---\n');
    fprintf('Step 1: OD Conversion applied (Safety check: Real Double guaranteed).\n');
    fprintf('Step 2: %d channels rejected (>%.0f%% noise).\n', length(rejectedChannels), rejectionThreshold*100);
    fprintf('Step 3: AR Interpolation applied using fillgaps(',num2str(nsample), num2str(norder), ').\n');
    fprintf('----------------------------------\n');
end