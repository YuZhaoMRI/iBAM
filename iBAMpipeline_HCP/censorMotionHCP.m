function censorMask = censorMotionHCP(regressorFile, fdThreshold, backFrames, forwardFrames, isPlotting, maxSpikeRatio)
% CENSORMOTIONHCP Computes Framewise Displacement (FD) and generates an adaptive censoring mask.
%
%   SYNTAX:
%       censorMask = censorMotionHCP(regressorFile)
%       censorMask = censorMotionHCP(regressorFile, fdThreshold, backFrames, forwardFrames)
%       censorMask = censorMotionHCP(regressorFile, 0.25, 1, 2, true, 0.025)
%
%   DESCRIPTION:
%       This function calculates Framewise Displacement (FD) using the 12-column
%       HCP 'Movement_Regressors.txt' file based on Power's method (Power et al., 2012).
%       
%       ADAPTIVE CRITERIA:
%       If the number of raw frames exceeding the absolute threshold (default: 0.25 mm) 
%       surpasses a user-defined ratio (default: 2.5% of the total scan length), 
%       the algorithm dynamically pivots to a relative thresholding strategy. It isolates 
%       exactly the top X% highest FD timepoints as primary motion spikes ('badFramesIdx') 
%       to protect data retention from systemic micro-movements, followed by temporal 
%       window dilation.
%
%   INPUT ARGUMENTS:
%       regressorFile  - String. Path to the HCP 'Movement_Regressors.txt' file.
%       fdThreshold    - Double. FD threshold in millimeters (Default: 0.25).
%       backFrames     - Integer. Number of volumes to censor prior to a spike (Default: 1).
%       forwardFrames  - Integer. Number of volumes to censor following a spike (Default: 2).
%       isPlotting     - Logical. Boolean flag to enable (true) or disable (false) 
%                        the creation of the diagnostic motion figure (Default: true).
%       maxSpikeRatio  - Double. Proportional guardrail limit and relative fallback 
%                        percentile expressed as a decimal (Default: 0.025 [2.5%]).
%
%   OUTPUT ARGUMENTS:
%       censorMask     - Logical column vector [NumFrames x 1]. 
%                        True (1) indicates clean volumes to retain.
%                        False (0) indicates contaminated volumes to scrub/discard.
%
%   Author: Yu Zhao
%   Date: June 2026
%   Version: 1.6

    %% 1. Handle Input Arguments and Defaults
    if nargin < 1
        error('CRITICAL ERROR: The path to the regressor file must be provided.');
    end
    if nargin < 2 || isempty(fdThreshold),   fdThreshold = 0.25; end % Updated to 0.25 mm requested default
    if nargin < 3 || isempty(backFrames),    backFrames = 1;     end
    if nargin < 4 || isempty(forwardFrames), forwardFrames = 2;    end
    if nargin < 5 || isempty(isPlotting),    isPlotting = false;  end 
    if nargin < 6 || isempty(maxSpikeRatio), maxSpikeRatio = 0.025; end % Adjustable boundary (e.g., 0.025 = 2.5%)

    %% 2. Load Regressors and Calculate FD (Power's Method)
    if ~exist(regressorFile, 'file')
        error('CRITICAL ERROR: File not found at path: %s', regressorFile);
    end

    allRegressors = load(regressorFile);
    motionParams = allRegressors(:, 1:6);

    % Transform rotational displacements from angular degrees to linear millimeters 
    motionParams(:, 4:6) = motionParams(:, 4:6) * (2 * pi * 50 / 360);

    % Compute backward temporal differences
    motionDiff = diff(motionParams);
    motionDiff = [zeros(1, 6); motionDiff]; 

    % Calculate FD as the instantaneous absolute sum of all 6 parameters
    FD = sum(abs(motionDiff), 2);
    numTotalFrames = length(FD);

    %% 3. Apply Adjustable Proportional Guardrail Logic
    badFramesIdx = find(FD > fdThreshold);
    isFallbackTriggered = false;
    effectiveThresholdStr = sprintf('%.2f mm (Absolute)', fdThreshold);
    
    % Define the maximum allowable raw spike count based on the adjustable parameter
    maxAllowedSpikes = numTotalFrames * maxSpikeRatio;

    % Trigger condition: If absolute bad frames exceed the user-defined percentage boundary
    if length(badFramesIdx) > maxAllowedSpikes
        isFallbackTriggered = true;
        
        % Calculate exactly how many frames represent the top percentile
        numTopSpikes = max(1, round(numTotalFrames * maxSpikeRatio)); 
        
        % Sort FD in descending order to capture the worst tail entries
        [~, sortedIdx] = sort(FD, 'descend');
        badFramesIdx = sortedIdx(1:numTopSpikes);
        
        % Re-sort indices chronologically to preserve safe downstream temporal window mapping
        badFramesIdx = sort(badFramesIdx, 'ascend');
        
        % Define the effective operational threshold line based on the relative cutoff boundary
        effectiveThresholdStr = sprintf('Top %.1f%% Relative (Cutoff: %.4f mm)', ...
            maxSpikeRatio * 100, FD(badFramesIdx(end)));
    end

    %% 4. Temporal Window Dilated Censoring
    censorMask = true(numTotalFrames, 1);

    % Linearly dilate the mask around the final validated set of spikes
    for i = 1:length(badFramesIdx)
        currFrame = badFramesIdx(i);
        
        startFrame = max(1, currFrame - backFrames);
        endFrame   = min(numTotalFrames, currFrame + forwardFrames);
        
        censorMask(startFrame:endFrame) = false;
    end

    % Calculate final data exclusion metrics
    remainingFrames = sum(censorMask);
    retentionRate   = (remainingFrames / numTotalFrames) * 100; 

    fprintf('\n================== MOTION SCRUBBING METRICS ==================\n');
    fprintf('Source File: %s\n', regressorFile);
    fprintf('Total Volumes: %d | Retained: %d | Scrubbed: %d\n', ...
        numTotalFrames, remainingFrames, numTotalFrames - remainingFrames);
    if isFallbackTriggered
        fprintf('NOTICE: Raw spikes exceeded %.1f%% of run. Relative mode activated.\n', maxSpikeRatio * 100);
    end
    fprintf('Effective Strategy: %s\n', effectiveThresholdStr);
    fprintf('Primary Spikes Isolated: %d\n', length(badFramesIdx));
    fprintf('Final Clean Data Retention Rate: %.2f%%\n', retentionRate);
    fprintf('==============================================================\n');

    %% 5. Diagnostics Visualization (Conditional Execution)
    if isPlotting
        figure('Color', 'w', 'Position', [150, 150, 1050, 420]);
        hold on;

        % Shade censored timepoints with light crimson background columns
        for t = 1:numTotalFrames
            if ~censorMask(t)
                fill([t-0.5, t+0.5, t+0.5, t-0.5], [0, 0, max(FD)*1.1, max(FD)*1.1], ...
                    [1 0.88 0.88], 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end

        % Plot continuous Framewise Displacement time-series
        plot(FD, 'LineWidth', 1.5, 'Color', [0.15 0.15 0.15], 'DisplayName', 'Framewise Displacement (FD)');

        % Draw dynamic boundary lines depending on the activated mode
        if ~isFallbackTriggered
            yline(fdThreshold, '--r', 'LineWidth', 1.5, ...
                  'Label', ['Absolute Threshold (' num2str(fdThreshold) ' mm)'], ...
                  'LabelHorizontalAlignment', 'right', 'FontSize', 10, 'FontWeight', 'bold');
        else
            if ~isempty(badFramesIdx)
                cutoffValue = FD(badFramesIdx(end));
                yline(cutoffValue, '--b', 'LineWidth', 1.5, ...
                      'Label', [sprintf('Adaptive Top %.1f%% Cutoff (', maxSpikeRatio * 100) num2str(cutoffValue, '%.3f') ' mm)'], ...
                      'LabelHorizontalAlignment', 'right', 'FontSize', 10, 'FontWeight', 'bold');
            end
        end

        % Axis styling and layout properties
        grid on;
        box on;
        set(gca, 'Layer', 'top', 'FontSize', 11, 'LineWidth', 1.1, 'TickDir', 'out');
        xlabel('Timepoints (Frames)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('FD (mm)', 'FontSize', 12, 'FontWeight', 'bold');

        titleStr = sprintf('Motion Scrubbing Profile (Data Retention Rate: %.2f%%)', retentionRate);
        if isFallbackTriggered, titleStr = [titleStr, sprintf(' [Adaptive %.1f%% Mode Active]', maxSpikeRatio * 100)]; end
        title(titleStr, 'FontSize', 13, 'FontWeight', 'bold');

        xlim([1, numTotalFrames]);
        ylim([0, max(FD)*1.1]);

        % Generate dummy entries to format a clean graphical Legend
        plot(nan, nan, 'Color', [1 0.88 0.88], 'LineWidth', 8, 'DisplayName', 'Censored Volumes (Neighbors Included)');
        legend('Location', 'northeast', 'Box', 'off');

        hold off;
    end
end