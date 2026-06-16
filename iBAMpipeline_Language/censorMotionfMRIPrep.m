function [censorMask, FD, badFramesIdx] = censorMotionfMRIPrep(confoundsFile, fdThreshold, backFrames, forwardFrames, isPlotting, maxSpikeRatio)
% CENSORMOTIONFMRIPREP Generates a censoring mask from fMRIPrep FD confounds.
%
%   SYNTAX:
%       censorMask = censorMotionfMRIPrep(confoundsFile)
%       censorMask = censorMotionfMRIPrep(confoundsFile, fdThreshold, backFrames, forwardFrames)
%       [censorMask, FD, badFramesIdx] = censorMotionfMRIPrep(confoundsFile, 0.25, 1, 2, true, 0.025)
%
%   DESCRIPTION:
%       This function reads the fMRIPrep desc-confounds_timeseries.tsv file
%       and uses the existing framewise_displacement column to identify
%       high-motion frames. It applies the same adaptive censoring strategy
%       used by censorMotionHCP:
%
%       1. Mark frames with FD above an absolute threshold.
%       2. If the raw spike count exceeds maxSpikeRatio of the run, switch
%          to a relative strategy and isolate the top maxSpikeRatio FD frames.
%       3. Dilate each selected spike by censoring neighboring frames.
%
%   INPUT ARGUMENTS:
%       confoundsFile  - String. Path to fMRIPrep desc-confounds_timeseries.tsv.
%       fdThreshold    - Double. FD threshold in millimeters (Default: 0.25).
%       backFrames     - Integer. Number of volumes to censor before a spike (Default: 1).
%       forwardFrames  - Integer. Number of volumes to censor after a spike (Default: 2).
%       isPlotting     - Logical. Enable/disable diagnostic figure (Default: false).
%       maxSpikeRatio  - Double. Proportional guardrail and relative fallback
%                        fraction (Default: 0.025 [2.5%]).
%
%   OUTPUT ARGUMENTS:
%       censorMask     - Logical column vector [NumFrames x 1].
%                        True indicates clean volumes to retain.
%                        False indicates volumes to scrub/discard.
%       FD             - Framewise displacement column used for censoring.
%       badFramesIdx   - Primary high-motion frames before neighbor dilation.
%
%   Notes:
%       fMRIPrep's first framewise_displacement value is often n/a/NaN
%       because there is no previous volume. This function sets NaN FD values
%       to 0 before thresholding.
%
%   Author: Yu Zhao; adapted for fMRIPrep confounds by Codex
%   Date: June 2026
%   Version: 1.0

    %% 1. Handle Input Arguments and Defaults
    if nargin < 1
        error('CRITICAL ERROR: The path to the fMRIPrep confounds file must be provided.');
    end
    if nargin < 2 || isempty(fdThreshold),   fdThreshold = 0.25;  end
    if nargin < 3 || isempty(backFrames),    backFrames = 1;      end
    if nargin < 4 || isempty(forwardFrames), forwardFrames = 2;   end
    if nargin < 5 || isempty(isPlotting),    isPlotting = false;  end
    if nargin < 6 || isempty(maxSpikeRatio), maxSpikeRatio = 0.025; end

    if fdThreshold < 0
        error('fdThreshold must be non-negative.');
    end
    if backFrames < 0 || forwardFrames < 0
        error('backFrames and forwardFrames must be non-negative.');
    end
    if maxSpikeRatio <= 0 || maxSpikeRatio > 1
        error('maxSpikeRatio must be in the interval (0, 1].');
    end

    %% 2. Load fMRIPrep Confounds and Read Existing FD
    if ~exist(confoundsFile, 'file')
        error('CRITICAL ERROR: File not found at path: %s', confoundsFile);
    end

    opts = detectImportOptions(confoundsFile, 'FileType', 'text', 'Delimiter', '\t');
    if isprop(opts, 'VariableNamingRule')
        opts.VariableNamingRule = 'preserve';
    end
    confoundsTable = readtable(confoundsFile, opts);

    fdColumnName = 'framewise_displacement';
    fdColumnIdx = find(strcmp(confoundsTable.Properties.VariableNames, fdColumnName), 1);
    if isempty(fdColumnIdx)
        error('CRITICAL ERROR: Column "%s" was not found in: %s', fdColumnName, confoundsFile);
    end

    rawFD = confoundsTable{:, fdColumnIdx};
    if isnumeric(rawFD)
        FD = rawFD;
    else
        FD = str2double(string(rawFD));
    end
    FD = FD(:);
    FD(isnan(FD)) = 0;

    numTotalFrames = length(FD);
    if numTotalFrames == 0
        error('CRITICAL ERROR: No FD values were found in: %s', confoundsFile);
    end

    %% 3. Apply Adjustable Proportional Guardrail Logic
    badFramesIdx = find(FD > fdThreshold);
    isFallbackTriggered = false;
    effectiveThresholdStr = sprintf('%.2f mm (Absolute)', fdThreshold);

    maxAllowedSpikes = numTotalFrames * maxSpikeRatio;

    if length(badFramesIdx) > maxAllowedSpikes
        isFallbackTriggered = true;

        numTopSpikes = max(1, round(numTotalFrames * maxSpikeRatio));
        [~, sortedIdx] = sort(FD, 'descend');
        badFramesIdx = sortedIdx(1:numTopSpikes);
        badFramesIdx = sort(badFramesIdx, 'ascend');

        cutoffValue = min(FD(badFramesIdx));
        effectiveThresholdStr = sprintf('Top %.1f%% Relative (Cutoff: %.4f mm)', ...
            maxSpikeRatio * 100, cutoffValue);
    end

    %% 4. Temporal Window Dilated Censoring
    censorMask = true(numTotalFrames, 1);

    for i = 1:length(badFramesIdx)
        currFrame = badFramesIdx(i);

        startFrame = max(1, currFrame - backFrames);
        endFrame   = min(numTotalFrames, currFrame + forwardFrames);

        censorMask(startFrame:endFrame) = false;
    end

    remainingFrames = sum(censorMask);
    retentionRate   = (remainingFrames / numTotalFrames) * 100;

    fprintf('\n================== MOTION SCRUBBING METRICS ==================\n');
    fprintf('Source File: %s\n', confoundsFile);
    fprintf('FD Source: fMRIPrep framewise_displacement column\n');
    fprintf('Total Volumes: %d | Retained: %d | Scrubbed: %d\n', ...
        numTotalFrames, remainingFrames, numTotalFrames - remainingFrames);
    if isFallbackTriggered
        fprintf('NOTICE: Raw spikes exceeded %.1f%% of run. Relative mode activated.\n', maxSpikeRatio * 100);
    end
    fprintf('Effective Strategy: %s\n', effectiveThresholdStr);
    fprintf('Primary Spikes Isolated: %d\n', length(badFramesIdx));
    fprintf('Final Clean Data Retention Rate: %.2f%%\n', retentionRate);
    fprintf('==============================================================\n');

    %% 5. Diagnostics Visualization
    if isPlotting
        yMax = max([FD; fdThreshold]) * 1.1;
        if yMax <= 0
            yMax = 1;
        end

        figure('Color', 'w', 'Position', [150, 150, 1050, 420]);
        hold on;

        for t = 1:numTotalFrames
            if ~censorMask(t)
                fill([t-0.5, t+0.5, t+0.5, t-0.5], [0, 0, yMax, yMax], ...
                    [1 0.88 0.88], 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end

        plot(FD, 'LineWidth', 1.5, 'Color', [0.15 0.15 0.15], ...
            'DisplayName', 'fMRIPrep Framewise Displacement (FD)');

        if ~isFallbackTriggered
            yline(fdThreshold, '--r', 'LineWidth', 1.5, ...
                  'Label', ['Absolute Threshold (' num2str(fdThreshold) ' mm)'], ...
                  'LabelHorizontalAlignment', 'right', 'FontSize', 10, 'FontWeight', 'bold');
        else
            if ~isempty(badFramesIdx)
                cutoffValue = min(FD(badFramesIdx));
                yline(cutoffValue, '--b', 'LineWidth', 1.5, ...
                      'Label', [sprintf('Adaptive Top %.1f%% Cutoff (', maxSpikeRatio * 100) num2str(cutoffValue, '%.3f') ' mm)'], ...
                      'LabelHorizontalAlignment', 'right', 'FontSize', 10, 'FontWeight', 'bold');
            end
        end

        grid on;
        box on;
        set(gca, 'Layer', 'top', 'FontSize', 11, 'LineWidth', 1.1, 'TickDir', 'out');
        xlabel('Timepoints (Frames)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('FD (mm)', 'FontSize', 12, 'FontWeight', 'bold');

        titleStr = sprintf('fMRIPrep Motion Scrubbing Profile (Data Retention Rate: %.2f%%)', retentionRate);
        if isFallbackTriggered
            titleStr = [titleStr, sprintf(' [Adaptive %.1f%% Mode Active]', maxSpikeRatio * 100)];
        end
        title(titleStr, 'FontSize', 13, 'FontWeight', 'bold');

        xlim([1, numTotalFrames]);
        ylim([0, yMax]);

        plot(nan, nan, 'Color', [1 0.88 0.88], 'LineWidth', 8, ...
            'DisplayName', 'Censored Volumes (Neighbors Included)');
        legend('Location', 'northeast', 'Box', 'off');

        hold off;
    end
end
