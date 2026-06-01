function data2D = applyBandpassFilter(data2D, TR, lowFreq, highFreq)
%APPLYBANDPASSFILTER Applies temporal band-pass filtering to 4D fMRI data.
%
%   dataFiltered = APPLYBANDPASSFILTER(data2D, TR, lowFreq, highFreq)
%   performs band-pass filtering on 2D fMRI data [time x voxels] to remove low-
%   frequency drift and high-frequency physiological noise while preserving
%   the phase of the BOLD signal.
%
%   The method combines:
%     - High-pass filtering using DCT basis (as in SPM's spm_filter) to
%       remove slow drifts.
%     - Zero-phase low-pass filtering using a 4th-order Butterworth filter
%       (via filtfilt) to avoid phase distortion.
%
%   Inputs:
%       data2D    - 2D fMRI data [time x voxels]
%       TR        - Repetition time in seconds (e.g., 2.0).
%       lowFreq   - Lower cutoff frequency in Hz (e.g., 0.01). Default: 0.01 Hz.
%       highFreq  - Upper cutoff frequency in Hz (e.g., 0.18). Default: 0.18 Hz.
%
%   Output:
%       dataFiltered - Band-pass filtered 4D data (same size as input).
%
%   Notes:
%       - Typical resting-state band: 0.01–0.1 Hz or 0.01–0.18 Hz.
%       - Zero-phase filtering (filtfilt) ensures no phase delay in BOLD signal.
%       - Requires SPM toolbox for spm_filter (high-pass DCT).
%
%   Example:
%       dataFiltered = applyBandpassFilter(vol4D, 0.8, 0.01, 0.18);
%
%   References:
%       - Common in CONN toolbox, DPARSF, and SPM-based pipelines.
%       - Power et al. (2014). Methods to detect, characterize, and remove motion artifact.
%
%   See also: spm_filter, butter, filtfilt, detrend.

    % Default frequency cutoffs (common in resting-state fMRI)
    if nargin < 3 || isempty(lowFreq),  lowFreq = 0.01; end   % Hz
    if nargin < 4 || isempty(highFreq), highFreq = 0.18; end  % Hz
    

    
    %% High-pass filtering using DCT (removes low-frequency drift)
    HParam = 1 / lowFreq;  % High-pass cutoff period in seconds
    K = struct('RT', TR, 'HParam', HParam, 'row', 1:size(data2D,1));
    K = spm_filter(K);     % Generate DCT high-pass filter matrix
    data2D = spm_filter(K, data2D);  % Apply high-pass
    
    %% Low-pass filtering using zero-phase Butterworth filter
    fs = 1 / TR;                   % Sampling frequency (Hz)
    nyquist = fs / 2;              % Nyquist frequency
    Wn = highFreq / nyquist;       % Normalized cutoff frequency
    [b, a] = butter(4, Wn, 'low'); % 4th-order low-pass Butterworth
    data2D = filtfilt(b, a, data2D);  % Zero-phase forward-backward filtering
    
end