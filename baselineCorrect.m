function [Yc, baseLineInfo] = baselineCorrect(Y, opts)
%ROBUSTBASELINECORRECT Performs robust baseline correction and noise normalization on fluorescence traces.
%
%   [Yc, normalization_info] = ROBUSTBASELINECORRECT(Y, opts) applies a robust
%   baseline estimation. This method is particularly
%   useful for noisy traces with sparse transients, as it avoids contamination
%   of the baseline by activity peaks.
%
%   The procedure:
%     1. Estimates per-trace noise level (sigma) using median absolute deviation
%        (MAD) on first differences (robust to outliers/transients).
%     2. Estimates baseline as the expected value of a low quantile (p_baseline)
%        under Gaussian noise assumption.
%     3. Subtracts the baseline and divides by sigma, yielding approximately
%        zero-mean, unit-variance signals (robust dF/sigma or z-scored traces).
%
%   Inputs:
%       Y         - fMRI matrix [time × traces] (T × N). 
%       opts      - (optional) Struct with options:
%                   .p_baseline  - Quantile fraction for baseline estimation
%                                  (default: 0.10, i.e., 10th percentile).
%                   .topk_frac   - Currently unused (reserved for future extensions,
%                                  e.g., excluding top fraction of values).
%
%   Outputs:
%       Yc               - Baseline-corrected and noise-normalized traces [T × N].
%       baseLineInfo - Struct containing:
%                            .baselines - Estimated baseline per trace [1 × N].
%                            .sigma     - Estimated noise std per trace [1 × N].
%
%   Notes:
%       - Noise estimation uses MAD on temporal differences, scaled for Gaussian
%         consistency (factor 1/(0.6745*sqrt(2))).
%       - Baseline uses the inverse error function to convert quantile to mean
%         under symmetric noise.
%       - Handles constant traces gracefully (avoids division by zero).
%       - Suitable for calcium imaging preprocessing before spike inference or
%         activity detection.
%
%
%   References:
%       - Common in calcium imaging pipelines (e.g., Suite2p, CaImAn).
%       - MAD-based noise estimation: standard robust statistic.
%


    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    
    p_baseline = getOpt(opts, 'p_baseline', 0.10);
    topk_frac  = getOpt(opts, 'topk_frac', 0.01);  % Reserved/unused in current version
    
    [T, N] = size(Y);
    
    %% 1) Robust noise estimation (sigma) using difference-MAD
    dY = diff(Y, 1, 1);                  % (T-1) × N temporal differences
    med_dY = median(dY, 1);              % Median per trace
    mad_diff = median(abs(dY - med_dY), 1);
    sigma = mad_diff / (0.6745 * sqrt(2));  % Scale to Gaussian std equivalent
    
    % Fallback for zero/constant sigma (e.g., dead pixels or flat traces)
    zero_mask = sigma == 0;
    if any(zero_mask)
        medY = median(Y, 1);
        sigma(zero_mask) = median(abs(Y(:, zero_mask) - medY(zero_mask)), 1) / 0.6745;
    end
    
    %% 2) Robust baseline estimation using low quantile
    Q_low = quantile(Y, p_baseline, 1);   % 1 × N
    z_low = sqrt(2) * erfinv(2 * p_baseline - 1);  % Theoretical z-score for quantile
    baselines = Q_low - sigma .* z_low;
    
    %% 3) Baseline-corrected signals
    Yc = Y - repmat(baselines, T, 1);
    
    %% 4) Normalize by estimated noise (unit variance)
    Yc = Yc ./ repmat(sigma, T, 1);

    
    % Package normalization parameters for later use
    baseLineInfo = struct('baselines', baselines, ...
                                'sigma', sigma);

end

% --- Helper function for options ---
function v = getOpt(s, name, def)
    if isfield(s, name)
        v = s.(name);
    else
        v = def;
    end
end


