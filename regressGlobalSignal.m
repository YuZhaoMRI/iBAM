function residual=regressGlobalSignal(brain_data,grayMatterMaskFile)
%REGRESSGLOBALSIGNAL Removes the global signal from fMRI data via regression.
%   residual = REGRESSGLOBALSIGNAL(data4D, brainMask) performs global
%   signal regression (GSR) on the input 2D fMRI data [time x voxels].
%   residual = REMOVE_GSIGNAL(brain_data) performs global signal regression
%   (GSR) on the input brain data. The global signal is calculated as the
%   average time series across all voxels (ignoring NaNs) in gray matter. It is then
%   detrended (mean removed) and regressed out from each voxel's time series.
%
%   Input:
%       brain_data  - 2D matrix of fMRI BOLD signals [time x voxels].                   
%       grayMatterMaskFile - gray matter mask.nii
%
%   Output:
%       residual    - Data after global signal removal [time x voxels].
%
%   Note:
%       - NaN values are ignored when computing the global signal.
%       - A constant intercept term is included in the regression model.
%       - Global signal regression can reduce physiological noise but is
%         controversial as it may introduce artifactual negative correlations.
%
%   See also: detrend, regressGlobalSignal.

    % Step 1: Compute the global signal (mean across voxels at each time point)
    masNii=load_nii(grayMatterMaskFile);
    GMMask3D=masNii.img;
    GMMask1D=reshape(GMMask3D,[prod(size(GMMask3D)) 1]);
    GMIndex=find(GMMask1D==1);
    global_signal = nanmean(brain_data(:,GMIndex), 2);  % [T x 1]
    
    % Optional (recommended): Demean the global signal to avoid constant bias
    global_signal = detrend(global_signal, 0);  % Remove mean (constant term)
    % global_signal = detrend(global_signal);   % Further remove linear trend (uncomment if needed)
    
    % Step 2: Regress out the global signal from each voxel
    % Linear model: data(t) = beta * global_signal(t) + intercept + residual
    % We want the residual (data with global signal removed)
    
    % Design matrix: include global signal and constant intercept
    X_reg = [global_signal, ones(size(global_signal, 1), 1)];  % [T x 2]
    
    % Perform linear regression for all voxels simultaneously
    beta = X_reg \ brain_data;         % Regression coefficients [2 x voxels]
    
    % Predicted component from global signal and intercept
    predicted = X_reg * beta;          % [T x voxels]
    
    % Residuals: data after removing global signal influence
    residual = brain_data - predicted; % [T x voxels]

end