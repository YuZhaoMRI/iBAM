
%   performStatisticalInference Non-parametric statistical inference engine for iBAM task activations.
%
%
%   DESCRIPTION:
%       This function serves as the standalone statistical inference engine for 
%       the iBAM (Implicit Brain Activity Mapping) framework. Instead of using 
%       direct functional inputs and outputs, it internally loads data matrices 
%       from hardcoded configuration paths and exports the final statistical maps 
%       directly to the disk.
%
%       It evaluates empirical Mode Intensity maps derived from task-fMRI against 
%       a non-parametric empirical null distribution constructed from null model Mode Intensity maps.
%
%       STATISTICAL INFERENCE & CORRECTION ROUTINE:
%       For each spatial node or voxel, a localized empirical cumulative distribution 
%       function (eCDF) is compiled from the null ensemble. The raw, non-parametric 
%       p-value is calculated based on the tail probability of the empirical task value. 
%       To strictly control the false positive rate across thousands of parallel 
%       statistical tests, a False Discovery Rate (FDR) correction is applied globally. 
%       The corrected probabilities are transformed into standard normal equivalent 
%       Z-scores, and all finalized metrics (p-values, Z-scores, and the binary 
%       activation mask) are saved automatically to the designated output directory.
%
%   CONFIGURATION (Specified within the function body):
%       - Task MI Map Path: Path to the '.mat' file containing the true empirical 
%                           task fMRI Mode Intensity maps [NumNodes x 1].
%       - Null MI Dir:      Path to the directory containing the compiled null model 
%                           Mode Intensity maps generated via stochastic resampling.
%       - Output Directory: Target folder where the computed pValues, zScores, and 
%                           activationMask arrays are saved as a unified workspace.
%       - Alpha Threshold:  Desired FDR significance level for the final mask (e.g., 0.05).
%
%
%   References:
%       Zhao, Y. et al. (2026). Implicit brain activity mapping reveals language 
%       as a higher-order strategic tool in human reasoning.
%
%   Author: Yu Zhao
%   Date: June 2026


clear;clc;

%% Define the target False Discovery Rate (FDR) control threshold (Adjustable parameter)
fdrAlpha = 0.05; % E.g., 0.05 (standard), 0.01 (strict), or 0.001 (highly conservative)


%% set up inputpath
inputpath='F:\Language_task_iBAM_results';
% Verify the directory's physical existence on disk
if exist(inputpath, 'dir')
    fprintf('>>> Success: Target input path validated and accessible:\n    %s\n', inputpath);
else
    % Throw a fatal error to halt the pipeline execution if the path is missing
    error('CRITICAL ERROR: Input directory not found at: \n    %s\nAborting pipeline execution to prevent downstream matrix read failures.', inputpath);
end
% outputpath
outputpath=inputpath;


%% load brain mask
maskFileName='../Masks/MNI152_T1_2mm_brain.nii';
[brainMask3D,brainMaskDim3D,brainIndex]=getBrainMaskInfo(maskFileName);


%% get the participant folders
% Get all contents in the specified path (both files and folders)
contents = dir(inputpath);
% Filter out folders (excluding '.' and '..')
sublist = contents([contents.isdir]);  % Get all folders
sublist = sublist(~ismember({sublist.name}, {'.', '..'}));  % Exclude current and parent directories


for m=1:length(sublist)
display(['Perform statistical inference for subject #',num2str(m)])
sub=char(sublist(m).name);
inputpath_sub=[inputpath,'\',sub];

%% load null model MI maps
niiFileName=[inputpath_sub,'\','sub',sub,'_SW_MI_NullModels.nii'];
niiData=load_nii(niiFileName);
vol4DNullModels=single(niiData.img);
[N1,N2,N3,numNullModels]=size(vol4DNullModels);
vol4DNullModels=reshape(vol4DNullModels,[N1*N2*N3 numNullModels]);
MI_NullModels=vol4DNullModels(brainIndex,:);
clear niiData

%% load task MI maps
niiFileName=[inputpath_sub,'\','sub',sub,'_SW_MI.nii'];
niiData=load_nii(niiFileName);
vol4DTask=single(niiData.img(:,:,:,1:end));
taskTimings=squeeze(vol4DTask(1,1,1,:));
[N1,N2,N3,numSlidingWindow]=size(vol4DTask);
vol4DTask=reshape(vol4DTask,[N1*N2*N3 numSlidingWindow]);
MI_task=vol4DTask(brainIndex,:);



% 1. Preallocate matrix to store empirical cumulative probabilities (Left-tail areas)
% Dimensions: [Number of Brain Nodes/Voxels  x  Number of Sliding Windows]
leftArea = zeros(length(brainIndex), numSlidingWindow); 

% 2. Execute parallelized voxel-wise empirical cumulative distribution mapping
parfor i = 1:length(brainIndex)

    % Extract the empirical null distribution of Mode Intensities for the current node/voxel
    currentNullDistribution = MI_NullModels(i, :);   

    % Generate a logical mask identifying statistical outliers within the null ensemble
    % Identify statistical outliers within the empirical null ensemble using a robust,
    % non-parametric Median Absolute Deviation (MAD) approach with a tuned ThresholdFactor of 2.
    nullOutlierMask = isoutlier(currentNullDistribution,'median', 'ThresholdFactor', 2);

    % Extract the linear indices of the valid, noise-structure-preserving null models
    validNullIndices = find(~nullOutlierMask);

    % Perform Kernel Density Estimation (KDE) to construct the empirical 
    % Cumulative Distribution Function (eCDF) using non-outlier null instances.
    % 'f' returns the cumulative probabilities, 'xi' returns the corresponding evaluation points.
    [f, xi] = ksdensity(MI_NullModels(i, validNullIndices), 'Function', 'cdf');
    
    % Map the true empirical task Mode Intensities (across all sliding windows) 
    % to the closest evaluation point 'xi' in the null distribution space.
    % This finds the point of minimum absolute distance for each window.
    [~, closestIdx] = min(abs(MI_task(i, :) - xi(:)), [], 1); 
    
    % Retrieve the corresponding empirical cumulative probability (left-tail area)
    leftArea(i, :) = f(closestIdx);
end

% Apply strict Winsorization/probability truncation thresholds.
% This prevents absolute mathematical zeros and ones, which would otherwise 
% yield infinite values (Inf / -Inf) during downstream standard normal Z-score conversion.
leftArea(leftArea < 0.00000001) = 0.00000001;
leftArea(leftArea > 0.99999999) = 0.99999999;




% Total number of voxels in the structural 3D bounding box
totalVoxelsNum = N1 * N2 * N3;

%  Preallocate 4D volumetric time-series matrices for Z-scores and p-values
zMapTimeSeries = zeros(totalVoxelsNum, numSlidingWindow);
pMapTimeSeries = zeros(totalVoxelsNum, numSlidingWindow);

% Map empirical cumulative probabilities to standard normal Z-scores (Right-tailed alternative)
% Under this mapping, leftArea close to 1 yields highly positive Z-scores (Task > Null).
zMapTimeSeries(brainIndex, :) = norminv(leftArea);

% 3. Compute corresponding raw p-values for the upper-tail hypothesis test
pMapTimeSeries(brainIndex, :) = 1 - leftArea;

% 4. Reshape the 2D flattened voxel arrays back into native 4D fMRI space
zMapTimeSeries = reshape(zMapTimeSeries, N1, N2, N3, numSlidingWindow);
pMapTimeSeries = reshape(pMapTimeSeries, N1, N2, N3, numSlidingWindow);
zMapTimeSeries_thresholded = zMapTimeSeries;
zMapTimeSeries_raw = zMapTimeSeries;

% 5. Execute False Discovery Rate (FDR) thresholding sequentially for each sliding window
for k = 1:numSlidingWindow
    
    % Extract and flatten the raw p-value map for the current time window
    currentPVector = reshape(pMapTimeSeries(:, :, :, k), [totalVoxelsNum, 1]);
    
    % Sort p-values in ascending order as strictly required by SPM's FDR engine
    sortedPVector = sort(currentPVector, 'ascend');
    
    % Calculate the dynamic critical Z-score threshold ('criticalZ') based on user-defined fdrAlpha
    % Arguments specified: 
    %   fdrAlpha   - Parameterized False Discovery Rate threshold (e.g., q <= 0.05).
    %   [1 Inf]    - Degrees of freedom (set to Inf for standard normal Z-distribution).
    %   'Z'        - Indicates statistical field type is Gaussian standard normal.
    %   1          - Number of contrasts being tested.
    %   sortedPVector - Ordered vector of empirical p-values within the search volume.
    %   1          - Topology indicator (voxel-wise FDR discovery control).
    criticalZ = spm_uc_FDR(fdrAlpha, [1 Inf], 'Z', 1, sortedPVector, 1);
    
    % Extract the 3D Z-score volume for the current time slice
    currentZMap3D = zMapTimeSeries(:, :, :, k);
    
    % Apply topological threshold masking: Set non-significant voxels (Z < criticalZ) to absolute zero
    currentZMap3D(currentZMap3D < criticalZ) = 0;
    
    % Write the thresholded, clean activation cluster matrix back into the 4D time-series array
    zMapTimeSeries_thresholded(:, :, :, k) = currentZMap3D;
end



%% label the corner of images with task timings;  
% -1 indicates neutual task-on periods; 
%  1 indicates fear task-on periods;and 0 indicates rest.
% for zMapTimeSeries_thresholded
taskLabels=find(taskTimings==1);
zMapTimeSeries_thresholded(:,1:2,1:2,taskLabels)=1;
zMapTimeSeries_thresholded(1:2,:,1:2,taskLabels)=1;
zMapTimeSeries_thresholded(1:2,1:2,:,taskLabels)=1;
% for zMapTimeSeries_raw
taskLabels=find(taskTimings==1);
zMapTimeSeries_raw(:,1:2,1:2,taskLabels)=1;
zMapTimeSeries_raw(1:2,:,1:2,taskLabels)=1;
zMapTimeSeries_raw(1:2,1:2,:,taskLabels)=1;



outputpath_sub=inputpath_sub;
%%  Export FDR-Thresholded Z-map time-series
niiData.img=niiData.img.*0;
niiData.img(:,:,:,1:numSlidingWindow)=zMapTimeSeries_thresholded;
niiZmapFileName_out = sprintf('%s\\sub%s_SW_Zmap_thresholded_FDR%.2f.nii', ...
    outputpath_sub, sub, fdrAlpha);
% Save the updated NIfTI architecture to the file system via the NIfTI toolbox
save_nii(niiData, niiZmapFileName_out);

%%  Export Raw Z-map time-series
niiData.img=niiData.img.*0;
niiData.img(:,:,:,1:numSlidingWindow)=zMapTimeSeries_raw;
niiZmapFileName_out=[outputpath_sub,'\','\sub',sub,'_SW_Zmap_raw.nii'];
save_nii(niiData,niiZmapFileName_out);

end


