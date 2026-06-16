
% GENERATENULLMODELS_SW Generate an empirical null distribution of Mode Intensities (MIs) for iBAM.
%
%   DESCRIPTION:
%       This script is specifically designed to construct an empirical null 
%       distribution of Mode Intensities (MIs). To determine whether an observed 
%       MI signifies genuine neural activation rather than stochastic fluctuations, 
%       the framework implements a formal empirical model selection procedure. 
%       Because MIs can emerge from either Gaussian or residual structured noises—
%       potentially leading to spurious detections—this pipeline tests the null 
%       hypothesis that the observed task-based MI is fully explained by the 
%       intrinsic noise components of the fMRI signal.
%
%       To estimate this empirical null distribution, the script constructs 
%       500 noise-structure-preserving null models by detecting Latent Signal 
%       Modes (LSMs) on a resting-state fMRI time series. Prior to null model 
%       construction, an ICA-based filtering strategy is employed to remove 
%       spontaneous resting-state neural activity and structured noise (NSSN-filtered), 
%       thereby preventing these components from artificially inflating the 
%       detected baseline model intensities. 
%
%       During execution, a 15-second temporal sliding window (consistent with 
%       the duration used in the task-fMRI analysis) is utilized to stochastically 
%       sample signals from the NSSN-filtered resting-state fMRI series. The 
%       Functional Topology-Preserving Modal Analysis (FT-PMA) is then implemented 
%       on these sampled signals to characterize their LSMs and MIs, ultimately 
%       yielding the empirical null distribution used for statistical inference.
%
%   INPUTS:
%       Path to rest fmri data - preprocessed rest fMRI time series in MNI 152 space.
%                                
%   OUTPUTS:
%       Path to MI maps of null models      - Mode Intensity (MI) maps.
%
%   REFERENCE:
%       For detailed methodology, particularly regarding the preprocessing pipeline 
%       and FSSW configuration, please refer to the Methods section of:
%       "Implicit brain activity mapping reveals language as a higher-order 
%       strategic tool in human reasoning."
%
%   Author: Yu Zhao
%   Date: May 2026
%   Version: 1.0



clear;clc;

%% get the participant folders
inputpath = 'E:\S1200_Retest\S1200\Emotion\';
% Get all contents in the specified path (both files and folders)
contents = dir(inputpath);
% Filter out folders (excluding '.' and '..')
sublist = contents([contents.isdir]);  % Get all folders
sublist = sublist(~ismember({sublist.name}, {'.', '..'}));  % Exclude current and parent directories

%% outputpath
outputpath='E:\HCP_Emotion_results_0605\Emotion\Test';
if ~isfolder(outputpath)
    mkdir(outputpath)
end

%% Set Parameters
%  Define the source directory where parameters were stored
parameterFile = fullfile(outputpath, 'pipeline_parameters.mat');
%  Verify the parameter file exists before attempting to load
if exist(parameterFile, 'file')
    % Load the variables directly into the current workspace
    load(parameterFile, 'freqBand', 'FWHM', 'NumDummyPoints', 'durHRF');   
    % RE-VALIDATION AND COUPLING CHECK:
    %   These parameters are loaded dynamically to enforce absolute methodological 
    %   consistency with the primary task-fMRI pipeline. 
    
    fprintf('Successfully reloaded pipeline parameters from: %s\n', parameterFile);
    fprintf('Current Configuration -> FWHM: %dmm, durHRF: %ds, Band: %.2f-%.2f Hz\n', ...
        FWHM, durHRF, freqBand(1), freqBand(2));
else
    % Throw a fatal error to halt execution if parameters are missing
    error('CRITICAL ERROR: Pipeline parameters file not found at %s. Aborting null model generation to prevent configuration mismatch.', parameterFile);
end
% Set other Parameters
numNullModels=500;       % Number of noise-structure-preserving null models to construct
durRest = 4 * 60;        % Duration of the resting-state fMRI scan (seconds)


maskFileName = '../Masks/MNI152_T1_2mm_brain.nii'; 
% CRITICAL WARNING:
% 1. The spatial resolution and matrix size of this mask MUST perfectly match the fMRI data.
% 2. DO NOT change or replace this mask file under normal circumstances, as it is strictly 
%    coupled/bound with the FSSW (Functional Structure-encoded Spatial Window).
% 3. If you absolutely must change this mask, you MUST recalculate and update the FSSW accordingly.

% Extract mask matrix, dimensions, and linear indices of brain voxels
[brainMask3D, brainMaskDim3D, brainIndex] = getBrainMaskInfo(maskFileName);

%% load the matrices that ecode  a functional structure-encoded spatial window (FSSW) predefined for each voxel.
load('../Masks/FSSW_K98.mat');

%% Initialize Parallel Computing Pool
% Delete any existing parallel pool to ensure a clean start
delete(gcp('nocreate')); 

% Open a new parallel pool with 30 workers for high-performance computing
% (Adjust the number of workers based on your computer's CPU hardware performance)
parpool(30);


%% generate null models for all subjects
display(['Number of subjects is ',num2str(length(sublist))]);
for m=1:length(sublist)

display(['generate null models for subject #',num2str(m)])
sub=char(sublist(m).name);
outputpath_sub=[outputpath,'\',sub];
if ~isfolder(outputpath_sub)
    mkdir(outputpath_sub)
end

% load neural activity and structured noise (NASN)-filtered rest fMRI volumes 
%% LR
%% Set flags to determine which experimental runs will be processed and analyzed
% load rest fMRI data
dataNii=load_nii(['E:\S1200_Retest\S1200\Rest_FIX_NS\',sub,'\rfMRI_REST1_LR\rfMRI_REST1_LR_NASN_filtered.nii.gz']);

%% We identify and exclude volumes with significant head motion contamination 
% (exceeding the specified FD threshold) along with their immediate temporal 
% neighbors. Additionally, the dummy scans are discarded 
% to eliminate T1-equilibration effects and ensure the remaining time-series 
% reflects a stable, steady-state fMRI signal.
motionLabels= censorMotionHCP(['E:\S1200_Retest\S1200\Rest_FIX_NS\',sub,'\rfMRI_REST1_LR\Movement_Regressors.txt']);
retainedFramesIdx = find(motionLabels);
retainedFramesIdx(retainedFramesIdx < 1+NumDummyPoints) = []; % remove the dummy scans 
fMRIData4D=double(dataNii.img(:,:,:,retainedFramesIdx));% remove the dummy and resting-state fMRI scans at the end 

[N1,N2,N3,N4]=size(fMRIData4D);
voxelsize=dataNii.original.hdr.dime.pixdim(2);
TR=dataNii.original.hdr.dime.pixdim(5);

%% spatial smoothinng 
for ii=1:N4
    fMRIData4D(:,:,:,ii)=imgaussfilt3(fMRIData4D(:,:,:,ii).*brainMask3D,FWHM/voxelsize/2.355);  %三维图像的三维高斯滤波
end

%% preprocessing of fMRI time courses, including (1) detrend (2) bandpass filter (3)normalize
fMRIData2D=reshape(fMRIData4D,[N1*N2*N3 N4]);
fMRIData2D=fMRIData2D(:,:)';
fMRIData2D=detrend(fMRIData2D,2);

%% performs global signal regression
grayMatterMaskFile='../Masks/GrayMatterMask.nii';
fMRIData2D=regressGlobalSignal(fMRIData2D,grayMatterMaskFile);


%% The combined temporal filtering strategy effectively removes noise without introducing phase delays into the BOLD signal.
fMRIData2D = applyBandpassFilter(fMRIData2D, TR,freqBand(1),freqBand(2));

%% normalize the fMRI signals
fMRIData2D_LR=zscoreTimeSeries(fMRIData2D);



%% RL
%% Set flags to determine which experimental runs will be processed and analyzed
% load rest fMRI data
dataNii=load_nii(['E:\S1200_Retest\S1200\Rest_FIX_NS\',sub,'\rfMRI_REST1_RL\rfMRI_REST1_RL_NASN_filtered.nii.gz']);

%% We identify and exclude volumes with significant head motion contamination 
% (exceeding the specified FD threshold) along with their immediate temporal 
% neighbors. Additionally, the dummy scans are discarded 
% to eliminate T1-equilibration effects and ensure the remaining time-series 
% reflects a stable, steady-state fMRI signal.
motionLabels= censorMotionHCP(['E:\S1200_Retest\S1200\Rest_FIX_NS\',sub,'\rfMRI_REST1_RL\Movement_Regressors.txt']);
retainedFramesIdx = find(motionLabels);
retainedFramesIdx(retainedFramesIdx < 1+NumDummyPoints) = []; % remove the dummy scans 
fMRIData4D=double(dataNii.img(:,:,:,retainedFramesIdx));

[N1,N2,N3,N4]=size(fMRIData4D);
voxelsize=dataNii.original.hdr.dime.pixdim(2);
TR=dataNii.original.hdr.dime.pixdim(5);

%% spatial smoothinng 
for ii=1:N4
    fMRIData4D(:,:,:,ii)=imgaussfilt3(fMRIData4D(:,:,:,ii).*brainMask3D,FWHM/voxelsize/2.355);  %三维图像的三维高斯滤波
end

%% preprocessing of fMRI time courses, including (1) detrend (2) bandpass filter (3)normalize
fMRIData2D=reshape(fMRIData4D,[N1*N2*N3 N4]);
fMRIData2D=fMRIData2D(:,:)';
fMRIData2D=detrend(fMRIData2D,2);

%% performs global signal regression
grayMatterMaskFile='../Masks/GrayMatterMask.nii';
fMRIData2D=regressGlobalSignal(fMRIData2D,grayMatterMaskFile);


%% The combined temporal filtering strategy effectively removes noise without introducing phase delays into the BOLD signal.
fMRIData2D = applyBandpassFilter(fMRIData2D, TR,freqBand(1),freqBand(2));

%% normalize the fMRI signals
fMRIData2D_RL=zscoreTimeSeries(fMRIData2D);
fMRIData2D_total=[fMRIData2D_LR;fMRIData2D_RL];


%% remove the baselines of fMRI data
infoFile=[outputpath_sub,'\','sub',sub,'_SW_baseLineInfo.mat'];
load(infoFile);
[fMRIData2D_total] = nullModelBaselineCorrect(fMRIData2D_total,baseLineInfo);
fMRIData2D_total(isnan(fMRIData2D_total))=0;
fMRIData2D_total(isinf(fMRIData2D_total))=0;

%% Define Sliding Window and Calculate Total Steps
slidingWindow=1:ceil(durHRF/TR); 
numWindFrame=size(fMRIData2D_total,1)-length(slidingWindow)+1;  %the number of the frames that can be calculated

% Validate sample size sufficiency against permutation parameters
% If the number of remaining clean timepoints (numWindFrame) is lower than the 
% target permutation count (numNullModels), the partition boundary is dynamically 
% truncated to match the data limitations, ensuring strict statistical validity.
if numWindFrame < numNullModels
    % Cache the original user setting for clear diagnostic reporting
    originalNullModels = numNullModels;
    
    % Force alignment between permutation count and available degrees of freedom
    numNullModels = numWindFrame;
    
    % Display a detailed warning tracking the adaptive parameter correction
    warning('iBAM:DataLimitation', ...
        ['\n[WARNING] The number of available clean frames (%d) is fewer than ' ...
         'the requested null models (%d).\nTo prevent matrix underflow, ' ...
         'the number of null models has been adjusted to: %d.\n'], ...
        numWindFrame, originalNullModels, numNullModels);
end

randWindFrame=randsample(1:numWindFrame,numNullModels);


fprintf('Sliding Window Calculation Initiated...\n');
prevPercent = -1; % Set to -1 so that 0% prints immediately at start
msgTemplate = 'Progress: %3d%%';
msgLength = 0;
FrameCounter=0;
ModeIntensity2D=[];
scaleFactor=10;


for n=randWindFrame 

% Display processing progress for each subject's null model MI mapping (Throttled at 1% increments)
    FrameCounter=FrameCounter+1;
%     display(['Frame #',num2str(FrameCounter)]);
    currentPercent = floor((FrameCounter / numNullModels) * 100);
    if currentPercent > prevPercent
        if msgLength > 0
            fprintf(repmat('\b', 1, msgLength)); 
        end
        newMsg = sprintf(msgTemplate, currentPercent);
        fprintf('%s', newMsg);
        msgLength = length(newMsg); 
        
        prevPercent = currentPercent;
    end

index_SlidingWindow=slidingWindow+n-1;
dataInSlidingWindow=fMRIData2D_total(index_SlidingWindow,:);
MI=zeros(length(brainIndex),1);

parfor ind=1:length(brainIndex)

    % The two variables, FSSWIndices_Dictionary and FSSWWeights_Dictionary, are loaded from the FSSW_K98.mat file.
    FSSWIndices  = FSSWIndices_Dictionary(:,ind);  % linear indices of FSSW voxels in mask
    FSSWWeights = FSSWWeights_Dictionary(:,ind);   
    if isempty(find(FSSWIndices ))
       continue; 
    end 

    % Find the Center Location (Maximum Weight) for Each FSSW
    [~,centerPosition]=maxk(FSSWIndices,1); %

    %% a functional topography informed principal mode analysis (FT-PMA) was implemented via singular value decomposition
    %  (SVD) to extract and quantify the LSM of fMRI signals within FSSW. Here, the conventional SVD
    %  is modified to allow the integration of the brain s functional topography encoded in FSSW, wherein
    %  the SVD was performed on time series weighted by the FSSW, effectively prioritizing signal
    %  contributions from voxels proximal to the central voxel, and thus sharpening the spatial specificity
    %  of brain activation detection. The first principal component of SVD in the temporal domain
    %  represents the LSM, while its corresponding singular value measures MI and serves as a functional
    %  metric for neural activation intensity
    X=dataInSlidingWindow(:,FSSWIndices);
    W=repmat(FSSWWeights',size(X,1),1); 
    X=X.*W;
    [U,S,V]=svds(double(X),1,'largest');

    %% SIGN AMBIGUITY RESOLUTION FOR LSM & MI:
    %   To resolve the inherent sign ambiguity in eigen-decomposition, a seamless
    %   polarity alignment is implemented. We first confine the distribution
    %   coefficient at the window center (i) to a positive value, where the
    %   localized spatial eigenvector v1 is adjusted as v1 = sign(v1(i)) * v1.
    %   Accordingly, the polarity of the dominant latent signal mode (LSM) is
    %   aligned as u1 = sign(v1(i)) * u1 to maintain local consistency. Following
    %   this local alignment, we calculate the sign of the sum of the first spatial
    %   eigenvector to determine the global polarity of the LSM, which is
    %   subsequently applied to the first eigenvalue to unambiguously determine
    %   the true direction and polarity of brain activation across the entire network.
    MI(ind)=S.*sign(sum(sign(V(centerPosition)).*U))/sqrt(size(X,1))/scaleFactor;

end
    ModeIntensity2D=[ModeIntensity2D,MI];
end

fprintf('\nCalculation of MI map is completed successfully.\n');

MI_save=zeros(prod(brainMaskDim3D),numNullModels);%
MI_save(brainIndex,:)=ModeIntensity2D;
MI_save=reshape(MI_save,[brainMaskDim3D,numNullModels]);




%% The calculated MI maps were exported into NIfTI format, using the previously loaded fMRI.nii as a spatial template to retain the original geometric metadata.
templateNii=dataNii;
templateNii.img=templateNii.img.*0;
% template_nii.img(:,:,:,1+dumypoints:wind_frame+dumypoints)=single(MI_save);
templateNii.img=(single(MI_save));
templateNii.hdr.dime.dim(5)=numNullModels;
templateNii.original.hdr.dime.dim(5)=numNullModels;
f_name=[outputpath_sub,'\','Emotion_sub',sub,'_SW_MI_NullModels'];
% save_nii(template_nii,[outputpath '\sub' f_name,'_original_98nei_sub_001_run_1.nii.gz'])
save_nii(templateNii,[f_name, '.nii'])



% save_nii(template_nii,[f_name,'.nii']);
result_dir=outputpath;
str=['D:/workbench/bin_windows64/wb_command -volume-to-surface-mapping ',[f_name,'.nii '],'../Masks/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii ',[result_dir,'/BrainL.shape.gii'],' -trilinear'];
eval(['!',str]);
str=['D:/workbench/bin_windows64/wb_command -volume-to-surface-mapping ',[f_name,'.nii '],'../Masks/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii ',[result_dir,'/BrainR.shape.gii'],' -trilinear'];
eval(['!',str]);
str=['D:/workbench/bin_windows64/wb_command -cifti-create-dense-timeseries ',[f_name,'.dtseries.nii '],'-left-metric ',[result_dir,'/BrainL.shape.gii'],' -roi-left  ../Masks/L.atlasroi.32k_fs_LR.shape.gii -right-metric ',[result_dir,'/BrainR.shape.gii'],' -roi-right ../Masks/R.atlasroi.32k_fs_LR.shape.gii -timestep 1'];
eval(['!',str]);
delete([result_dir,'/BrainL.shape.gii']);
delete([result_dir,'/BrainR.shape.gii']);

end
delete(gcp('nocreate'));

