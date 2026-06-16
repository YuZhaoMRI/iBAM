%genModeIntensityMap_SW Calculate Mode Intensity Maps using iBAM with a Sliding Window
%
% PURPOSE:
%   This function implements the implicit Brain Activity Mapping (iBAM) 
%   algorithm to probe implicit, non-task-timing-locked brain activities 
%   by calculating Mode Intensity (MI) maps using a temperal sliding window approach.
%
% ALGORITHM MECHANISM:
%   iBAM operates on a voxel-wise basis to detect implicit brain activity. 
%   Prior to executing this algorithm, a Functional Structure-encoded Spatial 
%   Window (FSSW) incorporating group-level functional topography must be 
%   pre-computed for each voxel. This FSSW structure is loaded as a spatial 
%   prior during runtime execution. Within these defined FSSWs, the algorithm 
%   extracts the Latent Signal Mode (LSM) from the preprocessed fMRI time 
%   series. By identifying these LSMs and quantifying their corresponding Mode 
%   Intensities (MIs), iBAM enables robust detection of implicit brain 
%   activations without the need for task-timing-based signal modeling.
%
% MANDATORY PREPROCESSING REQUIREMENTS:
%   The input fMRI time series MUST be rigorously preprocessed beforehand. 
%   The pipeline requires data that has undergone:
%     1. Standard minimal preprocessing via 'fMRIprep' or the 'HCP minimal 
%        preprocessing pipeline'.
%     2. Spatial normalization/registration to the standard MNI 152 space.
%     3. Advanced denoising via ICA-based filtering (see more details in
%        the Methods section in our paper) to remove structured spatial-temporal artifacts.
%
% INPUTS:
%   Path to task fmri data  - preprocessed fMRI time series in MNI 152 space.
%                
% OUTPUTS:
%   Path to MI maps    - Mode Intensity (MI) maps.
%
% REFERENCE:
%   For detailed methodology, particularly regarding the preprocessing pipeline 
%   and FSSW configuration, please refer to the Methods section of:
%   "Implicit brain activity mapping reveals language as a higher-order 
%   strategic tool in human reasoning."
%
% Author: Yu Zhao
% Date: May 2026
% Version: 1.0
%==========================================================================


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


%% set parameters
freqBand       = [0.01 0.18]; % Temporal bandpass filter range (Hz)
FWHM           = 4;           % Gaussian kernel FWHM for spatial smoothing (mm) 

% % Remove the Dummy Scans
% 1. MAGNETIZATION STEADY-STATE: Although prospective dummy scans are typically 
%    discarded automatically by the scanner console during acquisition, non-steady-state 
%    volumes (T1 effects) may still leak into the initial time series.
% 2. FILTERING EDGE EFFECTS: Temporal bandpass filtering can heavily propagate initial 
%    signal transient spikes into subsequent frames, inducing unwanted edge artifacts 
%    and transient signal fluctuations (ringing). This is especially severe if the 
%    first few volumes suffer from non-steady-state magnetization.
% 3. PARADIGM DESIGN RECOMMENDATION: For custom/in-house datasets, it is highly 
%    recommended to configure 10+ prospective dummy scans on the scanner, insert a 
%    30-second resting baseline prior to task onset, and manually discard the first 
%    5 to 10 volumes. Delaying the actual task onset guarantees that the filter's 
%    transient response settles completely within the baseline period, safeguarding 
%    the integrity of task-related implicit brain dynamics.
NumDummyPoints = 5;           % Number of dummy volumes removed for T1 equilibration

durRest        = 0 * 60;      % Total duration of resting-state fMRI acquisition (s)

% Duration of the canonical HRF (seconds). 
% Note: If increasing this value to improve the signal-to-noise 
% ratio (SNR) of MI maps, ensure it is updated concurrently in genModeIntensityMap_SW.m.
durHRF         = 15;         


% . Construct the fully qualified filename for data preservation
saveFileName = fullfile(outputpath, 'pipeline_parameters.mat');
%  Serialize and export the specific parameter variables
save(saveFileName, 'freqBand', 'FWHM', 'NumDummyPoints', 'durHRF', '-v7.3');
fprintf('Configuration parameters successfully exported to: %s\n', saveFileName);
% PARAMETER CONSISTENCY FOR NULL MODEL CONFIGURATION:
%   These saved parameters will be reloaded and implemented during the subsequent 
%   null model construction phase. Preserving these exact configuration values 
%   is mathematically imperative to guarantee strict methodological consistency 
%   between the empirical null distribution generation and the primary analysis.



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



% Iterate through each subject sequentially
for isub=1:length(sublist)

display(['generate MI maps for subject #',num2str(isub)])
sub=char(sublist(isub).name);

outputpath_sub=[outputpath,'\',sub];
if ~isfolder(outputpath_sub)
    mkdir(outputpath_sub)
end

%% load fMRI volumes (left-right phase-ecoding)
% for test dataset
dataNii=load_nii([inputpath ,sub,'\MNINonLinear\Results\tfMRI_EMOTION_LR\tfMRI_EMOTION_LR_SNfiltered.nii.gz']);
% for retest dataset
%dataNii=load_nii(['E:\S1200_Retest\Retest\Emotion\',sub,'\MNINonLinear\Results\tfMRI_EMOTION_LR\tfMRI_EMOTION_LR.nii.gz']);

TR=dataNii.original.hdr.dime.pixdim(5);
fMRIData4D=double(dataNii.img(:,:,:,1+NumDummyPoints:end));% remove the dummy scans. note 

[N1,N2,N3,N4]=size(fMRIData4D);

voxelsize=dataNii.original.hdr.dime.pixdim(2);

%% Spatial smoothinng 
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


%% load structured noise (SN)-filtered fMRI volumes (right-left phase-ecoding) RL
%data_nii=load_nii(['E:\S1200_Retest\S1200\Emotion\',sub,'\MNINonLinear\Results\tfMRI_EMOTION_RL\tfMRI_EMOTION_RL.nii.gz']);
 data_nii=load_nii([inputpath ,sub,'\MNINonLinear\Results\tfMRI_EMOTION_RL\tfMRI_EMOTION_RL_SNfiltered.nii.gz']);
TR=dataNii.original.hdr.dime.pixdim(5);
fMRIData4D=double(dataNii.img(:,:,:,1+NumDummyPoints:end));% remove the dummy and resting-state fMRI scans at the end 
[N1,N2,N3,N4]=size(fMRIData4D);
voxelsize=dataNii.original.hdr.dime.pixdim(2);
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

%% get the timings of fMRI tasks
eventFile_LR=[inputpath ,sub,'\MNINonLinear\Results\tfMRI_EMOTION_LR\EVs\'];
Descent_dur=10;% second,Descent duration of HRF
TaskTimings_LR = generateHCPTaskTimingVector(eventFile_LR, TR, Descent_dur, NumDummyPoints);
fMRIData2D_LR(size(TaskTimings_LR,2)+1:end,:)=[];
eventFile_RL=[inputpath ,sub,'\MNINonLinear\Results\tfMRI_EMOTION_RL\EVs\'];
Descent_dur=10;% second,Descent duration of HRF
TaskTimings_RL = generateHCPTaskTimingVector(eventFile_RL, TR, Descent_dur, NumDummyPoints);
fMRIData2D_RL(size(TaskTimings_RL,2)+1:end,:)=[];
TaskTimings=[TaskTimings_LR,TaskTimings_RL];

%% correct baselines 
fMRIData2D_Total=[fMRIData2D_LR;fMRIData2D_RL];
[fMRIData2D_Total, baseLineInfo] = baselineCorrect(fMRIData2D_Total);
fMRIData2D_Total(isnan(fMRIData2D_Total))=0;
baseLineInfoFile=[outputpath_sub,'\','sub',sub,'_SW_baseLineInfo.mat'];
save(baseLineInfoFile,'baseLineInfo')



%% Define Sliding Window and Calculate Total Steps
% Convert canonical HRF duration into the number of volumes/frames (Window Size)
slidingWindow = 1:ceil(durHRF / TR); 

% Calculate the total number of sliding windows (steps) across the time series
numTotalWindows  = size(fMRIData2D_Total, 1) -length(slidingWindow) + 1;
scaleFactor=10;
ModeIntensity2D=[];

fprintf('Sliding Window Calculation Initiated...\n');
prevPercent = -1; % Set to -1 so that 0% prints immediately at start
msgTemplate = 'Progress: %3d%%';
msgLength = 0;


for n=1:numTotalWindows  

% Display processing progress for each subject's MI mapping (Throttled at 1% increments)
    currentPercent = floor((n / numTotalWindows) * 100);
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
dataInSlidingWindow=fMRIData2D_Total(index_SlidingWindow,:);
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

MI_save=zeros(prod(brainMaskDim3D),numTotalWindows);%
MI_save(brainIndex,:)=ModeIntensity2D;
MI_save=reshape(MI_save,[brainMaskDim3D,numTotalWindows]);

%% label the corner of images with task timings;  
% -1 indicates neutual task-on periods; 
%  1 indicates fear task-on periods;and 0 indicates rest.
TaskTimings(numTotalWindows+1:end)=[];
taskLabels=find(TaskTimings==1);
MI_save(:,1:2,1:2,taskLabels)=1;
MI_save(1:2,:,1:2,taskLabels)=1;
MI_save(1:2,1:2,:,taskLabels)=1;
taskLabels=find(TaskTimings==-1);
MI_save(:,1:2,1:2,taskLabels)=-1;
MI_save(1:2,:,1:2,taskLabels)=-1;
MI_save(1:2,1:2,:,taskLabels)=-1;


%% The calculated MI maps were exported into NIfTI format, using the previously loaded fMRI.nii as a spatial template to retain the original geometric metadata.
templateNii=dataNii;
templateNii.img=templateNii.img.*0;
templateNii.img=(single(MI_save));
templateNii.hdr.dime.dim(5)=numTotalWindows;
templateNii.original.hdr.dime.dim(5)=numTotalWindows;
f_name=[outputpath_sub,'\','Emotion_sub',sub,'_SW_MI'];
save_nii(templateNii,[f_name, '.nii'])

%% Map volumetric data onto the brain surface via HCP Workbench
result_dir=outputpath_sub;
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

