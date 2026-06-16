%   DESCRIPTION:
%       This function calculates the voxel-wise activation ratio map separately 
%       for each experimental condition within the Human Connectome Project (HCP) 
%       Emotion Processing Task (Neutral and Fear blocks) at individual level.
%
%       THEORETICAL DEFINITION:
%       Following the iBAM statistical framework, the within-block activation rate 
%       is mathematically defined as the proportion of time points showing significant 
%       activation relative to the cumulative duration of all frames enclosed within 
%       the blocks of a specific task condition.

clear;clc;
dumypoints=5; % remove the first N points
inputpath='E:\HCP_Emotion_results_0605\Emotion\Test';
outputpath=inputpath;
fdrAlpha = 0.05;

%% load brain mask
maskFileName='../Masks/MNI152_T1_2mm_brain.nii';
[brainMask3D,brainMaskDim3D,brainIndex]=getBrainMaskInfo(maskFileName);


%% get the participant folders
contents = dir(inputpath);
sublist = contents([contents.isdir]);  % Get all folders
sublist = sublist(~ismember({sublist.name}, {'.', '..'}));  % Exclude current and parent directories


for m=1:length(sublist)
display(['Subject #',num2str(m)])
sub=char(sublist(m).name);

inputpath_sub=[inputpath,'\',sub];

niiZmapFileName = sprintf('%s\\Emotion_sub%s_SW_Zmap_thresholded_FDR%.2f.nii', ...
    inputpath_sub, sub, fdrAlpha);

if exist(niiZmapFileName, 'file') == 0
   continue;
end
niiData=load_nii(niiZmapFileName);
vol4D_observed=(single(niiData.img(:,:,:,1:end)));
TaskTimings=squeeze(vol4D_observed(1,1,1,:));
[N1,N2,N3,num_wid]=size(vol4D_observed);


[neutTaskTiming,fearTaskTiming]=seperateHCPTaskBlocks(TaskTimings);

vol4D_observed(vol4D_observed>0)=1;
taskInd=find(neutTaskTiming==1);
ActivationRatioMap=squeeze(sum(abs(vol4D_observed(:,:,:,taskInd)),4)/length(taskInd));
niiData.img=(single(ActivationRatioMap));
niiData.hdr.dime.dim(5)=1;
fname=[inputpath_sub,'\','Emotion_sub',sub,'_neut_AR_map.nii'];
save_nii(niiData,fname);



taskInd=find(fearTaskTiming==1);
ActivationRatioMap=squeeze(sum(abs(vol4D_observed(:,:,:,taskInd)),4)/length(taskInd));
ActivationRatioMap(ActivationRatioMap<0.5)=0;
niiData.img=single(ActivationRatioMap);
niiData.hdr.dime.dim(5)=1;
fname=[inputpath_sub,'\','Emotion_sub',sub,'_fear_AR_map.nii'];
save_nii(niiData,fname);


end