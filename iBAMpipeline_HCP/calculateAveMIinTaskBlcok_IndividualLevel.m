
%   Calculates the average Mode Intensity (MI) at the individual level.
%   AVG_MI = CALCULATEAVERAGEMI(MIIMAGE, TASKBLOCKS, CONDITION) takes an
%   individual's MI image, identifies all frames that fall within the task
%   blocks of a specified condition, and returns the average of these frames.

clear;clc;
dumypoints=5; % remove the first N points
inputpath='E:\HCP_Emotion_results_0605\Emotion\Test';
outputpath=inputpath;


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

outputpath_sub=[outputpath,'\',sub];
f_name=[outputpath_sub,'\','Emotion_sub',sub,'_SW_MI.nii'];
if exist(f_name, 'file') == 0
   continue;
end
data_nii=load_nii(f_name);
vol4D_observed=(single(data_nii.img(:,:,:,1:end)));
TaskTimings=squeeze(vol4D_observed(1,1,1,:));
[N1,N2,N3,num_wid]=size(vol4D_observed);


[neutTaskTimings,fearTaskTimings]=seperateHCPTaskBlocks(TaskTimings);

taskind=find(TaskTimings==-1);
cmap=squeeze(sum(abs(vol4D_observed(:,:,:,taskind)),4));
data_nii.img=(single(cmap/length(taskind)));
data_nii.hdr.dime.dim(5)=1;
f_name=[outputpath_sub,'\','Emotion_sub',sub,'_neut_aveMI_map'];
save_nii(data_nii,[f_name, '.nii']);



taskind=find(TaskTimings==1);
cmap=squeeze(sum(abs(vol4D_observed(:,:,:,taskind)),4));
data_nii.img=(single(cmap/length(taskind)));
data_nii.hdr.dime.dim(5)=1;
f_name=[outputpath_sub,'\','Emotion_sub',sub,'_fear_aveMI_map'];
save_nii(data_nii,[f_name, '.nii']);


end


