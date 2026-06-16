%   DESCRIPTION:
%    computes the group-level average MI map from a cohort of subjects.
%    Specifically, the function performs a two-level averaging process:
%     1. Individual Level: For each subject's MI maps, it identifies all 
%        frames that fall within the task blocks of the specified condition 
%        and averages these frames together.
%     2. Group Level: It then averages these individual-level mean maps 
%        across all subjects to obtain the final group-level average.

clear;clc;

inputpath='F:\Language_task_iBAM_results';
outputpath=inputpath;


%% load brain mask
maskFileName='../Masks/MNI152_T1_2mm_brain.nii';
[brainMask3D,brainMaskDim3D,brainIndex]=getBrainMaskInfo(maskFileName);


%% get the participant folders
contents = dir(inputpath);
sublist = contents([contents.isdir]);  % Get all folders
sublist = sublist(~ismember({sublist.name}, {'.', '..'}));  % Exclude current and parent directories



MI_group=zeros(brainMaskDim3D);
numFrame_task=0;

for m=1:length(sublist)
display(['Subject #',num2str(m)])
sub=char(sublist(m).name);

outputpath_sub=[outputpath,'\',sub];
f_name=[outputpath_sub,'\','sub',sub,'_SW_MI.nii'];
if exist(f_name, 'file') == 0
   continue;
end
data_nii=load_nii(f_name);
vol4D_observed=(single(data_nii.img(:,:,:,1:end)));
TaskTimings=squeeze(vol4D_observed(1,1,1,:));
[N1,N2,N3,num_wid]=size(vol4D_observed);


taskInd=find(TaskTimings==1);
temperalmap2=squeeze(sum((vol4D_observed(:,:,:,taskInd)),4));
MI_group=MI_group+temperalmap2;
numFrame_task=numFrame_task+length(taskInd);


end



MI_group=MI_group/numFrame_task;

niiData=load_nii('../Masks/niitemplate.nii');
niiData.img=single(MI_group);
niiData.hdr.dime.dim(5)=1;
niiData.original.hdr.dime.dim(5)=1;
fname=[outputpath,'\','GroupLevel_MI_map.nii'];
save_nii(niiData,fname);
