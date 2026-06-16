%   DESCRIPTION:
%    computes the group-level average MI map from a cohort of subjects.
%    Specifically, the function performs a two-level averaging process:
%     1. Individual Level: For each subject's MI maps, it identifies all 
%        frames that fall within the task blocks of the specified condition 
%        and averages these frames together.
%     2. Group Level: It then averages these individual-level mean maps 
%        across all subjects to obtain the final group-level average.

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


MI_neut_group=zeros(brainMaskDim3D);
MI_fear_group=zeros(brainMaskDim3D);


numFrame_neut=0;
numFrame_fear=0;
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



taskInd=find(TaskTimings==-1);
temperalmap1=squeeze(sum((vol4D_observed(:,:,:,taskInd)),4));
MI_neut_group=MI_neut_group+temperalmap1;
numFrame_neut=numFrame_neut+length(taskInd);

taskInd=find(TaskTimings==1);
temperalmap2=squeeze(sum((vol4D_observed(:,:,:,taskInd)),4));
MI_fear_group=MI_fear_group+temperalmap2;
numFrame_fear=numFrame_fear+length(taskInd);


end


MI_neut_group=MI_neut_group/numFrame_neut;
MI_fear_group=MI_fear_group/numFrame_fear;

niiData=load_nii('../Masks/niitemplate.nii');
niiData.img=single(MI_neut_group);
niiData.hdr.dime.dim(5)=1;
fname=[outputpath,'\','GroupLevel_neut_MI_map.nii'];
save_nii(niiData,fname);

niiData.img=single(MI_fear_group);
niiData.hdr.dime.dim(5)=1;
fname=[outputpath,'\','GroupLevel_fear_MI_map.nii'];
save_nii(niiData,fname);
