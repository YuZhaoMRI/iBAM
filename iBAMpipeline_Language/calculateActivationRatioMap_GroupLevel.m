
% DESCRIPTION:
%       This function calculates the voxel-wise activation ratio map separately 
%       for each experimental condition at group level.


clear;clc;

inputpath='F:\Language_task_iBAM_results';
outputpath=inputpath;
fdrAlpha = 0.05;

%% load brain mask
maskFileName='../Masks/MNI152_T1_2mm_brain.nii';
[brainMask3D,brainMaskDim3D,brainIndex]=getBrainMaskInfo(maskFileName);


AR_group=zeros(brainMaskDim3D);

%% get the participant folders
contents = dir(inputpath);
sublist = contents([contents.isdir]);  % Get all folders
sublist = sublist(~ismember({sublist.name}, {'.', '..'}));  % Exclude current and parent directories

numSub=0;

for m=1:length(sublist)
display(['Subject #',num2str(m)])
sub=char(sublist(m).name);
inputpath_sub=[inputpath,'\',sub];

%%
niiARmapFileName = [inputpath_sub,'\','sub',sub,'_AR_map.nii'];
if exist(niiARmapFileName, 'file') == 0
   continue;
end



%%
niiData=load_nii(niiARmapFileName);
AR_individual=single(niiData.img);


numSub=numSub+1;
AR_group=AR_group+AR_individual;


end


AR_group=AR_group/numSub;

niiData=load_nii('../Masks/niitemplate.nii');
niiData.img=single(AR_group);
niiData.hdr.dime.dim(5)=1;
niiData.original.hdr.dime.dim(5)=1;
fname=[inputpath,'\','GroupLevel_AR_map.nii'];
save_nii(niiData,fname);
% 

