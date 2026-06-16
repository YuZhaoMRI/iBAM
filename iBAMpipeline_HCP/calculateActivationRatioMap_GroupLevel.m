
% DESCRIPTION:
%       This function calculates the voxel-wise activation ratio map separately 
%       for each experimental condition within the Human Connectome Project (HCP) 
%       Emotion Processing Task (Neutral and Fear blocks) at group level.


clear;clc;
dumypoints=5; % remove the first N points
inputpath='E:\HCP_Emotion_results_0605\Emotion\Test';
outputpath=inputpath;
fdrAlpha = 0.05;

%% load brain mask
maskFileName='../Masks/MNI152_T1_2mm_brain.nii';
[brainMask3D,brainMaskDim3D,brainIndex]=getBrainMaskInfo(maskFileName);


AR_neut_group=zeros(brainMaskDim3D);
AR_fear_group=zeros(brainMaskDim3D);

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
niiZmapFileName_neut = [inputpath_sub,'\','Emotion_sub',sub,'_neut_AR_map.nii'];
if exist(niiZmapFileName_neut, 'file') == 0
   continue;
end

niiZmapFileName_fear = [inputpath_sub,'\','Emotion_sub',sub,'_fear_AR_map.nii'];
if exist(niiZmapFileName_fear, 'file') == 0
   continue;
end

%%
niiData=load_nii(niiZmapFileName_neut);
AR_neut_individual=single(niiData.img);
niiData=load_nii(niiZmapFileName_fear);
AR_fear_individual=single(niiData.img);

numSub=numSub+1;
AR_neut_group=AR_neut_group+AR_neut_individual;
AR_fear_group=AR_fear_group+AR_fear_individual;

end


AR_neut_group=AR_neut_group/numSub;
AR_fear_group=AR_fear_group/numSub;

niiData=load_nii('../Masks/niitemplate.nii');
niiData.img=single(AR_neut_group);
niiData.hdr.dime.dim(5)=1;
fname=[inputpath,'\','GroupLevel_neut_AR_map.nii'];
save_nii(niiData,fname);

niiData.img=single(AR_fear_group);
niiData.hdr.dime.dim(5)=1;
fname=[inputpath,'\','GroupLevel_fear_AR_map.nii'];
save_nii(niiData,fname);
