clear;clc;

%% outputpath
outputpath='F:\Language_task_results_v5\Logical';
if ~isfolder(outputpath)
    mkdir(outputpath)
end

%% load brain mask
maskFileName='../Masks/MNI152_T1_2mm_brain.nii';
[brainMask3D,brainMaskDim3D,brainIndex]=getBrainMaskInfo(maskFileName);


%% get the participant folders
inputpath = 'F:\Language_task_results_v5\Logical';
% Get all contents in the specified path (both files and folders)
contents = dir(inputpath);
% Filter out folders (excluding '.' and '..')
sublist = contents([contents.isdir]);  % Get all folders
sublist = sublist(~ismember({sublist.name}, {'.', '..'}));  % Exclude current and parent directories




num_sub=0;
cmap_total=0;
for m=1:42%length(sublist)
display(['iBAM for subject #',num2str(m)])
sub=char(sublist(m).name);

outputpath_sub=[outputpath,'\',sub];
if ~isfolder(outputpath_sub)
    mkdir(outputpath_sub)
end

f_name=[outputpath_sub,'\','logical_sub',sub,'_Cmap.nii'];
% f_name=[outputpath_sub,'\','Emotion_sub',sub,'_SW_MI.nii'];
data_nii=load_nii(f_name);
cmap=single(data_nii.img);
cmap_total=cmap_total+cmap;
num_sub=num_sub+1;



end



C=(single(cmap_total/num_sub));
C(C<0.0)=0;
data_nii.img=C;
data_nii.hdr.dime.dim(5)=1;
f_name=['Group_logical_Cmap'];
save_nii(data_nii,[f_name, '.nii']);