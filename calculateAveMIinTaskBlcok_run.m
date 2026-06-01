
clear;clc;
dumypoints=5; % remove the first N points
input_dir='F:\Language_task_results_v5\Visual';
test_label=[ ];
run_label='run-3';


%% outputpath
outputpath='F:\Language_task_results_v5\Visual';
if ~isfolder(outputpath)
    mkdir(outputpath)
end


%% load brain mask
maskFileName='../Masks/MNI152_T1_2mm_brain.nii';
[brainMask3D,brainMaskDim3D,brainIndex]=getBrainMaskInfo(maskFileName);


%% get the participant folders
inputpath = 'F:\Language_task_results_v5\Visual';
% Get all contents in the specified path (both files and folders)
contents = dir(inputpath);
% Filter out folders (excluding '.' and '..')
sublist = contents([contents.isdir]);  % Get all folders
sublist = sublist(~ismember({sublist.name}, {'.', '..'}));  % Exclude current and parent directories


total_MI=0;
num_sub=0;
for m=1:length(sublist)
display(['iBAM for subject #',num2str(m)])
sub=char(sublist(m).name);



outputpath_sub=[outputpath,'\',sub];
if ~isfolder(outputpath_sub)
    mkdir(outputpath_sub)
end

f_name=[outputpath_sub,'\',run_label,'_sub',sub,'_SW_MI_visual.nii'];
if exist(f_name, 'file') == 0
   continue;
end
data_nii=load_nii(f_name);
vol4D_observed=abs(single(data_nii.img(:,:,:,1:end)));
TaskTimings=squeeze(vol4D_observed(1,1,1,:));
[N1,N2,N3,num_wid]=size(vol4D_observed);



taskind=find(TaskTimings==1);
cmap=squeeze(sum(vol4D_observed(:,:,:,taskind),4));
data_nii.img=(single(cmap/length(taskind)));
data_nii.hdr.dime.dim(5)=1;
f_name=[outputpath_sub,'\','Visual_sub',sub,'_aveMI_map'];
save_nii(data_nii,[f_name, '.nii']);

total_MI=total_MI+cmap;
num_sub=num_sub+length(taskind);
end


data_nii.img=total_MI/num_sub;
data_nii.hdr.dime.dim(5)=1;
f_name=['Group_visual_sub',sub,'_aveMI_map'];
save_nii(data_nii,[f_name, '.nii']);



