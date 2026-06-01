clear;clc;

%% outputpath
outputpath='F:\Language_task_results_v5\Reference';
if ~isfolder(outputpath)
    mkdir(outputpath)
end

%% load brain mask
maskFileName='../Masks/MNI152_T1_2mm_brain.nii';
[brainMask3D,brainMaskDim3D,brainIndex]=getBrainMaskInfo(maskFileName);


%% get the participant folders
inputpath = 'F:\Language_task_results_v5\Reference';
% Get all contents in the specified path (both files and folders)
contents = dir(inputpath);
% Filter out folders (excluding '.' and '..')
sublist = contents([contents.isdir]);  % Get all folders
sublist = sublist(~ismember({sublist.name}, {'.', '..'}));  % Exclude current and parent directories



delete(gcp('nocreate'));
parpool(64);

total_frame=0;
total_cmap=0;
for m=1:length(sublist)
display(['iBAM for subject #',num2str(m)])
sub=char(sublist(m).name);

outputpath_sub=[outputpath,'\',sub];
if ~isfolder(outputpath_sub)
    mkdir(outputpath_sub)
end
f_name=[outputpath_sub,'\','run-3_sub',sub,'_SW_MI_NullModels.nii'];
% f_name=[outputpath_sub,'\','Emotion_sub',sub,'_SW_MI_NullModels.nii'];
data_nii=load_nii(f_name);
vol4D_null=single(data_nii.img);
[N1,N2,N3,num_nullmodels]=size(vol4D_null);
vol4D_null=reshape(vol4D_null,[N1*N2*N3 num_nullmodels]);
Null_models=vol4D_null(brainIndex,:);
clear data_nii

f_name=[outputpath_sub,'\','run-3_sub',sub,'_SW_MI_reference.nii'];
% f_name=[outputpath_sub,'\','Emotion_sub',sub,'_SW_MI.nii'];
data_nii=load_nii(f_name);
vol4D_observed=single(data_nii.img(:,:,:,1:end));
% Tasktimings=squeeze(vol4D_observed(1,1,1,:));
[N1,N2,N3,num_wid]=size(vol4D_observed);
vol4D_observed=reshape(vol4D_observed,[N1*N2*N3 num_wid]);
Observed_data=vol4D_observed(brainIndex,:);

left_area=zeros(length(brainIndex),num_wid);
parfor i=1:length(brainIndex)
%       [f,xi]=ksdensity(Null_models(i,:),'Function','cdf', 'Bandwidth',0.002,'NumPoints',num_nullmodels);
%       [f,xi]=ksdensity(Null_models(i,:),'Function','cdf', 'Bandwidth',0.005,'NumPoints',1*num_nullmodels);
    null_data=Null_models(i,:);      
    outliers=isoutlier(null_data);
    null_idx=find(outliers==0);
    [f,xi]=ksdensity(Null_models(i,null_idx),'Function','cdf');
    [~,index]=min(abs(1*Observed_data(i,:)-xi(:)),[],1);
    left_area(i,:)=f(index);
end
% left_area(left_area<0.0000001)=0.0000001;
% left_area(left_area>0.9999999)=0.9999999;




Z_map=zeros(N1*N2*N3,num_wid);
Z_map(brainIndex,:)=norminv(left_area);
Z_map=reshape(Z_map,N1,N2,N3,num_wid);

P_map=zeros(N1*N2*N3,num_wid);
P_map(brainIndex,:)=1-left_area;
P_map=reshape(P_map,N1,N2,N3,num_wid);


for k=1:num_wid
P_vector=reshape(P_map(:,:,:,k),[N1*N2*N3 1]);
P_vector=sort(P_vector,'ascend');
u=spm_uc_FDR(0.05,[1 Inf],'Z',1,P_vector,1);
Z_map_3D=Z_map(:,:,:,k);
Z_map_3D(Z_map_3D<u)=0;
Z_map(:,:,:,k)=Z_map_3D;
end

% %
% taskLabels=find(Tasktimings==1);
% Z_map(:,1:2,1:2,taskLabels)=1;
% Z_map(1:2,:,1:2,taskLabels)=1;
% Z_map(1:2,1:2,:,taskLabels)=1;
% taskLabels=find(Tasktimings==-1);
% Z_map(:,1:2,1:2,taskLabels)=-1;
% Z_map(1:2,:,1:2,taskLabels)=-1;
% Z_map(1:2,1:2,:,taskLabels)=-1;
% %


% Z_map = removeSmallObjects3D(Z_map, 10);
data_nii.img=data_nii.img.*0;
data_nii.img(:,:,:,1:num_wid)=Z_map;
f_name=[outputpath_sub,'\','sub',sub,'_SW_Zmap'];
save_nii(data_nii,[f_name, '.nii']);

% mTasktimings=Tasktimings;
% for k=1:length(Tasktimings)-1
% if Tasktimings(k)==0&&Tasktimings(k+1)==1
%    mTasktimings(k+1)=0;
% end
% end
% for k=1:length(Tasktimings)-1
% if Tasktimings(k)==1&&Tasktimings(k+1)==0
%    mTasktimings(k+1)=0;
% end
% end

% % neut C map
% cZ_map=Z_map;
% cZ_map(cZ_map>0)=1;
% taskind=find(mTasktimings==1);
% cmap=squeeze(sum(cZ_map(:,:,:,taskind),4));
% C=(single(cmap/length(taskind)));
% C(C<0.0)=0;
% data_nii.img=C;
% data_nii.hdr.dime.dim(5)=1;
% f_name=[outputpath_sub,'\','Logical_sub',sub,'_Cmap'];

% save_nii(data_nii,[f_name, '.nii']);
% total_cmap=total_cmap+cmap;
% total_frame=total_frame+length(taskind);

end


% 
% C=(single(total_cmap/total_frame));
% C(C<0.2)=0;
% data_nii.img=C;
% data_nii.hdr.dime.dim(5)=1;
% f_name=['Group_logical_sh2Cmap'];
% save_nii(data_nii,[f_name, '.nii']);
% C=(single(total_cmap/total_frame));
% 
% C(C<0.5)=0;
% data_nii.img=C;
% data_nii.hdr.dime.dim(5)=1;
% f_name=['Group_logical_sh5Cmap'];
% save_nii(data_nii,[f_name, '.nii']);
% 
% 
% %%
% f_name=['Visual_sub2025082302_Cmap.nii'];
% % f_name=[outputpath_sub,'\','Emotion_sub',sub,'_SW_MI.nii'];
% data_nii=load_nii(f_name);
% img=data_nii.img;
% img(img<0.2)=0;
% data_nii.img=img;
% data_nii.hdr.dime.dim(5)=1;
% f_name=['Visual_sub2025082302_sh2Cmap'];
% save_nii(data_nii,[f_name, '.nii']);
% 
% img(img<0.5)=0;
% data_nii.img=img;
% data_nii.hdr.dime.dim(5)=1;
% f_name=['Visual_sub2025082302_sh5Cmap'];
% save_nii(data_nii,[f_name, '.nii']);