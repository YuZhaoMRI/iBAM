clear;clc;
dumypoints=5; % remove the first N points
input_dir='F:\MIQ_zhaoyu\Codes_for_HCP_task\Results';
test_label=[ ];

% %% outputpath
% outputpath='E:\HCP_retest_results\Emotion\results_sliding_window';
% if ~isfolder(outputpath)
%     mkdir(outputpath)
% end


%% outputpath
outputpath='E:\Language_task_results_v3\Logical';
if ~isfolder(outputpath)
    mkdir(outputpath)
end


%% load brain mask
 maskBrain_nii=load_nii('../Masks/MNI152_T1_2mm_brain.nii');
%maskBrain_nii=load_nii('../Masks/BN_Atlas_246_2mm.nii');
maskBrain3D=double(maskBrain_nii.img);
dim3DImag=maskBrain_nii.hdr.dime.dim(2:4);
maskBrain1D=reshape(maskBrain3D,[prod(dim3DImag) 1]);
maskBrain3D(maskBrain3D>0)=1;
voxelsize=2;
indexBrain=find(maskBrain1D>1);
numNeighbors=2*7.^2;% 3-7



 sublist=["601127"];
for m=1:1%length(sublist)

sub=char(sublist(m));

% temp_file=fullfile([input_dir '\' sub],'**\*null_MI.nii');
% temp_file=dir( temp_file);
% data_nii=load_nii([temp_file.folder, '\',temp_file.name]);

% data_nii=load_nii('E:\HCP_retest_results\Emotion\results_sliding_window\sub601127500null_MI.nii');
data_nii=load_nii('E:\Language_task_results_v3\Logical\run-1_sub2025082302_SW_MI_NullModelp25.nii');

vol4D_null=single(data_nii.img);
[N1,N2,N3,num_nullmodels]=size(vol4D_null);
vol4D_null=reshape(vol4D_null,[N1*N2*N3 num_nullmodels]);
Null_models=vol4D_null(indexBrain,:);
clear data_nii

% temp_file=fullfile([input_dir '\' sub],'**\*SW_MI.nii');
% temp_file=dir( temp_file);
% data_nii=load_nii([temp_file.folder, '\',temp_file.name]);


% data_nii=load_nii('E:\HCP_retest_results\Emotion\results_sliding_window\neut_sub601127tSW_MI.nii');
data_nii=load_nii('E:\Language_task_results_v3\Logical\run-1_sub2025082302_SW_MI_logicalq25.nii');

vol4D_observed=single(data_nii.img(:,:,:,1:end));
Tasktimings=squeeze(vol4D_observed(1,1,1,:));

[N1,N2,N3,num_wid]=size(vol4D_observed);
vol4D_observed=reshape(vol4D_observed,[N1*N2*N3 num_wid]);
Observed_data=vol4D_observed(indexBrain,:);
% t=sum(Observed_data,1);
% 
% figure
% plot(t)
% 
bl=zeros(size(Observed_data,2),1);
for f=1:size(Observed_data,2)
    bl(f)=quantile(abs(Observed_data(:,f)),0.5);
end
figure
plot(bl)

left_area=zeros(length(indexBrain),num_wid);
parfor i=1:length(indexBrain)
%       [f,xi]=ksdensity(Null_models(i,:),'Function','cdf', 'Bandwidth',0.002,'NumPoints',num_nullmodels);
%       [f,xi]=ksdensity(Null_models(i,:),'Function','cdf', 'Bandwidth',0.005,'NumPoints',1*num_nullmodels);
null_data=Null_models(i,:);      
outliers=isoutlier(abs(null_data),'mean');
null_idx=find(outliers==0);
[f,xi]=ksdensity(Null_models(i,null_idx),'Function','cdf','Bandwidth',0.005,'NumPoints',1*500);


     [~,index]=min(abs(1*Observed_data(i,:)-xi(:)),[],1);
     left_area(i,:)=f(index);
end
% left_area(left_area<0.0000001)=0.0000001;
% left_area(left_area>0.9999999)=0.9999999;


Z_map=zeros(N1*N2*N3,num_wid);
P_map=1-left_area;

Z_map(indexBrain,:)=norminv(left_area);
Z_map=reshape(Z_map,N1,N2,N3,num_wid);
P_vector=sort(P_map(:),'ascend');
u=spm_uc_FDR(0.05,[1 Inf],'Z',1,P_vector,1);
Z_map(Z_map<u)=0;

% Z_map = removeSmallObjects3D(Z_map, 10);

data_nii.img=data_nii.img.*0;
data_nii.img(:,:,:,1:num_wid)=Z_map;
f_name=[outputpath,'\','run1_sub',sub,'_SW2sp25_Zmap'];
% save_nii(template_nii,[outputpath '\sub' f_name,'_original_98nei_sub_001_run_1.nii.gz'])
save_nii(data_nii,[f_name, '.nii']);

% % save_nii(template_nii,[f_name,'.nii']);
% result_dir=outputpath;
% str=['D:/workbench/bin_windows64/wb_command -volume-to-surface-mapping ',[f_name,'.nii '],'../Masks/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii ',[result_dir,'/BrainL.shape.gii'],' -trilinear'];
% eval(['!',str]);
% str=['D:/workbench/bin_windows64/wb_command -volume-to-surface-mapping ',[f_name,'.nii '],'../Masks/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii ',[result_dir,'/BrainR.shape.gii'],' -trilinear'];
% eval(['!',str]);
% str=['D:/workbench/bin_windows64/wb_command -cifti-create-dense-timeseries ',[f_name,'.dtseries.nii '],'-left-metric ',[result_dir,'/BrainL.shape.gii'],' -roi-left  ../Masks/L.atlasroi.32k_fs_LR.shape.gii -right-metric ',[result_dir,'/BrainR.shape.gii'],' -roi-right ../Masks/R.atlasroi.32k_fs_LR.shape.gii -timestep 1'];
% eval(['!',str]);
% delete([result_dir,'/BrainL.shape.gii']);
% delete([result_dir,'/BrainR.shape.gii']);



cZ_map=Z_map;
cZ_map(cZ_map>0)=1;

for k=1:length(Tasktimings)-1
    if Tasktimings(k)==1 && Tasktimings(k+1)==0
       Tasktimings(k)=0;Tasktimings(k-1)=0;
    end
end

taskind=find(Tasktimings==1);
cmap=squeeze(sum(cZ_map(:,:,:,taskind),4));
data_nii.img=(single(cmap/length(taskind)));
data_nii.hdr.dime.dim(5)=1;
f_name=[outputpath,'\','run1_sub',sub,'_SW2sv3p2_Cmap'];
% save_nii(template_nii,[outputpath '\sub' f_name,'_original_98nei_sub_001_run_1.nii.gz'])
save_nii(data_nii,[f_name, '.nii']);


end

119/182

