clear;clc;

outputpath='F:\Language_task_results_v5\Logical';
if ~isfolder(outputpath)
    mkdir(outputpath)
end

%% get the participant folders
inputpath = 'F:\Language_task';
% Get all contents in the specified path (both files and folders)
contents = dir(inputpath);
% Filter out folders (excluding '.' and '..')
sublist = contents([contents.isdir]);  % Get all folders
sublist = sublist(~ismember({sublist.name}, {'.', '..'}));  % Exclude current and parent directories


%% Set flags to determine which experimental runs will be processed and analyzed
run_label='run-1'; % run-1: logical; run-2: Semantic​​; run-3: Visual
if run_label=='run-1' 
   Events_file='logical_grammer.txt';      nii_name='Language_filtered_func_data_run-1_clean_noise_second_tr.nii.gz';  
end
if run_label=='run-3'
   Events_file='geometrical.txt';          nii_name='Language_filtered_func_data_run-3_clean_noise_second.nii.gz';  
end


%% Set Parameters
numNullModels=500;
freqBand = [0.01 0.18];  % Frequency band for temporal filtering (Hz)
FWHM = 4;                % Full width at half maximum for spatial smoothing (mm)
dumypoints = 5;          % Number of initial time points to remove
durRest = 0 * 60;        % Duration of the resting-state fMRI scan (seconds)
durHRF = 15;             % Duration of the canonical HRF (seconds)
scaleFactor=10;



%% load brain mask
maskFileName='../Masks/MNI152_T1_2mm_brain.nii';
[brainMask3D,brainMaskDim3D,brainIndex]=getBrainMaskInfo(maskFileName);

%% load the matrix that ecode the functional neighborhoods
load('../Masks/nei_ind_and_w_matrix_sum_cosine_sim_funcmat_819_maxk_72_cenw_1_MNI152_T1_2mm_brain.mat');

delete(gcp('nocreate'));
parpool(30);
%% iBAM for all subjects
display(['Number of subjects is ',num2str(length(sublist))]);
for m=1:length(sublist)

display(['iBAM for subject #',num2str(m)])
sub=char(sublist(m).name);
outputpath_sub=[outputpath,'\',sub];
if ~isfolder(outputpath_sub)
    mkdir(outputpath_sub)
end
% load fMRI volumes 

%% Set flags to determine which experimental runs will be processed and analyzed
nii_name='Language_filtered_func_data_run-3_clean_noise_second.nii.gz';  
nii_file=fullfile(sublist(m).folder, sublist(m).name, 'run-3', nii_name);
dataNii=load_nii(nii_file);
TR=dataNii.original.hdr.dime.pixdim(5);
NumRestFrame=ceil(durRest/TR);
fMRIData4D=double(dataNii.img(:,:,:,end-600:end));% remove the dummy and resting-state fMRI scans at the end 
[N1,N2,N3,N4]=size(fMRIData4D);
voxelsize=dataNii.original.hdr.dime.pixdim(2);
%% spatial smoothinng 
for ii=1:N4
    fMRIData4D(:,:,:,ii)=imgaussfilt3(fMRIData4D(:,:,:,ii).*brainMask3D,FWHM/voxelsize/2.355);  %三维图像的三维高斯滤波
end

%% preprocessing of fMRI time courses, including (1) detrend (2) bandpass filter (3)normalize
fMRIData2D=reshape(fMRIData4D,[N1*N2*N3 N4]);
fMRIData2D=fMRIData2D(:,:)';
fMRIData2D=detrend(fMRIData2D,2);

%% performs global signal regression
grayMatterMaskFile='../Masks/GrayMatterMask.nii';
fMRIData2D=regressGlobalSignal(fMRIData2D,grayMatterMaskFile);


%% The combined temporal filtering strategy effectively removes noise without introducing phase delays into the BOLD signal.
fMRIData2D = applyBandpassFilter(fMRIData2D, TR);

%% normalize the fMRI signals
fMRIData2D=zscoreFMRI(fMRIData2D);



%% remove the baselines of fMRI data
[fMRIData2D, baseLineInfo] = baselineCorrect(fMRIData2D);


%% 取滑动窗口的起始位置
slidingWindow=1:ceil(durHRF/TR); 
numWindFrame=size(fMRIData2D,1)-length(slidingWindow)+1;  %the number of the frames that can be calculated
randWindFrame=randsample(1:numWindFrame,numNullModels);
FrameCounter=0;
ModeIntensity2D=[];
for n=randWindFrame 
FrameCounter=FrameCounter+1;
display(['Frame #',num2str(FrameCounter)])
index_window=slidingWindow+n-1;
dataInWindow=fMRIData2D(index_window,:);
MI=zeros(length(brainIndex),1);

parfor ind=1:length(brainIndex)
    neighbIndex = nei_ind_sum(:,ind); % load from nei_ind_and_w_matrix_sum_cosine_sim_funcmat_819_maxk_98_cenw_1_MNI152_T1_2mm_brain.mat
    neigbWeight = nei_w_sum(:,ind);   % load from nei_ind_and_w_matrix_sum_cosine_sim_funcmat_819_maxk_98_cenw_1_MNI152_T1_2mm_brain.mat
    if isempty(find(neighbIndex))
       continue; 
    end 
    [~,centerPosition]=maxk(neigbWeight,1);
    X=dataInWindow(:,neighbIndex);
    W=repmat(neigbWeight',size(X,1),1); 
    X=X.*W;
    [U,S,V]=svds(double(X),1,'largest');
    MI(ind)=S.*sign(sum(sign(V(centerPosition)).*U))/sqrt(size(X,1))/scaleFactor;

end
ModeIntensity2D=[ModeIntensity2D,MI];
end
% save([sub,'_MI_wind.mat'],'MI_2D','-v7.3');

MI_save=zeros(prod(brainMaskDim3D),numNullModels);%
MI_save(brainIndex,:)=ModeIntensity2D;
MI_save=reshape(MI_save,[brainMaskDim3D,numNullModels]);




% 时间序列的转换，注意命令和后缀
templateNii=dataNii;
templateNii.img=templateNii.img.*0;
% template_nii.img(:,:,:,1+dumypoints:wind_frame+dumypoints)=single(MI_save);
templateNii.img=(single(MI_save));
templateNii.hdr.dime.dim(5)=numNullModels;
templateNii.original.hdr.dime.dim(5)=numNullModels;

f_name=[outputpath_sub,'\',run_label,'_sub',sub,'_SW_MI_RestReference'];
% save_nii(template_nii,[outputpath '\sub' f_name,'_original_98nei_sub_001_run_1.nii.gz'])
save_nii(templateNii,[f_name, '.nii'])


% save_nii(template_nii,[f_name,'.nii']);
result_dir=outputpath_sub;
str=['D:/workbench/bin_windows64/wb_command -volume-to-surface-mapping ',[f_name,'.nii '],'../Masks/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii ',[result_dir,'/BrainL.shape.gii'],' -trilinear'];
eval(['!',str]);
str=['D:/workbench/bin_windows64/wb_command -volume-to-surface-mapping ',[f_name,'.nii '],'../Masks/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii ',[result_dir,'/BrainR.shape.gii'],' -trilinear'];
eval(['!',str]);
str=['D:/workbench/bin_windows64/wb_command -cifti-create-dense-timeseries ',[f_name,'.dtseries.nii '],'-left-metric ',[result_dir,'/BrainL.shape.gii'],' -roi-left  ../Masks/L.atlasroi.32k_fs_LR.shape.gii -right-metric ',[result_dir,'/BrainR.shape.gii'],' -roi-right ../Masks/R.atlasroi.32k_fs_LR.shape.gii -timestep 1'];
eval(['!',str]);
delete([result_dir,'/BrainL.shape.gii']);
delete([result_dir,'/BrainR.shape.gii']);

end
delete(gcp('nocreate'));