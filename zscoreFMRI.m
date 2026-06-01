function fMRIData2D=zscoreFMRI(fMRIData2D)
%ZSCOREFMRI  Z-score normalizes fMRI voxel time series (zero mean, unit variance).
%
%   dataNorm = ZSCOREFMRI(data) standardizes each voxel's BOLD time series
%   to have zero mean and unit standard deviation (z-scoring). 
%
%   Input:
%       data      - fMRI data, - 2D matrix [time x voxels]
%
%   Output:
%       dataNorm  - Z-scored data, same size and shape as input.
%
%   Notes:
%       - Voxels with zero standard deviation are assigned std = 1 to avoid
%         division by zero (typically constant background voxels).
%       - Operates along the time dimension.
%
%   Example:
%       dataZ = zscoreFMRI(fmriData4D);
    mean_voi=mean(fMRIData2D);
    fMRIData2D=fMRIData2D-repmat(mean_voi,size(fMRIData2D,1),1);
    std_voi=std(fMRIData2D,1);
    std_voi(std_voi==0)=1;
    fMRIData2D=fMRIData2D./repmat(std_voi,size(fMRIData2D,1),1);
end