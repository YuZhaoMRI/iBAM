function [brainMask3D,brainMaskDim3D,brainIndex]=getBrainMaskInfo(fileName)
%   INPUT:
%       fileName support (brain mask).nii file extension.
%
%   OUTPUT:
%       brainMask3D    
%       brainMaskDim3D  
%       brainIndex     
%
mask_nii=load_nii(fileName);
brainMask3D=double(mask_nii.img);
brainMaskDim3D=mask_nii.hdr.dime.dim(2:4);
brainMask1D=reshape(brainMask3D,[prod(brainMaskDim3D) 1]);
brainMask3D(brainMask3D>0)=1;
brainIndex=find(brainMask1D>0);