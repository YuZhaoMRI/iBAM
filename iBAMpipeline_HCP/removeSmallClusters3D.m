function img_clean_T = removeSmallClusters3D(img, minSize, varargin)
%REMOVESMALLOBJECTS3D  Remove 3D connected components smaller than minSize voxels
%   Works on any MATLAB version (R2014b+), supports 6/18/26 connectivity
%   Preserves original intensity values of large objects
%
%   img_clean = removeSmallObjects3D(img, minSize)
%   img_clean = removeSmallObjects3D(img, minSize, 'Connectivity', 26)
%
%   Example:
%       cleaned = removeSmallObjects3D(label3D, 10);
%       cleaned = removeSmallObjects3D(seg, 20, 'Connectivity', 6);

%% Input parsing
p = inputParser;
addRequired(p,'img');
addRequired(p,'minSize',@(x) isscalar(x) && x>0);
addParameter(p,'Connectivity',6,@(x) ismember(x,[6 18 26]));
parse(p,img,minSize,varargin{:});

conn = p.Results.Connectivity;

%% Method that works in ALL MATLAB versions (including very old ones)
[N1,N2,N3,N4]=size(img);
img_clean_T=zeros(size(img));
for j=1:N4
img3D=squeeze(img(:,:,:,j));
BW = img3D > 0;                                      % binary mask

CC = bwconncomp(BW, conn);                         % 3D connected components (supports 3D from very early versions)

stats = regionprops(CC, 'Area', 'PixelIdxList');   % get size of each component

% Find components that are large enough
largeComponents = find([stats.Area] >= minSize);

% Create output image
img_clean = zeros([N1,N2,N3], 'like', img);         % preserve original data type
for i = 1:numel(largeComponents)
    img_clean(CC.PixelIdxList{largeComponents(i)}) = img3D(CC.PixelIdxList{largeComponents(i)});
end
    img_clean_T(:,:,:,j)=img_clean;
end

end