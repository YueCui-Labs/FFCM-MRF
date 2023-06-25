function L=LUT2label(im,LUT)
% Create a label image using a look-up-table (LUT) obtained from clustering
% of image intensities using either 'FastFCMeans' or 'FastCMeans' functions.  
%
% INPUT:
%   - im  : N-dimensional grayscale image in integer format. 
%   - LUT : L-by-1 array that specifies the intensity-class relations,
%           where L is the dynamic intensity range of the input image. 
%           Specifically, LUT(1) corresponds to the label assigned to 
%           min(im(:)) and LUT(L) corresponds to the label assigned 
%           to max(im(:)). 
%
% OUTPUT:
%   - L   : label image of the same size as the input image. 
%


% Check image format
if isempty(strfind(class(im),'int'))
    error('Input image must be specified in integer format (e.g. uint8, int16)')
end

% Intensity range
Imin=min(im(:));
Imax=max(im(:));
I=Imin:Imax;

% Create label image
L=zeros(size(im),'uint8');
idx=unique(LUT);
for k=1:numel(idx)
    
    % Intensity range for k-th class
    i=find(LUT==idx(k));
    i1=i(1);
    if numel(i)>1
        i2=i(end);
    else
        i2=i1;
    end
    
    % Map the intensities in the range [I(i1),I(i2)] to class k 
    bw=im>=I(i1) & im<=I(i2);
    L(bw)=idx(k);
    
end

