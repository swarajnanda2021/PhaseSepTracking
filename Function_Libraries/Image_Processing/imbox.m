function [ imb, ix , iy , int ] = imbox( imd , box , typ ) % ,std
%imbox Image box filter to change the resolution of images by integer
%   amount.
%   
%   Input
%       imd: image data
%       box: box filter
%       typ: direct, average, or integrated value
%
%   Box filtering allows to recude the image resolution by n and increase
%   the image resolution by 1/n.

% correct input
if nargin<3
    typ='int'; % integrate resolution
end

% make indexed ref grid
[ry,rx]=meshgrid(1:size(imd,2),1:size(imd,1));

% make indexed grid
[iy,ix]=meshgrid((ceil(box(2))+1)/2:box(2):size(imd,2)-(ceil(box(2))-1)/2,...
                 (ceil(box(1))+1)/2:box(1):size(imd,1)-(ceil(box(1))-1)/2);

% first do convolution (faster than for loop)
switch typ
    case 'dir'
        A=1;
        h=1;
    case 'avg'
        A=1;
        h=ones(ceil(box))/prod(ceil(box));
    case 'int' 
        A=prod(box);
        h=ones(ceil(box))/prod(ceil(box));
end
int=A*imfilter(imd,h,'symmetric'); % 'symmetric' only needed for int and not imb

% interpolate to lower resolution
imb=zeros(size(ix,1),size(ix,2),size(imd,3));
for i=1:size(imb,3)
    imb(:,:,i)=interp2(int(:,:,i)',ix',iy','cubic',0)';%floor(
end

end