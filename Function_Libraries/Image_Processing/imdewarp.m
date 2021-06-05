function [ dew ] = imdewarp( x , y , imd , xx , yy , int , extrap )
%imdewarp Dewarp a distorted image to a rectified one on grid
%   Inputs
%       - x y   reference configuration
%       - I     image intesity values
%       - xx yy interpolation in reference configuration
%       - int   interpolation method when interpolating to grid
%           interpolation asways return a square grid for image processing

%-- 0 Correct input

% no interpolate to grid?
if nargin <6
    
    int='linear'; % make empty, thus dont interpolate
    
end
if nargin < 7
    extrap=nan;
end
if strcmp(extrap,'bval')
    % averaging kernel
    hx=ones(ceil(size(imd,1)/32),3,size(imd,3));
    hx=hx/sum(hx(:));
    hy=ones(3,ceil(size(imd,2)/32),size(imd,3));
    hy=hy/sum(hy(:));
    
    % average image boundary
    imbx=imd;
    imbx(:,1:3,:)=imfilter(imbx(:,1:3,:),hx,'symmetric');
    imbx(:,end-2:end,:)=imfilter(imbx(:,end-2:end,:),hx,'symmetric');
    imbx(1:3,:,:)=imfilter(imbx(1:3,:,:),hy,'symmetric');
    imbx(end-2:end,:,:)=imfilter(imbx(end-2:end,:,:),hy,'symmetric');
    
    imby=imd;
    imby(1:3,:,:)=imfilter(imby(1:3,:,:),hy,'symmetric');
    imby(end-2:end,:,:)=imfilter(imby(end-2:end,:,:),hy,'symmetric');
    imby(:,1:3,:)=imfilter(imby(:,1:3,:),hx,'symmetric');
    imby(:,end-2:end,:)=imfilter(imby(:,end-2:end,:),hx,'symmetric');
    
    imd=(imbx+imby)/2;
    
    % lower boundaries image
    bpad=[ceil((min(x(:))-min(xx(:)))/(x(2,1)-x(1,1)))*(min(xx(:))<min(x(:)))...
        ceil((min(y(:))-min(yy(:)))/(y(1,2)-y(1,1)))*(min(yy(:))<min(y(:)))];
    imd = padarray(imd,bpad+1,'replicate','pre'); % supress image data edges
    x=2*padarray(x,bpad+1,'replicate','pre')...
        -padarray(x,bpad+1,'symmetric','pre');
    y=2*padarray(y,bpad+1,'replicate','pre')...
        -padarray(y,bpad+1,'symmetric','pre');
    
    % remove dubplicate padding
    imd(bpad(1)+1,:,:)=[]; imd(:,bpad(2)+1,:)=[];
    x(bpad(1)+1,:)=[]; x(:,bpad(2)+1)=[];
    y(bpad(1)+1,:)=[]; y(:,bpad(2)+1)=[];
    
    % upper boundaries image
    bpad=[ceil((max(xx(:))-max(x(:)))/(x(2,1)-x(1,1)))*(max(xx(:))>max(x(:)))...
        ceil((max(yy(:))-max(y(:)))/(y(1,2)-y(1,1)))*(max(yy(:))>max(y(:)))];
    imd = padarray(imd,bpad+1,'replicate','post'); % supress image data edges
    x=2*padarray(x,bpad+1,'replicate','post')...
        -padarray(x,bpad+1,'symmetric','post');
    y=2*padarray(y,bpad+1,'replicate','post')...
        -padarray(y,bpad+1,'symmetric','post');
    
    % remove dublicate padding
    imd(end-bpad(1),:,:)=[]; imd(:,end-bpad(2),:)=[];
    x(end-bpad(1),:)=[]; x(:,end-bpad(2))=[];
    y(end-bpad(1),:)=[]; y(:,end-bpad(2))=[];
    
    % overwrite with standard extrapval
    extrap=0;
end

%-- 1 Interpolate to grid
dew=zeros(size(xx,1),size(xx,2),size(imd,3));
for i=1:size(dew,3)
    dew(:,:,i)=interp2(y,x,imd(:,:,i),yy,xx,int,extrap);
end

end