function [fgr,bgr,segp,segm] = imrelscl( imd , fobj, iter )
%image relavant scale by imposed filter UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% correct input
if nargin<3
    iter=0;
end

% size filter and range
switch class(fobj)
    case 'double'
        
        siz=fobj(1:2);
        sig=sqrt(prod(siz))/2/sqrt(2)*2/3;
        
    case 'strel'
        siz=[size(fobj.Neighborhood,1) size(fobj.Neighborhood,2)];
        sig=sqrt(prod(siz))/2/sqrt(2)*2/3;
end

% filters from decomposed mexican hat 
G=fspecial('gaussian',3*siz,sig);%3*(siz-1)+1
LoG=-1/2*fspecial('log',3*siz,sig)*sig^2; % scaled to same intensity
B=G-LoG; % to retrieve the background filter 

% integral value, this is import when maximizing the contrast, stablelizes the object identification!
scl=max(G(:));
G=G./scl;
B=B./scl;

%figure; surf(G); view(2)
%figure; surf(LoG); view(2)
%figure; surf(B); view(2)

% define objects to interest
fgr=imfilter(imd,G,'symmetric');%./imfilter(imd,G,'symmetric');%medfilt3((scl-bgr)./bgr,[1 1 3]); % object[ive] scale
bgr=imfilter(imd,B,'symmetric');%./imfilter(imd,G,'symmetric');%medfilt3((scl-bgr)./bgr,[1 1 3]); % object[ive] scale

% define noise level
seg=imd>bgr;
devp=3*std(imd(seg)-fgr(seg));
devm=3*std(imd(~seg)-fgr(~seg));

for i=1:iter % fixed s.t. no blow up time
    % segment before identification
    segp=(imd-bgr)>devp;
    segm=(imd-bgr)<-devm;
    seg=segm|segp;%imdilate(,ones(3));
    
    % explicit removal data from image
    for k=1:size(bgr,3)
        bgr(:,:,k)=regionfill(imd(:,:,k),seg(:,:,k));
    end
    
    bgr=imfilter(bgr,B,'symmetric'); % filter noise to scale
    
end 

% segment before identification
segp=(imd-bgr)>devp;
segm=(imd-bgr)<-devm; % imdilate(,fobj);%;
    
end
