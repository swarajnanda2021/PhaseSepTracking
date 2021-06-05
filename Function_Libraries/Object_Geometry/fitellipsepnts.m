function [xm,ax,ang,c] = fitellipsepnts(x,typ)
%ellipsoidpoints Fit ellipsoid to 2D points, find the midpoint, use pca,
%   and use standard deviation time sqrt(3) to find uniform spread axis.
%   
%   Input
%       x spatio temporal points to fit to
%       typ type of ellipsoid axis ( free [best rotation] | fixed [no rotation] )
%       
%   Output
%       xm midpoint
%       ax axis ellipsoid
%       ang anglular position ellipsoid
%       quadric corresponding to the ellipsoid
%       
%   Note: 3D points use fitellipsoidpnts

% correct input
if nargin < 2
    typ='free';
end

% error message
if size(x,1)~=2
    error('Incorrect dimensional data fitellipsepnts.m')
end

%x=rand(2,2)*normrnd(0,2,2,10000)+50*rand(2,1); close all
%figure; plot(x(1,:),x(2,:),'.'); axis equal

% compute mean
xm=mean(x,2);

% center
x=x-xm;

%figure; plot(x(1,:),x(2,:),'.'); axis equal

% compute principle components using svd
switch typ
    case 'free'
        [~,~,R]=svd(x'); % find optimal rotation
        ang=rotm2eula([R zeros(2,1); zeros(1,2) 1],'ZXY')';
        ang=ang(1);
    case 'fixed'
        ang=0;
        R=eye(2);
end

% rotate
x=R\x;

%figure; plot(x(1,:),x(2,:),'.'); axis equal

% find rescaling
ax=mean(abs(x),2)+sqrt(3)*std(abs(x),[],2); % 

% rescale
x=x./ax;

%figure; plot(x(1,:),x(2,:),'.'); axis equal

% assemble quadric
if nargout ==4
    [ c ] = eshp2cvec( xm,ax,ang ); % 3*
end

%x=R*(ax.*x)+xm
%figure; plot(x(1,:),x(2,:),'.'); axis equal

%xe=cvec2pnts(contrans(c,3,'resize'),1000);
%hold on; plot(xe(1,:),xe(2,:),'.')

end

