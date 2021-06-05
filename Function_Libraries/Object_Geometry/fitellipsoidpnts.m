function [xm,ax,ang,q] = fitellipsoidpnts(x,typ)
%ellipsoidpoints Fit ellipsoid to 3D points, find the midpoint, use pca,
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
%   Note: 2D points use ellipsepoints

% correct input
if nargin < 2
    typ='free';
end

% error message
if size(x,1)~=3
    error('Incorrect dimensional data fitellipsoidpnts.m')
end

%x=rand(3,3)*normrnd(0,2,3,10000)+50*rand(3,1); close all
%figure; plot3(x(1,:),x(2,:),x(3,:),'.'); axis equal

% compute mean
xm=mean(x,2);

% center
x=x-xm;

%figure; plot3(x(1,:),x(2,:),x(3,:),'.'); axis equal

% compute principle components using svd
switch typ
    case 'free'
        [~,~,R]=svd(x'); % find optimal rotation
        ang=rotm2eula(R,'ZXY')';
    case 'fixed'
        ang=zeros(3,1);
        R=eye(3);
end

% rotate
x=R\x;

%figure; plot3(x(1,:),x(2,:),x(3,:),'.'); axis equal

% find rescaling
ax=mean(abs(x),2)+sqrt(3)*std(abs(x),[],2); % 

% rescale
x=x./ax;

%figure; plot3(x(1,:),x(2,:),x(3,:),'.'); axis equal

% assemble quadric
if nargout ==4
    [ q ] = eshp2qvec( xm,ax,ang ); % 3*
end

%x=R*(ax.*x)+xm
%figure; plot3(x(1,:),x(2,:),x(3,:),'.'); axis equal

%xe=qvec2pnts(quadtrans(q,3,'resize'),10000);
%hold on; plot3(xe(1,:),xe(2,:),xe(3,:),'.')

end

