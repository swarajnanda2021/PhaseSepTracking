function [Y,val,Xbnds] = testmedian(X,Nsigma,Nseg)
%testmedian Segment and validate data based oncomponent wise median 
%   statistics 
%
%	Author: Abel-John Buchner / Koen
%
%   Input
%       X: Input data array, for example 3xN when 3D (e.g. velocity vectors)
%       Nsigma: Number of times sigma to set dynamic threshold
%       Nseg: Number of components to be invalid before segmentation [1 or 1xN]
%
%   Output
%       Y: Output data from median filter
%       val: Valid data flag which passed the median, logical array (1x1 | 1xN)
%       Xbnds: X bounds
%
%   Note: NaN values are ignored in case there

% Correct input
if nargin<2
    Nsigma  = 2; % discrete closest to 1.96 , can be any value
end
if nargin<3
    Nseg=1; % all have to be with median
end

%figure; quiver3(zeros(1,size(X,2)),zeros(1,size(X,2)),zeros(1,size(X,2)),X(1,:),X(2,:),X(3,:),0); axis equal

% Median
Xmed=nanmedian(X,2);

%hold on; quiver3(0,0,0,Xmed(1),Xmed(2),Xmed(3),0,'r'); axis equal

% Median deviation
Xsigma = nanmedian(abs(X-Xmed),2);
Xbnds=[Xmed-Nsigma*Xsigma Xmed+Nsigma*Xsigma ];

% Total segmentation components
val = ~( sum( abs(X - Xmed) > (Nsigma*Xsigma) ,1)>=Nseg ) ; %size(X,1)

% Segment data
Y = X(:,val);

%quiver3(zeros(1,size(Y,2)),zeros(1,size(Y,2)),zeros(1,size(Y,2)),Y(1,:),Y(2,:),Y(3,:),0); axis equal
%quiver3(0,0,0,Xmed(1),Xmed(2),Xmed(3),0,'r'); axis equal; hold off

end

