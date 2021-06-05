function [Y,val,Xbnds] = testmeanstd(X,Nsigma,Nseg)
%testmeansts Segment and validate data based oncomponent wise mean and 
%   standard deviation statistics 
%
%	Author: Abel-John Buchner / Koen
%
%   Input
%       X: Input data array, for example 3xN when 3D (e.g. velocity vectors)
%       Nsigma: Number of times sigma to set dynamic threshold
%       Nseg: Number of components to be invalid before segmentation [1x1 | 1xN]
%
%   Output
%       Y: Output data from median filter
%       val: Valid data flag which passed the median, logical array (1xN)
%       Xbnds: X bounds
%
%   Note: NaN values are ignored in case there

% Correct input
if nargin<2
    Nsigma  = 1.96; % 95[%] confidence range
end
if nargin<3
    Nseg=1; % all have to be with median
end

%figure; quiver3(zeros(1,size(X,2)),zeros(1,size(X,2)),zeros(1,size(X,2)),X(1,:),X(2,:),X(3,:),0); axis equal

% Mean
Xmean=nanmean(X,2);

%hold on; quiver3(0,0,0,Xmean(1),Xmean(2),Xmean(3),0,'r'); axis equal

% Standard deviation
Xsigma = nanstd(abs(X-Xmean),[],2);
Xbnds=[Xmean-Nsigma*Xsigma Xmean+Nsigma*Xsigma ]; % bounds

% Total segmentation components
val = ~( sum( abs(X - Xmean) > (Nsigma*Xsigma) ,1)>=Nseg ) ; %size(X,1)

% Segment data
Y = X(:,val);

%quiver3(zeros(1,size(Y,2)),zeros(1,size(Y,2)),zeros(1,size(Y,2)),Y(1,:),Y(2,:),Y(3,:),0); axis equal
%quiver3(0,0,0,Xmean(1),Xmean(2),Xmean(3),0,'r'); axis equal; hold off

end

