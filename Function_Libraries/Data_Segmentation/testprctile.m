function [Y,val,IPR] = testprctile(X,IPRlim,Nseg)
%testpercentile Segment and validate data based oncomponent wise quantile 
%   statistics 
%
%	Author: Abel-John Buchner / Koen
%
%   Input
%       X: Input data array, for example 3xN when 3D (e.g. velocity vectors)
%       IQRlim: limits inter quantile range
%       Nseg: Number of components to be invalid before segmentation [1x1 | 1xN]
%
%   Output
%       Y: Output data from median filter
%       val: Valid data flag which passed the median, logical array (1xN)
%       IPR: Inter Percentile Range
%
%   Note: NaN functionality not tested..

% Correct input
if nargin<2
    IPRlim= [25 75]; % default values
end
if nargin<3
    Nseg=1; % all have to be with median
end

%figure; quiver3(zeros(1,size(X,2)),zeros(1,size(X,2)),zeros(1,size(X,2)),X(1,:),X(2,:),X(3,:),0); axis equal

% Interquantile range
IPR = prctile(X,IPRlim,2); % this call is SLOW in a loop.... so use IQRlims=[0 1] where possible

% Total segmentation components
val = ~( sum( X<IPR(:,1) | X>IPR(:,2) ,1)>=Nseg ) ; %size(X,1)

% Segment data
Y = X(:,val);

%hold on; quiver3(zeros(1,size(Y,2)),zeros(1,size(Y,2)),zeros(1,size(Y,2)),Y(1,:),Y(2,:),Y(3,:),0); axis equal

end

