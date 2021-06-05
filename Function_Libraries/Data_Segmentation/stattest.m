function [dat,val,bnds] = stattest(dat,outl,conf,iter)
%stattest Statistical test on outlier data for segmentation purpose
%
%   Input
%       dat Some data
%       outl Outlier test
%       conf Confidence parameters specific test
%       iter Max number iteration (0-inf)
%
%   Output
%       Segmented data
%       Segmentation flag
%       Statistical bounds
%

% correct input
if nargin<2
    outl='meanstd'; % simple mean std model
end
if nargin<3
    conf=3; % safe confidence
end
if nargin<4
    iter=1; % no iteration
end

% iterate outlier removal
cnt=1; % initiate iteration count
siz=0; % initiate size valid
val=true(1,size(dat,2)); % initiate valid set
while cnt<=iter && siz~=nnz(val) % if below iteration and if remaining to change
    
    % check size valid
    siz=nnz(val);
    
    % select test
    switch outl
        case 'median'
            [~,V,bnds] = testmedian(dat(:,val),conf);
        case 'meanstd'
            [~,V,bnds] = testmeanstd(dat(:,val),conf);
        case 'quantile'
            [~,V,bnds] = testquantile(dat(:,val),conf);%[.10 .90]
        case 'prntile'
            [~,V,bnds] = testprctile(dat(:,val),conf);%[10 90]
        case 'none'
            V=val;
            bnds=[-inf inf];
        otherwise
            error('stattest.m: Not listed outlier test')
    end % test
    
    % update valid set
    val(val)=V;
    
    % count iterations
    cnt=cnt+1;
    
end % cnt

% segment data
dat=dat(:,val);

end

