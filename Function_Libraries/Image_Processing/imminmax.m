function corrected = imminmax(original,fobj,trs)
%improc_minmax Image min-max filter based on pwMinMax by J. Westerweel
%
% Input
%   Orginal Image
%   Windowed Filter Size (Circular)
%   Level Filter Tresholds
%
% Output
%   Corrected Image
%
% Edits by Koen Muller to pwMinMax
% 14-7-2016 symm. boundaries
% 17-09-2016 disk neighborhood
% 01-05-2017 relative segmentation and included mask

%size filter object
siz=size(fobj.Neighborhood);

% find range
rng=range(reshape(original,[],1));
lvl=trs*rng;

% Find min and max
low = imerode(original,fobj); % equiv. local min
upp = imdilate(original,fobj); % equiv. local max

% Smooth min and max
low = imfilter(low,double(fobj.Neighborhood),'symmetric')/sum(reshape(fobj.Neighborhood,[],1));
upp = imfilter(upp,double(fobj.Neighborhood),'symmetric')/sum(reshape(fobj.Neighborhood,[],1));

% Find contrast
contrast = upp-low;

% Put lower limit on difference between upp and low
cont = immultiply(contrast > lvl ,contrast-lvl) + lvl;

% Filtered image
corrected = imdivide(original-low,cont);

% mask half boundary points
% corrected(1:ceil(siz(1)/2),:,:)=0;
% corrected(:,1:ceil(siz(2)/2),:)=0;
% corrected(end-floor(siz(1)/2):end,:,:)=0;
% corrected(:,end-floor(siz(2)/2):end,:)=0;

end