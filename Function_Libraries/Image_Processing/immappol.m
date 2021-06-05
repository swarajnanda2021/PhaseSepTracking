function [ lset ] = immappol( siz , pol )
%map shape(s!) by closed polygon coordinates UNTITLED Summary of this function goes here
%   Detailed explanation goes here: returns a lset that describes shape and
%   overlap
% slow by for looping singles :b1 = poly2mask(x1, y1, rows, columns)

% correct sie polygon
pol=reshape(pol,size(pol,2)*2,size(pol,3));
pol=pol';

% initiate msk
lset=zeros(siz);

% msk shapes by polygon contours
lset = sum(insertShape(lset','FilledPolygon',pol),3)';

end

