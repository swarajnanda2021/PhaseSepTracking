function [ img ] = imevallset( x , y , fun , cdat )
%reconstruct image from list of levelsets, some notes for develop
%   images are positive number so are their levelsets
%   images are 2d thus loop outside here for moving levelsets
%   images are additive so are the levelsets
%
%   input: image grid, levelset function, coefficient data

% initiate image
img=zeros(size(x));

% loop mapping the levelsets
for i=1:size(cdat,2)
    
    imn=fun(cdat(:,i)',x,y);
    
    imn=imn.*(imn>0);
    
    img=img+imn;
    
end

end