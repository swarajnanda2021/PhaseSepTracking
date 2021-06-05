function [ xx,yy,XX,YY] = boxgrid( x,y,X,Y ,box)
%boxgrid box a image meshgrid including its dewarped grid accompanying the
%   imbox function

% make indexed ref grid
[ry,rx]=meshgrid(1:size(y,2),1:size(x,1));

% make indexed box grid
[iy,ix]=meshgrid((ceil(box(2))+1)/2:box(2):size(y,2)-(ceil(box(2))-1)/2,...
                 (ceil(box(1))+1)/2:box(1):size(x,1)-(ceil(box(1))-1)/2);

% interpolate rectilinear grid
if ~isempty(x) && ~isempty(y)
    xx=interp2(ry,rx,x,iy,ix,'linear',0);
    yy=interp2(ry,rx,y,iy,ix,'linear',0);
else
    xx=[];
    yy=[];
end

% interpolate deformed grid
if ~isempty(x) && ~isempty(y) && ~isempty(X) && ~isempty(Y)
    XX=interp2(y,x,X,yy,xx,'cubic',0);
    YY=interp2(y,x,Y,yy,xx,'cubic',0);
else
    XX=[];
    YY=[];
end

end

