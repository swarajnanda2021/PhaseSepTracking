function def_cmap
%def_cmap Define nummeric camera mappings from the obtained distortion map
%   for fast image warping and interpolation.
%
%   Input:
%       Dmap.mat Distortion mapping 
%       Dpar.mat Distortion parameters
%       Kmat.mat Camera matrix 
%
%   Output:
%       wrp_c.mat Numeric camera warping data

%% get globals
global folder date cal prop ctrl Dmap

%% load parameters
Dpar=importdata([folder date cal vsl 'Dpar.mat']);
Kmat=importdata([folder date cal vsl 'Kmat.mat']);

%% solve

% define image grid
[y,x]=meshgrid(1:prop.res(2),1:prop.res(1));

% loop cameras
for c=1:prop.res(4)
    
    % define mapping
    map=matlabFunction(subs(sym(Dmap),sym('b',[1 length(Dpar{c}.map)]),Dpar{c}.map),'Optimize',false);
    
    % camera matrix as homography mapping
    H=Dpar{c}.H\Kmat{c};
    
    % evaluate full warping
    [X,Y,J,xx,yy,XX,YY]=imtrans(x,y,map,H,ctrl.optl);
    
    % write
    wrp.x=x;
    wrp.y=y;
    wrp.xx=xx;
    wrp.yy=yy;
    wrp.X=X;
    wrp.Y=Y;
    wrp.XX=XX;
    wrp.YY=YY;
    wrp.J=J;
    
    % save for each view
    save([folder date cal vsl 'wrp_',num2str(c),'.mat'],'wrp')
    
end

end