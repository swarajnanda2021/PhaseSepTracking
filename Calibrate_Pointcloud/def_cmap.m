function def_cmap
%eval_cmap get numeric camera mapping for using them

% globals
global folder date cal prop ctrl Dmap

% load
Dpar=importdata([folder date cal vsl 'Dpar.mat']);
Kmat=importdata([folder date cal vsl 'Kmat.mat']);

% solve
for c=1:prop.res(4)
    
%     % image data
%     imd=log(double(import_frames({folder date [cal vsl 'Calibration' vsl 'camera' num2str(c)]},prop.ext,1,1))+1); % board scatters
    
    % local grid
    [y,x]=meshgrid(1:prop.res(2),1:prop.res(1));%size(imd,2),1:size(imd,1)); % meshgrid image
    
    % define mapping
    map=matlabFunction(subs(sym(Dmap),sym('b',[1 length(Dpar{c})]),Dpar{c}),'Optimize',false);
    
    % define matrix
    H=Kmat{c};
    
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
    
    save([folder date cal vsl 'wrp_' num2str(c) '.mat'],'wrp')
end

end

