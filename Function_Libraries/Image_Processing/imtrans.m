function [X,Y,J,xx,yy,XX,YY] = imtrans( x , y , map , H , tol)
%evaldmap Evaluate distortion mapping for image manipulation
%
%   Input
%       xy grid
%       distortion mapping
%       homography mapping
%
%   Output1
%       X Y     deformed from reference config
%       J       jacobian for image correction in reference config
%       XX YY   square resample grid in def config
%       xx yy   deformed to reference config
%           computetional intesive due optimization

%-- correct input
if nargin<4
    H=eye(3);
end
if nargin<5
    tol=10^-3;
end

%-- Part 1: Forward displacement

%-- 0 initiate coordinates to manipulate
X=[x(:)'
    y(:)'];

%-- 1 Distortion correction
if ~isempty(map)
    
    % distortion dewarp coordinates
    X=map(X(1,:),X(2,:));
    
end

%-- 2 linear image transform after distortion correction

% camera dewarp coordinates
X=homc2inhc(H\inhc2homc(X));

% meshgrid
Y=reshape(X(2,:),size(y));
X=reshape(X(1,:),size(x));

%-- 3 Intensity correction
if nargout>2
    
    % fit local poly-surface
    Xcoef=SGfilt_coef(X,strel('disk',1,0),[1 1 0],1); % flat, disk, same as central scheme / sobel filter
    Ycoef=SGfilt_coef(Y,strel('disk',1,0),[1 1 0],1); 
    
    % compute derivatives
    dXdx=SGfilt_eval(Xcoef,[1 0 0]);
    dXdy=SGfilt_eval(Xcoef,[0 1 0]);
    dYdx=SGfilt_eval(Ycoef,[1 0 0]);
    dYdy=SGfilt_eval(Ycoef,[0 1 0]);
    
    % volume deformation
    J=dXdx(:).*dYdy(:)-dXdy(:).*dYdx(:);
    
    % mesh to grid
    J=reshape(J,size(x))./mean(J(:)); % regularize by mean value to avoid small and large numbers
    
else
    
    J=[];
    
end

%-- Part 2: nummeric inversion mapping
if nargout>3
    
    %-- 4 New square grid
    limX=[min(X(:)) max(X(:))];
    limY=[min(Y(:)) max(Y(:))];
    dX=min(range(X(:))/size(x,1),range(Y(:))/size(y,2)); % interpolate at highest resolution
    [YY,XX]=meshgrid(limY(1):dX:limY(2),limX(1):dX:limX(2));
    
    %-- 6 invert linear image transform first
    
    % linear inversion coordinates
    xx=homc2inhc(H*inhc2homc([XX(:)'
                              YY(:)']));
    
    %-- 5 Optimization problem for distortion inversion
    if ~isempty(map)
        
        % objective
        ObjF=matlabFunction(defobjf(sym('X',[2,1])-sym(map),[]),'Optimize',false);
        
        % variabes
        vars=sort(symvar(sym(ObjF)));
        
        % numbering knowns
        nk=(1:2)';
        
        % numbering unknowns
        nu=(1:length(vars))';
        nu(nk)=[];
        
        % find grad and hess
        [gObjF,hObjF]=difobjf(ObjF,vars(nu));
        gObjF=matlabFunction(gObjF,'Optimize',false);
        hObjF=matlabFunction(hObjF,'Optimize',false);
        
        %-- 6 nummeric inversion
        
        % initiate square solution grid undistorted configuration
        limx0=[min(x(:)) max(x(:))];
        limy0=[min(y(:)) max(y(:))];
        [yy0,xx0]=meshgrid(linspace(limy0(1),limy0(2),size(YY,2)),...
            linspace(limx0(1),limx0(2),size(XX,1)));
        
        % reshape
        xx0=[xx0(:)'
            yy0(:)'];
        
        % initiate solution write
        xxsol=zeros(size(xx0));
        
        % loop point seperate for stable local convergence
        for k=1:size(xx,2)
            
            % numering assembly
            asn=reshape((1:numel(xx(:,k)))',2,[]); % assembly numbering
            
            % get solution vec
            x0=reshape(xx0(:,k),[],1); % xx
            
            % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
            dat=[zeros(size(nk)) nk xx(:,k)
                ones(size(nu)) nu asn];
            
            % solve global optimization (exact solution)
            [xsol,res] = optobjf( [], ObjF, gObjF, hObjF , x0 , dat , tol ,'newton','no');%1/tol,10,1);
            
            % write
            xxsol(:,k)=reshape(xsol,2,[]);
            
            % message
            disp(['inverted image point ',num2str(k),' out of ',num2str(size(xx,2)),' at disparity ',num2str(res)])
            
        end
        
        % [over]write initial inversion
        xx=reshape(xxsol,2,[]);
        
    end
    
    % mesh to grid
    yy=reshape(xx(2,:),size(YY));
    xx=reshape(xx(1,:),size(XX));
    
end

end