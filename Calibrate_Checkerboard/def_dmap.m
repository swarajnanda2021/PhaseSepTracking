function [ Dmap ] = def_dmap
%def_dmap Define distortion mapping of choose
%   
%   Input: processing setting
%
%	Output:
%       Dmap.mat distortion mapping as symbolic matlab function

%% get globals
global folder date cal ctrl

%% Define mapping
% image coordinates
x=sym('x',[2 1]);

% available distortion mappings
switch ctrl.dmap
    case 'scalar-polymap' % mapping from scalar function
        
        % order poly-coefficient
        [oy,ox]=meshgrid(0:ctrl.mord+1,0:ctrl.mord+1);
        T=(ox+oy)<=ctrl.mord+1 & (ox+oy)>2; % remove linear terms
        ox=ox(T);
        oy=oy(T);
        T=reshape(T,[],1);
        
        % coefficients b for selected terms
        b=sym('b',[1 nnz(T)]);
        
        % define mapping
        f=b*((x(1).^ox).*(x(2).^oy)); % i.e. first + h.o.t.
        Dmap=x+gradient(f,x); % still arbitrary for translation etc.
        
    case 'vector-polymap' % vector mapping
        
        % order poly-coefficient
        [oy,ox]=meshgrid(0:ctrl.mord,0:ctrl.mord);
        T=(ox+oy)<=ctrl.mord & (ox+oy)>1;
        ox=ox(T);
        oy=oy(T);
        T=reshape(T,[],1);
        
        % coefficients b for selected terms
        b=sym('b',[1 2*nnz(T)]);
        
        % define mapping
        Dmap=x+reshape(b,2,[])*((x(1).^ox).*(x(2).^oy)); % normalized affine transform
        
    case 'radial-polymap' % polynomial radial distortion
        
        % order poly-coefficient
        [oy,ox]=meshgrid(0:ctrl.mord+1,0:ctrl.mord+1);
        T=(ox+oy)<=ctrl.mord+1 & (ox+oy)>2; % remove linear terms
        T=flipud(diag([repmat([0 1],1,(ctrl.mord+1)/2),0])==0) &T; % remove unnecesary terms
        ox=ox(T);
        oy=oy(T);
        T=reshape(T,[],1);
        
        % coefficients b for selected terms
        b=sym('b',[1 nnz(T)]);
        
        % define mapping
        f=b*((x(1).^ox).*(x(2).^oy)); % i.e. first + h.o.t.
        Dmap=x+gradient(f,x); % still arbitrary for translation etc.
        
    case 'division-mod' % Division model FitzGibbon 2001
        
        % coefs.
        b=sym('b',[1 7]); % b1 scale (fix), b2..8 homography, b9 nonlinearity
        
        % map
        Dmap= 1 / ( abs(b(7))* ...
            ( (        x(1) + b(1) * x(2) + b(2)   )^2   + ...
              (           abs(b(3))* x(2) + b(4)   )^2 ) + ...
              ( b(5) * x(1) + b(6) * x(2) + 1      )^2 )  ...
            * [        x(1) + b(1) * x(2) + b(2)
                          abs(b(3))* x(2) + b(4) ];
                  
    case 'interface-mod' % Interface model
        
        % coefs.
        b=sym('b',[1 7]); % b1 scale (fix), b2..8 homography, b9 nonlinearity
        
        % image mapping
        H=[ 1        b(1)  b(2)
            0    abs(b(3)) b(4)
            b(5)     b(6)  1 ];
        xl=H*[x
            1];
        
        % map
        Dmap= xl(1:2) / sqrt( abs(b(7))* ( xl(1)^2   + xl(2) ^2 ) + xl(3)^2 );
        
    case 'interface-cor' % Interface model with correction
        
        % order poly-coefficient
        [oy,ox]=meshgrid(0:2,0:2);
        T=(ox+oy)<=2 & (ox+oy)>1; % 2
        ox=ox(T);
        oy=oy(T);
        T=reshape(T,[],1);
        
        % order poly-coefficient
        [oyp,oxp]=meshgrid(0:ctrl.mord,0:ctrl.mord);
        Tp=(oxp+oyp)<=ctrl.mord & (oxp+oyp)>1;
        oxp=oxp(Tp);
        oyp=oyp(Tp);
        Tp=reshape(Tp,[],1);
        
        % coefs.
        b=sym('b',[1 7+nnz(T)+2*nnz(Tp)]); % 7+ 1*
        
        % image mapping
        H=[ 1     b(1)  b(2)
            0 abs(b(3)) b(4)
            b(5)  b(6)  1 ]; % no splitting Haff/Hproj promotes coupling
        xl=H*[ x
               (1+b(8:10)*((x(1).^ox).*(x(2).^oy))) ]; % 1 plus the correction - remove sqrt for impr. resul.
        
        % map
        x= xl(1:2) / sqrt( abs(b(7)) * ( xl(1)^2   + xl(2) ^2 ) + xl(3)^2 );
        
        % define mapping
        Dmap=x+reshape(b(11:end),2,[])*((x(1).^oxp).*(x(2).^oyp));
        
end

% export matlab function
Dmap=matlabFunction(Dmap,'Optimize',false);

%% save mapping function
save([folder date cal vsl 'Dmap.mat'],'Dmap')% save mapping

end