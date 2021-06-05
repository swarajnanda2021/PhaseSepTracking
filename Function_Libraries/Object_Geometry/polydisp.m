function [ coef,fdat ] = polydisp( tdat, ord , crt )
%polydisp Polynomial displacement field from a set of trajectories
%
%   Input
%       tdat data matrix [index timestamp X Y Z] 
%          trajectories all include same start frame 0
%       ord maximum polynomial order [X Y Z T]
%       crt: order of cross terms
%          none: 0
%          ord dim: 1 | 2 | 3 | 4 (time)
%          all: 5
%          array: custom
%
%   Output
%       coef coefficients [oX oY oZ oT cX cY cZ]
%       fdat fit data trajectories from initial point [Ni Xi fit-res]
%
%   Note: 
%       Supporting a grid indexing could impose simultaneous solution of
%       multiple vectors and impose consistentency contraints on the nodes.
%
%   We fit from the reference configuration to the current at time 0
%   configuration at time +/-
%
%   When initial frame misses in the trajectory data, the fit requires a
%   nonlinear optimization to find it itself, therefore not implemented.

% correct input
if nargin<2 % linear velocity field characteristics
    ord=[1 1 1 2]; % dudx dudy dudz dudt (u=dx/dt)
end
if length(ord)==3 % 2D
    ord=[ord(1:2) 0 ord(3)]; % add z=0 plane
end
if nargin<3
    if ord(4)==0 % spatial crossterms
        crt=1;
    else
        crt=4; % temporal crossterms
    end
end
if size(tdat,1)==4 % 2D
    tdat=[tdat
        zeros(1,size(tdat,2))]; % add z=0 plane
end

% terms polynomial model
[ox,oy,oz,ot]=ndgrid(0:ord(1),0:ord(2),0:ord(3),0:ord(4)); 

% polynomial cross terms
if numel(crt)==1
    if crt==0 || crt==5 % selected terms to multivariate polynomial fit
        T=ones(ord+1)*(crt~=0);
    else
        T=(ox+oy+oz+ot)<=ord(crt);
    end
    T(:,1,1,1)=1; T(1,:,1,1)=1; T(1,1,:,1)=1; T(1,1,1,:)=1;
    
    %constaint d(x,y,z,0)=0, no initial displacement by fitting to traj
    T(ot==0)=0;
    
else % custom
    T=crt;
end

% select polynomial order each term
Ord=[reshape(ox(T),1,[])
    reshape(oy(T),1,[])
    reshape(oz(T),1,[])
    reshape(ot(T),1,[])];

% define unique indexing
[ui,~,i]=unique(tdat(1,:));

% data size
sizX=size(tdat,1)-2;

% assemble n
N=tdat(2,:); % data

% assemble X
X=tdat(3:end,:); % data

% initial condition
G_dat=sparse(i,1:length(i),ones(size(i)),length(ui),length(i));
% [~,ind]=max(G_dat*spdiags(max(N)+1-N',0,length(N),length(N)),[],2); % find index first frame
Xi=X(:,N==0)*G_dat;
% Ni=N(ind)*G_dat;

% Compute basis function expansion
phi=prod( permute([Xi' N'],[1 3 2]).^permute(Ord,[3 2 1]) , 3 ) ; % ... - prod( permute([Xi' Ni'],[1 3 2]).^permute(Ord,[3 2 1]).*permute([0;0;0;1],[3 2 1]) , 3 );

% Assemble blockdiagonal for each physical dimension
Phi=kron(speye(sizX),phi);

% concatenate dimension data
x=reshape(X'-Xi',[],1); % T-I-X seq

% solve normal equations
M=Phi'*Phi; % vd Monde matrix
c=M\(Phi'*x); % solution

% compute residual
r=sqrt(sum(reshape((x-Phi*c).^2,[],sizX),2)');

% write coefficients
coef=[Ord
    reshape(c,[],sizX)'];

% write fit data
fdat=[N
    Xi
    r]; 

end

