function [ coef ,res] = polytraj( dat, ord , bias )%
%polytraj Polynomial trajectory to set of points 2, 3 and N dimensional
%   fit many object simultaneously
%   min_{c^k}sum_c||x^n_i-sum_k(c^k_i*(t^n)**k)||_2^2 -> ..||_F -> ..||_2
%   
%   Input:
%       dat data matrix [index timestamp X Y Z]
%       ord maximum polynomial order
%       bias to order, bias less order to points to regress
%   
%   Output
%       coef coefficients [index coefnumber cX cY cZ] unlike timestamp
%           coef is always complete and zeros padded for missing
%           information
%       res residual per point on the trajectory
%

% correct input
if nargin<3
    bias=0;
end

% define unique indexing
[ui,~,ki]=unique(dat(1,:));
[ut,~,kt]=unique(dat(2,:));

% data size
sizO=ord+1;
sizI=length(ui);
sizT=length(ut);
sizX=size(dat,1)-2;

%-- #) Assemble data

% assemble t
t=zeros(1,sizT,sizI); % zero padding
t(sub2ind(size(t),ones(size(kt)),kt,ki))=dat(2,:); % data

% assemble X
X=zeros(sizX,sizT,sizI); % zero padding
for i=1:sizX
    X(sub2ind(size(X),i*ones(size(kt)),kt,ki))=dat(2+i,:); % data
end

% construct validity data - is data present at that frame?
v=zeros(1,sizT,sizI); % zeros
v(sub2ind(size(t),ones(size(kt)),kt,ki))=1; % indicator / weigth

% construct order
o=repmat((1:sizO),1,1,sizI);

% construct index
ind=zeros(1,1,sizI);
ind(:)=ui;
ind=repmat(ind,1,sizO,1);

%-- #) Preconditioning to avoid ill-conditioning for inversion

% precondition time
prec.Tmed=mean(t,2); % middle position, but mean results better
t=(t-prec.Tmed); % shift positions and scale (!)
prec.Tran=1+std(abs(t),[],2); % rang in -2*medvar to 2*medvar, but std +1 results better
t=t./prec.Tran; % rescale per track (note decouples)

% % precondition space (not strictly neccesary, but for completeness)
% prec.Xmean=sum(X,2)./sum(v,2); % middle position
% X=(X-prec.Xmean); % shift positions
% prec.Xstd=eps+sqrt(sum(X.^2,2)./sum(v,2)); % domain in -sigma to ~sigma
% X=X./sqrt(sum(prec.Xstd.^2,1)); % rescale per track (note decouples, global solution not affected)

%-- #) Assemble inversion problem

% construct time matrix
V=permute(o,[2 1 3])<=(sum(v,2)-(sum(v,2)~=1)*bias) & repmat(v,sizO,1,1);
T=repmat(t,sizO,1,1).^repmat((0:ord)',1,sizT,sizI); % [ord frame track]

% assemble sparse matrices
i=repmat((1:sizO)',1,sizT)+(sizO*(reshape(1:sizI,1,1,[])-1));
j=repmat((1:sizT),sizO,1)+(sizT*(reshape(1:sizI,1,1,[])-1));
T=sparse(i(:),j(:),T(:),sizO+sizO*(sizI-1),sizT+sizT*(sizI-1)); % overwrite
V=sparse(i(:),j(:),V(:),sizO+sizO*(sizI-1),sizT+sizT*(sizI-1)); % overwrite

% validity
vT=sum(V,1)>0.5; % valid k in solution vector
vO=sum(V,2)'>0.5; % valid data in assembly

% further assemble dimension
T=kron(speye(sizX),T); % block diagonal for each physical dimension
vT=repmat(vT,1,sizX);
vO=repmat(vO,1,sizX);

% concatenate dimension data
b=reshape(permute(X,[2 3 1]),[],1); % T-I-X seq

%-- #) Solve trajectories fit(s)

% initiate coeffients and residuals
coef=zeros(size(vO')); % coef data
res=zeros(size(vT')); % residual data

% solve normal equations
M=T(vO,vT)*T(vO,vT)'; % vd Monde matrix
coef(vO)=M\(T(vO,vT)*b(vT)); % solution

%figure; spy(M); condest(M)

% % solve normal equations
% M=T(vO,vT)*T(vO,vT)'; % vd Monde matrix
% [P,R,C] = equilibrate(M); ! slow
% B = R*P*M*C;
% coef(vO)=C*(B\(R*P*T(vO,vT)*b(vT))); % solution
% 
% %figure; spy(B); condest(B)

% compute residuals for each trajectory point and dimension
res(vT)=b(vT)-T(vO,vT)'*coef(vO);

% reshape coeffient data
coef=permute(reshape(coef,sizO,sizI,sizX),[3 1 2]);

% reshape total residual data
res=sqrt(sum(permute(reshape(res,sizT,sizI,sizX),[3 1 2]).^2,1));

%-- #) Rescale and write output

% binomial coefficients for shifting the trajectorie time stamps
Cnk=zeros(sizO,sizO); % zeros supress unnnecesary terms
nk=tril((0:ord)' - (0:ord)); % lower triangle is valid
for n=1:sizO
    for k=1:n
        Cnk(n,k)=nchoosek(n-1,k-1);
    end
end

% temporal rescaling coeffient data
coef=coef./(prec.Tran.^(o-1)); % rescale
coef=permute(sum( ...
    repmat(permute(coef,[1 3 2]),1,1,1,sizO) ...
    .*repmat(permute(Cnk,[3 4 1 2]),sizX,sizI,1,1) ...
    .*( repmat(permute(-prec.Tmed,[1 3 2]),sizX,1,sizO,sizO) ...
    .^repmat(permute(nk,[3 4 1 2]),sizX,sizI,1,1) )...
    ,3),[1 4 2 3]); % shift using the binomial expansion

% % spatial rescaling coeffient data
% coef=coef.*sqrt(sum(prec.Xstd.^2,1)); % rescale
% coef(:,1,:)=coef(:,1,:)+prec.Xmean; % shift center

% format output
coef=[reshape(ind,1,[])
    reshape(o,1,[])
    reshape(coef,sizX,[])];% [index coefnumber cX cY cZ]

% % spatial rescaling residual data
% res=res.*sqrt(sum(prec.Xstd.^2,1));

% write reprojection data
res=reshape(res(sub2ind(size(t),ones(size(kt)),kt,ki)),1,[]) ;

end

