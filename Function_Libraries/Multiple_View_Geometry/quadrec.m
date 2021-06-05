function [ q, il, ir ] = quadrec( c, P , t , typ)
%quadrec Quadric reconstruction using the direct (linear) method on the
%   adjoint / inverse quadric. After cross and zisserman 1998
%   reconstruction from dual space geometry.
%
%   When no triangulation vector is supported we reconstruct the quadric by
%   the conic outlines and the triangulation vector can be derived. In this
%   case we need at least 3 views.
%
%   When a triangulation vector is supported we reconstruct a hybrid
%   quadric which is constrained to the triangulation vector. In this case
%   we need at least 2 views.
%
%   Input
%       c conic vectors from multiple views
%           (6xN matrix )xC cell or 3xCxN matrix
%       P projection matrices from different views (3x4)xC cell
%       t prescribe triangulation vector 3xN matrix
%       typ type of reconstruction when triangulation prescribe
%           free | diag(onal) | spher(ical)
%
%   Output
%       q computed quadric vector 10xN matrix
%       l scaling for the inverse conic for each view CxN matrix
%       r residual inverse light intensity solution 1xN matrix
%
%   Remark: nan values in conic and triangulation vector will be ignored in
%           reconstruction.
%
%   Note: When prescribing the triangulation vector and constraining the
%         quadric spherical ensures a non-ruled levelset.

% correct input
if iscell(c) % reorder data
    c=permute(cat(3,c{:}),[1 3 2]); % ordered as [con cam obj]
end
if nargin<3 || isempty(t) % define nan triangulation vector
    t=nan*ones(3,size(c,3));
end
t=reshape(t,3,1,[]);
if nargin <4 % if not stated spherical always use free
    typ='free';
end

% data size
sizC=size(c,2); % number of cameras
sizN=size(c,3); % number of objects

% mark valid data
v=~isnan(c);       % valid data
v=( sum(v,2)>=( 3 - double( sum(~isnan(t),1)==3 ) ) ) & v; % by triangulation vector and number cameras

% compute inverse conic
ic=nan(size(c)); % initiate
ic(v)=contrans(c(v),[],'invert');

% compute projection matrices inverse quadric vector
Ac=[1 2 4
    2 3 5
    4 5 6]; % conic vector assembly matrix
Dc=pinv(full(sparse(Ac(:)',1:numel(Ac),ones(1,numel(Ac)),max(Ac(:)),numel(Ac)))); % conic vector duplication matrix
Aq=[1 2 3 7
    2 4 5 8
    3 5 6 9
    7 8 9 10]; % quadric vector assembly matrix
Dq=pinv(full(sparse(Aq(:)',1:numel(Aq),ones(1,numel(Aq)),max(Aq(:)),numel(Aq)))); % quadric vector duplication matrix
Pm=cellfun(@(P) pinv(Dc)*kron(P,P)*Dq ,P,'UniformOutput',false);

% split solution by supported triangulation vector
g_val=squeeze(sum(sum(v,1),2)>0)';
g_unc=g_val & squeeze(isnan(sum(t,1)))'; % unconstrained part
g_con=g_val & ~g_unc; % constrained part

% define constraint prescribe scale=1 or (weigted) translation vector t
iq10=[zeros(9,sizN) ; double(g_unc)]; % impose scale
iqt=qmat2qvec([t.*permute(t,[2 1 3]) t; permute(t,[2 1 3]) double(reshape(g_con,1,1,sizN))]); % impose translation vector
iq0=reshape(sum(cat(3,iqt,iq10),3,'omitnan'),[],1); % forcing vector inversion

% assemble data matrices for projection and inverse conic matrices
Pm=kron(speye(sizN),cat(1,Pm{:}));
pm=Pm*iq0;
ic=reshape(ic,6,1,sizN*sizC);
i=repmat((1:size(ic,1))',1,size(ic,2))+(size(ic,1)*(reshape(1:size(ic,3),1,1,[])-1));
j=repmat((1:size(ic,2)),size(ic,1),1)+(size(ic,2)*(reshape(1:size(ic,3),1,1,[])-1));
ic=sparse(i(:),j(:),ic(:),size(ic,1)*size(ic,3),size(ic,2)*size(ic,3));

% define data validity to block nan values
V=reshape(v,6,1,sizN*sizC);
i=repmat((1:size(V,1))',1,size(V,2))+(size(V,1)*(reshape(1:size(V,3),1,1,[])-1));
j=repmat((1:size(V,2)),size(V,1),1)+(size(V,2)*(reshape(1:size(V,3),1,1,[])-1));
V=sparse(i(:),j(:),V(:),size(V,1)*size(V,3),size(V,2)*size(V,3));
Vd=sum(V,2)'>0; % valid data in assembly

% initiate inverse solution
iq=reshape(iq0,10,sizN); % inverse quadric to add solution to
il=nan(sizC,sizN); % inverse scaling
ir=nan(1,sizN); % inverse residual

iq(:,~g_val)=nan; % force

%-- 1) solve unconstrained problem

if nnz(g_unc)>0
    
    % define unknowns
    iqu=reshape([ones(9,1) ; zeros(1,1)] & g_unc,[],1);
    Punc=Pm(:,iqu);
    Pq=Vd*abs(Punc)~=0; % valid q in solution vector
    Pd=Pq*abs(Punc')~=0 & Vd;
    Pl=Pd*abs(ic)~=0 & ~isnan(Pd*abs(ic));
    
    % solve linear equations inverse quadric
    A=cat(2,ic(Pd,Pl),-Punc(Pd,Pq));
    B=A'*A;
    b=-A'*pm(Pd);
    y=-B\b; % solution vector
    r=(A*y-pm(Pd)).^2; % residual per data point
    
    % intiate \\ camera number can vary
    res=nan*ones(size(Pd));
    sol=nan*ones(size(Pq));
    lam=nan*ones(size(Pl)); % nan*ones(size(Pl));
    
    % write
    res(Pd)=r;
    sol(Pq)=y(nnz(Pl)+1:end);
    lam(Pl)=y(1:nnz(Pl));
    
    % resize
    lam=reshape(lam,sizC,[]);
    res=sqrt(sum(reshape(res,6*sizC,[]),1,'omitnan'));
    
    % split solution and write the unknowns
    iq(iqu)=iq(iqu)+sol';%y(nnz(Pl)+(1:nnz(Pq)));
    il(:,g_unc)=lam(:,g_unc);
    ir(g_unc)=res(g_unc); %3 6
    
end

%-- 2) solve constrained triangulation vector

if nnz(g_con)>0
    
    % switch constraint cases
    switch typ
        case 'free'
            
            % define unknowns
            iqu=reshape([ ones(6,1) ; zeros(4,1)] & g_con,[],1);
            Pcon=Pm(:,iqu);
            
        case 'diag'
            
            % define unknowns
            iqu=reshape([ [1 0 0 1 0 1]' ; zeros(4,1)] & g_con,[],1);
            Pcon=Pm(:,iqu);
            
        case 'spher'
            
            % define unknowns
            iqu=kron(spdiags(g_con',0,length(g_con),length(g_con)),[ [1 0 0 1 0 1]' ; zeros(4,1)]);%reshape([ [1 0 0 1 0 1]' ; zeros(4,1)] & g_con,[],1);
            iqu=iqu(:,sum(iqu,1)>0);
            Pcon=Pm*double(iqu);
            
            % redefine
            iqu=reshape([ [1 0 0 1 0 1]' ; zeros(4,1)] & g_con,[],1);
            
    end
    Pq=Vd*abs(Pcon)~=0; % valid q in solution vector
    Pd=Pq*abs(Pcon')~=0 & Vd;
    Pl=Pd*abs(ic)~=0 & ~isnan(Pd*abs(ic));
    
    % solve linear equations inverse quadric
    A=cat(2,ic(Pd,Pl),-Pcon(Pd,Pq));
    B=A'*A;
    b=-A'*pm(Pd);
    y=-B\b; % solution vector
    r=(A*y-pm(Pd)).^2; % residual per data point
    
    % intiate
    res=nan*ones(size(Pd));
    sol=nan*ones(size(Pq));
    lam=nan*ones(size(Pl));
    
    % write
    res(Pd)=r;
    sol(Pq)=y(nnz(Pl)+1:end);
    lam(Pl)=y(1:nnz(Pl));
    
    % replicate spherical case
    if strcmp(typ,'spher')
        sol=reshape(repmat(sol,3,1),1,[]);
    end
    
    % resize
    lam=reshape(lam,sizC,[]);
    res=sqrt(sum(reshape(res,6*sizC,[]),1,'omitnan'));
    
    % split solution and write the unknowns
    iq(iqu)=iq(iqu)+sol';%y(nnz(Pl)+(1:nnz(Pq)));
    il(:,g_con)=lam(:,g_con);
    ir(g_con)=res(g_con); %3 6
    
end

% correct
il(:,isnan(sum(iq,1)))=nan;
ir(:,isnan(sum(iq,1)))=nan;
iq(:,isnan(sum(iq,1)))=nan;

% compute quadric
v=~isnan(iq);
q=nan(size(iq));
q(v)=quadtrans(reshape(iq(v),10,[]),[],'invert');

end

% i=repmat((1:size(C,1))',1,size(C,2))+(size(C,1)*(reshape(1:size(C,3),1,1,[])-1));
% j=repmat((1:size(C,2)),size(C,1),1)+(size(C,2)*(reshape(1:size(C,3),1,1,[])-1));
% iC=inv(sparse(i(:),j(:),C(:),size(C,1)*size(C,3),size(C,2)*size(C,3)));
% iC=reshape(full(iC(sub2ind(size(iC),i(:),j(:)))),size(C,1),size(C,2),size(C,3));

% i=repmat((1:size(iQ,1))',1,size(iQ,2))+(size(iQ,1)*(reshape(1:size(iQ,3),1,1,[])-1));
% j=repmat((1:size(iQ,2)),size(iQ,1),1)+(size(iQ,2)*(reshape(1:size(iQ,3),1,1,[])-1));
% Q=inv(sparse(i(:),j(:),iQ(:),size(iQ,1)*size(iQ,3),size(iQ,2)*size(iQ,3)));
% Q=reshape(full(Q(sub2ind(size(Q),i(:),j(:)))),size(iQ,1),size(iQ,2),size(iQ,3));

%     if 0
%
%     % select the constraint for the trace
% %     s=s(nnz(Pl)+1:end);
%
%     % define constraint vector
% %     a=iqu;%kron(speye(nnz(g_con)),[ [1 0 0 1 0 1]' ; zeros(4,1)]);%[zeros(nnz(Vl),1) ; double(iqu)];
%
%     %-- 3) solve constrained quadric with spherical volume
%
%     a=kron(speye(nnz(g_con)),[ [1/2 0 0 1/2 0 -1]' ; zeros(4,1)]);%[zeros(nnz(Vl),1) ; double(iqu)];
%
%     s=zeros(nnz(g_con),1);
%
%     % define unknowns
%     iqu=reshape([ ones(6,1) ; zeros(4,1)] & g_con,[],1);
%     Pcon=P(:,iqu);
%     Pq=Vd*abs(Pcon)~=0; % valid q in solution vector
%     Pd=Pq*abs(Pcon')~=0;
%     Pl=Pd*abs(ic)~=0;
%
%     % solve linear equations inverse quadric
%     A=cat(2,ic(Pd,Pl),-Pcon(Pd,Pq));
%     B=A'*A;
%     b=-A'*p(Pd);
%     H=[B a; a' zeros(numel(s))];
%     f=[b ; -s];%/3
%     y=-H\f; % solution vector
%     y=y(1:length(b)); % trash langrange multiplier
%     r=(A*y-p(Pd)).^2; % residual per data point
%
%     end
%
