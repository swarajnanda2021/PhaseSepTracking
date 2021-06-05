function [ H ] = dirlintrans( x1 , x2 , Ceq0 , Cnorm )
%hom_graph compute camera homography of a set of points on a plane
%   compact linear independent comp.

% data normalization
[T1,x1]=datnorm(x1);
[T2,x2]=datnorm(x2);% formulate this as a mapping s.t. constraint can be indep. ?

% formulate coefficient matrix
A=sparse(...
    [zeros(size(x1')) repmat(-x2(3,:)',1,size(x1,1)).*x1' repmat(x2(2,:)',1,size(x1,1)).*x1'
    repmat(x2(3,:)',1,size(x1,1)).*x1' zeros(size(x1')) repmat(-x2(1,:)',1,size(x1,1)).*x1'
    repmat(-x2(2,:)',1,size(x1,1)).*x1' repmat(x2(1,:)',1,size(x1,1)).*x1' zeros(size(x1'))]...
    );

% Constraints
if nargin==4 && ~isempty(Cnorm)
    [~,D,cV]=svd(Cnorm); % something
    r=rank(D);
    A=A*cV;
    A1=A(:,1:r);
    A2=A(:,r+1:end);
    D1=D(1:r,1:r);
    A=(A2*pinv(A2)-eye(size(A2*pinv(A2))))*A1*inv(D1);
end
if nargin==3 && ~isempty(Ceq0)
    if size(Ceq0,1)<size(Ceq0,2)
        Ceq0=[Ceq0
            zeros(size(Ceq0,2)-size(Ceq0,1),size(Ceq0,2))];
    end
    [~,D,Cper]=svd(Ceq0);
    d=diag(D);
    Cper(:,d~=0)=0;
    A=A*sparse(Cper);
end

% solve l.s. problem on alg. dist, DLT solution
[~,~,V] = svds(sparse(A),1,'smallest');

% solve h
h=V(:,end);

% dewarp constraints
if nargin==4 && ~isempty(Cnorm)
    h1=D1\h;
    h2=-pinv(A2)*A1*h1;
    h=[h1
        h2];
    h=cV*h;
end
if nargin==3 && ~isempty(Ceq0)
    h=Cper*h;
end

% in matrix form
H=reshape(h,size(x1,1),size(x2,1))';

% data denormalization
H=T2\H*T1;

% normalize matrix scale, which is free and positive 
H=H/H(end);

end