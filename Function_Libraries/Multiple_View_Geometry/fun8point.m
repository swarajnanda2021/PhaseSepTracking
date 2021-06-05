function [ F ] = fun8point( x1 , x2 )
%fun_mat Compute the fundamental matrix by camera coordinates
% x2'Fx1

% data normalization
[T1,x1]=datnorm(x1);
[T2,x2]=datnorm(x2);

% transpose coordinates
x1=x1';
x2=x2';

% formulate coefficient matrix (dyadic structure, improve automation)
A=[ x1(:,1).*x2(:,1) x1(:,2).*x2(:,1) x1(:,3).*x2(:,1) ...
    x1(:,1).*x2(:,2) x1(:,2).*x2(:,2) x1(:,3).*x2(:,2) ...
    x1(:,1).*x2(:,3) x1(:,2).*x2(:,3) x1(:,3).*x2(:,3) ];

% solve l.s. problem on alg. dist, DLT solution
[~,~,V] = svd(A);

% solve h
f=V(:,end);

% in matrix form
F=reshape(f,3,[])';

% make singular by frobenius norm
[U,S,V]=svd(F);
S(end,end)=0;
F=U*S*V';

% data denormalization (here?)
F=T2'*F*T1;

end