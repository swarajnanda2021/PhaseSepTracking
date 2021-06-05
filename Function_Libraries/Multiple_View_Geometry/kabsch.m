function [ R,t, S] = kabsch( x1 , x2 )
%kabsch compute the best rotation (and translation) by use of the Kabsch
%   algorithm. Up to computer precision.

% compute midpoint data
x1m=mean(x1,2);
x2m=mean(x2,2);

% center data
x1n=x1-x1m;
x2n=x2-x2m;

% define covariance matrix
C=x2n*x1n'; % dot(...,dim)

% compute svd
[U,S,V]=svd(C);

% constrain proper
D=eye(size(S));
D(end)=det(U*V');

% define rotation
R=U*D*V'; 

% get translation
t=mean(x2-R*x1,2);

end