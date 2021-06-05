function [ R , S ] = frotm( Q )
%frotm fit best rotation matrix to input matrix
%   much like kabsch but different definition for input, therefore not
%   centering. Up to computer precision

% compute signular value dec.
[U,S,V]=svd(Q);

% constrain to be proper
D=eye(size(S));
D(end)=det(U*V');

% best rotation
R=U*D*V';

end

