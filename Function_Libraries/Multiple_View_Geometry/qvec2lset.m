function [ H,g,s,L ] = qvec2lset( q,x )
%qvec2lset conic vector to level set
%
%   Input:
%       q Conic vector 10xN
%       x Point 3xN
%
%   Output:
%       H Hessian 3x3xN
%       g Gradient 3x1xN
%       s Scale 1x1xN
%       L Levelset value 1xN

% Reshape conic parameters in vector def list
q=reshape(q,10,1,[]);

% Hessian matrix
H=[2*q(1,:,:) q(2,:,:) q(3,:,:)
    q(2,:,:) 2*q(4,:,:) q(5,:,:)
    q(3,:,:) q(5,:,:) 2*q(6,:,:)];

% gradient vector
g=[q(7,:,:)
    q(8,:,:)
    q(9,:,:)];

% offset
s=q(10,:,:);

% evaluate
if nargin>1
    q=reshape(q,10,[]);
    x=[x(1,:).^2
        x(1,:).*x(2,:)
        x(1,:).*x(3,:)
        x(2,:).^2
        x(2,:).*x(3,:)
        x(3,:).^2
        x(1,:)
        x(2,:)
        x(3,:)
        ones(size(x(1,:)))];
    L=sum(q.*x,1);
end

end

