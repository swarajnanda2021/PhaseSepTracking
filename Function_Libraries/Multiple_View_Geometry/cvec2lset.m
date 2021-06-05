function [ H,g,s,L ] = cvec2lset( c,x )
%cvec2lset conic vector to level set
%
%   Input:
%       c Conic vector 6xN
%       x Point 2xN
%
%   Output:
%       H Hessian 2x2xN
%       g Gradient 2x1xN
%       s Scale 1x1xN
%       L Levelset value 1xN

% Reshape conic parameters in vector def list
c=reshape(c,6,1,[]);

% Hessian matrix
H=[2*c(1,:,:) c(2,:,:)
    c(2,:,:) 2*c(3,:,:)];

% gradient vector
g=[c(4,:,:)
    c(5,:,:)];

% offset
s=c(6,:,:);

% evaluate
if nargin>1
    c=reshape(c,6,[]);
    x=[x(1,:).^2
        x(1,:).*x(2,:)
        x(2,:).^2
        x(1,:)
        x(2,:)
        ones(size(x(1,:)))];
    L=sum(c.*x,1);
end

end

