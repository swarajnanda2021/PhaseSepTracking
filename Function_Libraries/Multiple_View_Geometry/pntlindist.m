function [ D , y ] = pntlindist( x , l, typ )
%pntlindist Point to nearest point on line distance, compute this distance
%   as a covariance between x in row and l in column space.
%       input,
%           x point
%           l line

% correct input
if size(x,1)==2
    x=inhc2homc(x);
end
if size(l,1)==2
    l=inhc2homc(l);
end
if nargin < 3
    typ='matrix';
end

% compute distance
switch typ
    case 'matrix'
        D=x'*l./sqrt(sum(l(1:2,:).^2,1));
    case 'vector'
        D=sum(x.*l,1)./sqrt(sum(l(1:2,:).^2,1));
end

% compute coordinate on line
if nargout>1
    
    switch typ
        case 'matrix'
            y=cat(3,(l(2,:).*(x(1,:)'*l(2,:)-x(2,:)'*l(1,:))-l(1,:).*l(3,:))./sum(l(1:2,:).^2,1),...
                (l(1,:).*(-x(1,:)'*l(2,:)+x(2,:)'*l(1,:))-l(2,:).*l(3,:))./sum(l(1:2,:).^2,1));
            y=permute(y,[3 1 2]);
        case 'vector'
            y=cat(3,(l(2,:).*(x(1,:).*l(2,:)-x(2,:).*l(1,:))-l(1,:).*l(3,:))./sum(l(1:2,:).^2,1),...
                (l(1,:).*(-x(1,:).*l(2,:)+x(2,:).*l(1,:))-l(2,:).*l(3,:))./sum(l(1:2,:).^2,1));
            y=permute(y,[3 2 1]);
    end
end

end