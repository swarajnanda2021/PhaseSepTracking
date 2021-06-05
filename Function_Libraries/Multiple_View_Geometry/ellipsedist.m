function [ D ] = ellipsedist( C , Xp , typ )
%Cfunc Cost function between coordinate

% correct input
if nargin < 3
    typ='matrix';
end

% compute center coordinate ellipse
Xc=cvec2eshp(C);
    
% evaluate polynomial basis
xc=[Xc(1,:).^2
    Xc(1,:).*Xc(2,:)
    Xc(2,:).^2
    Xc(1,:)
    Xc(2,:)
    ones(size(Xc(1,:)))];
xp=[Xp(1,:).^2
    Xp(1,:).*Xp(2,:)
    Xp(2,:).^2
    Xp(1,:)
    Xp(2,:)
    ones(size(Xp(1,:)))];

% Compute distance by ellipse bodylength
switch typ
    case 'matrix'
        D=sqrt( abs( 1 - C'*xp ./ sum(C'.*xc',2) ) );
    case 'vector'
        D=sqrt( abs( 1 - sum(C.*xp,1) ./ sum(C.*xc,1) ) );
end

end

