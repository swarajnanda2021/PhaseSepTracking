function [ D ] = quaddist( Q , Xp , typ )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% correct input
if nargin < 3
    typ='matrix';
end

% compute center coordinate ellipse
Xc=qvec2eshp(Q);

% evaluate polynomial basis
xc=[Xc(1,:).^2
    Xc(1,:).*Xc(2,:)
    Xc(1,:).*Xc(3,:)
    Xc(2,:).^2
    Xc(2,:).*Xc(3,:)
    Xc(3,:).^2
    Xc(1,:)
    Xc(2,:)
    Xc(3,:)
    ones(size(Xc(1,:)))];
xp=[Xp(1,:).^2
    Xp(1,:).*Xp(2,:)
    Xp(1,:).*Xp(3,:)
    Xp(2,:).^2
    Xp(2,:).*Xp(3,:)
    Xp(3,:).^2
    Xp(1,:)
    Xp(2,:)
    Xp(3,:)
    ones(size(Xp(1,:)))];

% Compute distance by ellipsoid bodylength
switch typ
    case 'matrix'
        D=sqrt( abs( 1 - Q'*xp ./ sum(Q'.*xc',2) ) );
    case 'vector'
        D=sqrt( abs( 1 - sum(Q.*xp,1) ./ sum(Q.*xc,1) ) );
end

end

