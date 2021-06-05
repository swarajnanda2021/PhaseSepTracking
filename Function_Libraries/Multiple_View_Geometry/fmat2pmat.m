function [ P1,P2 ] = fmat2pmat( F )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% matrix 1, reference
P1=eye(3,4);

% matrix, relative
[U,~,~]=svd(F);
e2=U(:,end);

E2=[0 -e2(3) e2(2) 
    e2(3) 0 -e2(1)
    -e2(2) e2(1) 0];
P2=[E2*F,e2]; % here there was a fix by transpose

end