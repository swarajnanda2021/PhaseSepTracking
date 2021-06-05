function [ F , e2 , M ] = pmat2fmat( P1,P2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% camera centre
C=null(P1);

% pseudoinverse
piP1=pinv(P1);

% epipole 2
e2=P2*C;
E2=[0 -e2(3) e2(2)
    e2(3) 0 -e2(1)
    -e2(2) e2(1) 0];

% submatrix
M=P2*piP1;

% fundamental matrix
F=E2*M;

end