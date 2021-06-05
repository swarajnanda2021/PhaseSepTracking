function [ Q ] = qvec2qmat( q )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Reshape conic parameters in vector def list
q=reshape(q,10,1,[]);

% Conic matrix
Q=[q(1,:,:) q(2,:,:)/2 q(3,:,:)/2 q(7,:,:)/2
    q(2,:,:)/2 q(4,:,:) q(5,:,:)/2 q(8,:,:)/2
    q(3,:,:)/2 q(5,:,:)/2 q(6,:,:) q(9,:,:)/2
    q(7,:,:)/2 q(8,:,:)/2 q(9,:,:)/2 q(10,:,:)];

end

