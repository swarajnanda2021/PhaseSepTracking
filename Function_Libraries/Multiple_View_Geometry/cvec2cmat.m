function [ C ] = cvec2cmat( c )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Reshape conic parameters in vector def list
c=reshape(c,6,1,[]);

% Conic matrix
C=[c(1,:,:) c(2,:,:)/2 c(4,:,:)/2
    c(2,:,:)/2 c(3,:,:) c(5,:,:)/2
    c(4,:,:)/2 c(5,:,:)/2 c(6,:,:)];

end