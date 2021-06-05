function [ c ] = lset2cvec( H,g,c0 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

c=cmat2cvec(cat(1,cat(2,H/2,g/2),cat(2,permute(g/2,[2 1 3]),c0)));

end

