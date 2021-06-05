function [ c ] = cmat2cvec( C )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% get parameters, like a voight type vector
c=squeeze([C(1,1,:)
    C(1,2,:)*2
    C(2,2,:)
    C(1,3,:)*2
    C(2,3,:)*2
    C(3,3,:)]); % store as list

end