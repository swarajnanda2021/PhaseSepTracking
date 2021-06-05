function [ q ] = qmat2qvec( Q )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% get parameters, like a voight type vector
q=squeeze([Q(1,1,:)
    Q(1,2,:)*2
    Q(1,3,:)*2
    Q(2,2,:)
    Q(2,3,:)*2
    Q(3,3,:)
    Q(1,4,:)*2
    Q(2,4,:)*2
    Q(3,4,:)*2
    Q(4,4,:)]); % store as list

end