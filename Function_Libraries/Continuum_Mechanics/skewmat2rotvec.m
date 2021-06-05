function rotvec = skewmat2rotvec(M)
%skewmat2rotvec Skew matrix to rotation vector (use for spin or vorticity)
%
%   Input
%       M Skew symmetric (2x2) 3x3xN Matrix
%
%   Output
%       rotvec Rotationvector notation 3xN

% correct input
if size(M,1)==2
    M=[M
        zeros(1,size(M,2))]; % zeropad 3
end
if size(M,2)==2
    M=[M zeros(size(M,1),1)]; % zeropad 3
end

% rotation vector
rotvec=[-M(2,3,:)
    M(1,3,:)
    -M(1,2,:)];

end

