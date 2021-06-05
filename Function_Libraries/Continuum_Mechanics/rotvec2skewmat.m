function M = rotvec2skewmat(rotvec)
%skewmat2rotvec Skew matrix to rotation vector (use for spin or vorticity)
%
%   Input
%       rotvec Rotationvector notation 3xN
%
%   Output
%       M Skew symmetric (2x2) 3x3xN Matrix

% correct input
if size(rotvec,1)==2
    rotvec=[rotvec
        zeros(1,size(rotvec,2))]; % zeropad 3
end

% reshape data for matrix
rotvec=reshape(rotvec,size(rotvec,1),1,[]);

% rotation matrix
M=[  zeros(1,1,size(rotvec,3)) -rotvec(3,:,:)              rotvec(2,:,:)
     rotvec(3,:,:)              zeros(1,1,size(rotvec,3)) -rotvec(1,:,:)
    -rotvec(2,:,:)              rotvec(1,:,:)              zeros(1,1,size(rotvec,3))];

end

