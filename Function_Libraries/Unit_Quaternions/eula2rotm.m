function [ R ] = eula2rotm( ang, seq)
%eula2rotm rotation matrix from euler angles, can be used for symbolic
%   input ang in radians

if nargin <2
    seq='ZYX';
end

switch seq
    case 'ZXZ' % classical
        R=Rz(ang(1))*Rx(ang(2))*Rz(ang(3)); % ZYX
    case 'ZYX' % tait bryan angles
        R=Rz(ang(1))*Ry(ang(2))*Rx(ang(3)); % ZYX
    case 'ZXY' % modified tait bryan angles
        R=Rz(ang(1))*Rx(ang(2))*Ry(ang(3)); % ZYX
    case 'ZYZ' % experiments usefull
        R=Rz(ang(1))*Ry(ang(2))*Rz(ang(3)); % ZYZ
    otherwise
        error('non-listed rotation sequence..')
end

end

function R=Rx(ang) % can be replace by rotmx etc..
R = [1 0         0
     0 cos(ang) -sin(ang)
     0 sin(ang)  cos(ang)];
end

function R=Ry(ang)
R = [cos(ang) 0 sin(ang)
     0        1  0
     -sin(ang) 0  cos(ang)];
end

function R=Rz(ang)
R = [cos(ang) -sin(ang) 0
     sin(ang)  cos(ang) 0
     0         0        1];
end