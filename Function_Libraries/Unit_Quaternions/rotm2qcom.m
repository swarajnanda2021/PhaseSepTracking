function q = rotm2qcom( R )
%rotm2quat convert a rotation matrix to a quaternion versor, suitable for 
%   symbolic input, note that a best rotation matrix can be obtain
%   beforehand using a svd or kabsch algortihm. Koen Muller

% make sure it is a rotation matrix
if size(R,3)==1
    R=frotm(R);
end

% component-wise from rotation formalism
qr=1/2*sqrt(1+R(1,1,:)+R(2,2,:)+R(3,3,:));
qi=(R(3,2,:)-R(2,3,:))./(4*qr);
qj=(R(1,3,:)-R(3,1,:))./(4*qr);
qk=(R(2,1,:)-R(1,2,:))./(4*qr);

% define total quaternion
q=[qr
    qi
    qj
    qk];

end