function [ ang ] = rotm2eula( bas, seq )
%rotm2eula euler angles from rotm matrix
%   if not defined 2 basii, compute to rh cartesian basis
%   two vector basii can be propagated by this function
%       appereantly this is the mapping back to reference definition..
%       note: bas : {t,n,b}_local vs ref : {e1,e2,e3}_global

%reference
ref=eye(3);

% angle sequence
switch seq
    case 'ZXZ' % classical euler definition
        % first normal projections
        Z1=dot(bas(:,3),ref(:,1));
        Z2=dot(bas(:,3),ref(:,2));
        Z3=dot(bas(:,3),ref(:,3));
        X3=dot(bas(:,1),ref(:,3));
        Y3=dot(bas(:,2),ref(:,3));
        
        %euler angles
        ang=zeros(1,3);
        ang(1)=atan2(Z1,-Z2);
        ang(2)=acos(Z3);
        ang(3)=atan2(X3,Y3);
    case 'ZYX' % tait bryan angles
        % first normal projections
        X1=dot(bas(:,1),ref(:,1));
        X2=dot(bas(:,1),ref(:,2));
        X3=dot(bas(:,1),ref(:,3));
        Y3=dot(bas(:,2),ref(:,3));
        Z3=dot(bas(:,3),ref(:,3));
        
        %tait bryan angles
        ang=zeros(1,3);
        ang(1)=atan2(X2,X1);
        ang(2)=asin(-X3);
        ang(3)=atan2(Y3,Z3);
        
    case 'ZXY' % modified tait bryan angles
        % first normal projections
        Y1=dot(bas(:,2),ref(:,1));
        Y2=dot(bas(:,2),ref(:,2));
        Y3=dot(bas(:,2),ref(:,3));
        X3=dot(bas(:,1),ref(:,3));
        Z3=dot(bas(:,3),ref(:,3));
        
        % mod tait bryan angles
        ang=zeros(1,3);
        ang(1)=atan2(-Y1,Y2);
        ang(2)=asin(Y3);
        ang(3)=atan2(-X3,Z3);
    case 'XYX' % experimental classical definition
        % first normal projections
        X1=dot(bas(:,1),ref(:,1));
        X2=dot(bas(:,1),ref(:,2));
        X3=dot(bas(:,1),ref(:,3));
        Y1=dot(bas(:,2),ref(:,1));
        Z1=dot(bas(:,3),ref(:,1));
        
        %euler angles
        ang=zeros(1,3);
        ang(1)=atan2(X2,-X3); %atan2(y,x)
        ang(2)=acos(X1);
        ang(3)=atan2(Y1,Z1);
    otherwise
        error('non-listed rotation sequence..')
end

end