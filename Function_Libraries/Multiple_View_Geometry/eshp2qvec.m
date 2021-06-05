function [ q ] = eshp2qvec( x,ax,ang )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% translation vector
t=reshape(x,3,1,[]);

% rotation matrix
R=nan(3,3,size(t,3));
R(:)=cell2mat(reshape(cellfun(@(x)eula2rotm(x,'ZXY'),num2cell(ang,1),'UniformOutput',false),1,1,[]));

% quadric shape
Q0=nan(4,4,size(t,3));
Q0(:)=cell2mat(reshape(cellfun(@(x)[diag(x.^(-2)) zeros(3,1); zeros(1,3) -1],num2cell(ax,1),'UniformOutput',false),1,1,[]));

% homography
H=nan(4,4,size(t,3));
H(:)=cell2mat(reshape(cellfun(@(R,t)[R t ; zeros(1,3) 1],num2cell(R,[1 2]),num2cell(t,[1 2]),'UniformOutput',false),1,1,[]));

% quadric
Q=nan(4,4,size(t,3));
Q(:)=cell2mat(cellfun(@(H,Q0)inv(H)'*Q0*inv(H),num2cell(H,[1 2]),num2cell(Q0,[1 2]),'UniformOutput',false));

% vector
q=qmat2qvec(Q);

end

