function [ c ] = eshp2cvec( x,ax,ang )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% translation vector
t=reshape(x,2,1,[]);

% rotation matrix
R=nan(2,2,size(t,3));
R(:)=cell2mat(reshape(cellfun(@(x)[cos(x) -sin(x) ; sin(x) cos(x)],num2cell(ang,1),'UniformOutput',false),1,1,[])); % single rotation according eula2rotm ZYX

% quadric shape
C0=nan(3,3,size(t,3));
C0(:)=cell2mat(reshape(cellfun(@(x)[diag(x.^(-2)) zeros(2,1); zeros(1,2) -1],num2cell(ax,1),'UniformOutput',false),1,1,[]));

% homography
H=nan(3,3,size(t,3));
H(:)=cell2mat(reshape(cellfun(@(R,t)[R t ; zeros(1,2) 1],num2cell(R,[1 2]),num2cell(t,[1 2]),'UniformOutput',false),1,1,[]));

% conic
C=nan(3,3,size(t,3));
C(:)=cell2mat(cellfun(@(H,Q0)inv(H)'*Q0*inv(H),num2cell(H,[1 2]),num2cell(C0,[1 2]),'UniformOutput',false));

% vector
c=cmat2cvec(C);

end

