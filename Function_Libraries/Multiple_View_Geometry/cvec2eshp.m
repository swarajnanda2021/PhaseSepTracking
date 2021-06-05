function [ x,ax,ang ] = cvec2eshp( c )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   converts a conis vector of any levelset to an ellipse shape
%   saddles will be flipped axis
% would be nice if this function conserves input data dimensions to best
% extent

% select valid part
g_val=~isnan(sum(c,1));

% remove invalid
c=c(:,g_val);

% initiate
x=nan*zeros(2,length(g_val));
ax=nan*zeros(2,length(g_val));
ang=nan*zeros(1,length(g_val));

% vector to levelset
[ H,g,s ] = cvec2lset(c);

% fast triangulation : performance!
if ~isempty(s)
    
    % create sparse numbering
    [r,c,b]=size(H);
    i=repmat((1:r)',1,c)+(r*(reshape(1:b,1,1,[])-1));
    j=repmat((1:c),r,1)+(c*(reshape(1:b,1,1,[])-1));
    A=sparse(i(:),j(:),H(:),max(i(:)),max(j(:)));
    
    % fast midpoint triangulation
    x(:,g_val)=reshape(-A\g(:),2,[]); % x=squeeze(cell2mat(cellfun(@(H,g)-H\g,num2cell(H,[1 2]),num2cell(g,[1 2]),'UniformOutput',false)));

% else
%     
%     % dimensional correction when empty
%     x=zeros(2,0);
    
end

% decompose shape
if nargout>1
    
    % initiate variables
%     ax=zeros(size(x));
%     theta=zeros(1,size(x,2));
    gind=find(g_val);
    for i=1:size(H,3) % loop data
        
        % decompose[R,D,~]=svd(H(:,:,i)); % 
        [R,D]=schur(H(:,:,i)/2); % we dont care about V flipping axis, we want to vis. a neighborhood
        d=s(:,:,i)+1/2*g(:,:,i)'*x(:,gind(i));
        D=-D/d;%
        [D,j]=sort(diag(D),'descend'); % sort
        R=[R(:,j) zeros(2,1); zeros(1,2) 1];
        
        % write
        ang(gind(i))=[1 0 0]*rotm2eula(R,'ZXY')'; % pitch yaw roll
        ax(:,gind(i))=sqrt(1./D);
    end
    
%     % dimensional correction when empty
%     if isempty(ax)
%         ax=zeros(2,0);
%     end
%     if isempty(theta)
%         theta=zeros(1,0);
%     end
    
end

end

% % determinant
% detH=c(1,:).*c(3,:)-1/4*c(2,:).^2; % determinant hessian in conic
% trH=c(1,:)+c(3,:);
% discr=sqrt(trH.^2-4*detH);
% 
% % midpoint
% x=[(2*c(3,:).*c(4,:)-c(2,:).*c(5,:))./(-4*detH)
%     (2*c(1,:).*c(5,:)-c(2,:).*c(4,:))./(-4*detH)];
% 
% % axis
% detC=detH.*c(6,:)+1/4*(c(2,:).*c(5,:).*c(4,:)...
%     -c(3,:).*c(4,:).^2-c(1,:).*c(5,:).^2); % real ellipse when < 0
% a=[sqrt(-1/2*detC.*(trH+discr))./(detH)
%     sqrt(-1/2*detC.*(trH-discr))./(detH)];
% 
% % orientation
% thet=atan((c(3,:)-c(1,:)-discr)./c(2,:));
% thet(c(2,:)==0 & c(1,:)<c(3,:))=0;
% thet(c(2,:)==0 & c(1,:)>c(3,:))=pi/2;
