function [ x, ax, ang ] = qvec2eshp( q )
%UNTITLED Summary of this function goes here
%   converts a quadric vector of any levelset to an ellipse shape
%   saddles will be flipped axis

% select valid part
g_val=~isnan(sum(q,1));

% remove invalid
q=q(:,g_val);

% initiate
x=nan*zeros(3,length(g_val));
ax=nan*zeros(3,length(g_val));
ang=nan*zeros(3,length(g_val));

% vector to levelset
[ H,g,s ] = qvec2lset(q);

% fast triangulation : performance!
if ~isempty(s) 
    
    % create sparse numbering
    [r,c,b]=size(H);
    i=repmat((1:r)',1,c)+(r*(reshape(1:b,1,1,[])-1));
    j=repmat((1:c),r,1)+(c*(reshape(1:b,1,1,[])-1));
    A=sparse(i(:),j(:),H(:),max(i(:)),max(j(:)));
    
    % fast midpoint triangulation
    x(:,g_val)=reshape(-A\g(:),3,[]); % x=squeeze(cell2mat(cellfun(@(H,g)-H\g,num2cell(H,[1 2]),num2cell(g,[1 2]),'UniformOutput',false)));
    
% else
%     
%     % dimensional correction when empty
%     x=zeros(3,0);
    
end

% decompose shape
if nargout>1
    
    % initiate variables
%     ax=zeros(size(x));
%     ang=zeros(size(x));
    gind=find(g_val);
    for i=1:size(H,3)
        
        % decompose
        [R,D]=schur(H(:,:,i)/2); %[R,D,~]=svd(H(:,:,i)); % we dont care about V flipping axis, we want to vis. a neighborhood
        d=s(:,:,i)+1/2*g(:,:,i)'*x(:,gind(i));
        D=-D/d;%-
        [D,j]=sort(diag(D),'descend'); % sort
        R=R(:,j);
        
        % write
        ang(:,gind(i))=rotm2eula(R,'ZXY')'; % pitch yaw roll
        ax(:,gind(i))=sqrt(1./D);
    end
    
%     % dimensional correction when empty
%     if isempty(ax)
%         ax=zeros(3,0);
%     end
%     if isempty(ang)
%         ang=zeros(3,0);
%     end
    
end

end

% % positive definite: better do this before data into function
% P=zeros(1,size(H,3));
% for i=1:size(H,3)
%     [~,P(:,i)]=chol(H(:,:,i));
% end
% H=H(:,:,P==0);
% g=g(:,:,P==0);
% s=s(:,:,P==0);
% % midpoint triangulation
% x=squeeze(cell2mat(cellfun(@(H,g)-H\g,num2cell(H,[1 2]),num2cell(g,[1 2]),'UniformOutput',false)));