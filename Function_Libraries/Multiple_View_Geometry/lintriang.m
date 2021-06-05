function [ X ] = lintriang( x, P )
%lintriang solve linear triangulation to homogeneous coordinates
%   here we use the inhomogeneous solution from Hartley and Zisserman 2004
%   to be able to vectorize the normal equations instead of looping and svd
%   for the homogeneous equation, note that this means we cannot solve 
%   points at infinity.
%
%  	Note: this is the fastest triangulation method in terms of linearity
%         and number unknowns;

% correct input
P=reshape(P,[],1);
x=reshape(cellfun(@(x)reshape(x,3,1,[]),x,'UniformOutput',false),[],1);
P=cellfun(@(P,x)repmat(P,1,1,size(x,3)),P,x,'UniformOutput',false);

% make coefficient matrix A
A=cellfun(@(x,P) x([1 2 1],1,:).*P([3 3 2],:,:)...
    -x([3 3 2],1,:).*P([1 2 1],:,:),x,P,'UniformOutput',false);
A=cell2mat(A);

% split coefficient matrix for inhomogenous solution
a=A(:,4,:);
A=A(:,1:3,:);

% create sparse numbering
[r,c,b]=size(A);
i=repmat((1:r)',1,c)+(r*(reshape(1:b,1,1,[])-1));
j=repmat((1:c),r,1)+(c*(reshape(1:b,1,1,[])-1));
A=sparse(i(:),j(:),A(:),r+r*(b-1),c+c*(b-1));%max(i(:)),max(j(:)));

% create normal equation
M=A'*A;
b=A'*a(:);

% solve normal equations
X=-M\b;
X=inhc2homc(reshape(X,3,[]));

end

% x=reshape(num2cell(reshape(cat(1,x{:}),3,1,[]),[1 2]),length(P),[]);
% P=repmat(reshape(P,[],1),1,size(x,2));
% A=cellfun(@(x,P) x([1 2 1]).*P([3 3 2],:)-x([3 3 2]).*P([1 2 1],:),x,P,'UniformOutput',false);
% A=reshape(A,size(A,1),1,[]);
% A=num2cell(cell2mat(A),[1 2]);
% M=cellfun(@(A)sparse(A(:,1:3)'*A(:,1:3)),A,'UniformOutput',false);
% b=cellfun(@(A)sparse(A(:,1:3)'*A(:,4)),A,'UniformOutput',false);
% 
% M=blkdiag(M{:});
% b=cat(1,b{:});

% % solve svd problem
% tic
% X=cellfun(@(A)leftsinvec(A),A,'UniformOutput',false);
% X=cat(2,X{:});
% toc

% % % solve normal equation
% X=cellfun(@(A)-inv(A(:,1:3)'*A(:,1:3))*A(:,1:3)'*A(:,4),A,'UniformOutput',false);
% X=cat(2,X{:});

% 
% % solve svd problem
% tic
% X=cellfun(@(A)leftsinvec(A),A,'UniformOutput',false);
% X=cat(2,X{:});
% toc
% 
% tic
% % coefficient matrix for svd
% X=zeros(4,size(x{1},2));
% for i=1:size(x{1},2)
%     A=cell(size(x));
%     for j=1:length(x)
%         % coefficient matrix
%         A{j}=[x{j}(1,i)'*P{j}(3,:) - x{j}(3,i)'*P{j}(1,:)
%             x{j}(2,i)'*P{j}(3,:) - x{j}(3,i)'*P{j}(2,:)
%             x{j}(1,i)'*P{j}(2,:) - x{j}(2,i)'*P{j}(1,:)];
%         % vectorize j
%     end
%     A=vertcat(A{:});
%     
%     [~,~,V]=svd(A);
%     
%     X(:,i)=V(:,end)/(V(end,end)); % select more solution in vec
%     
%     % speed up svd only by coding direct solution and vectorization to that
% end
% toc

% function [ v ] = leftsinvec( X )
% %UNTITLED3 Summary of this function goes here
% %   Detailed explanation goes here
% 
% [~,~,V]=svd(X);
% v=V(:,end)/V(end,end);
% 
% end