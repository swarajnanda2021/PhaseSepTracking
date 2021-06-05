function [ t,R,s,d ] = procrust( x1,x2, typ ) %
%procrust Procrustes analysis alternative to matlab inbuilt
%   transformation sequence
%       translation
%       rotation
%       scaling
%
%   different scaling
%       best rigid body (no scaling)
%       best uniform scaling
%       best component scaling
%       best matrix scaling (when decomposed in R*D*R' R' rotate x1 and x2)
%
% speed boost here would be to vectorize the LSQ solutions in components

% correct input
if nargin<3
    typ='uniform'; % standard procrusters analysis
end

%figure; plot3(x2(1,:),x2(2,:),x2(3,:),'.'); hold on; plot3(x1(1,:),x1(2,:),x1(3,:),'.')

% compute midpoint data
x1m=mean(x1,2);
x2m=mean(x2,2);

% compute the translation vector
t=x2m-x1m;

%figure; plot3(x2(1,:),x2(2,:),x2(3,:),'.'); hold on; plot3(reshape(x1(1,:,:)+t(1,:,:),1,[]),reshape(x1(2,:,:)+t(2,:,:),1,[]),reshape(x1(3,:,:)+t(3,:,:),1,[]),'.')

% center data
x1n=x1-x1m;
x2n=x2-x2m;

%figure; plot3(x2n(1,:),x2n(2,:),x2n(3,:),'.'); hold on; plot3(x1n(1,:),x1n(2,:),x1n(3,:),'.')

% initiate variables
x1r=zeros(size(x1n));
x1s=zeros(size(x1n));
R=zeros(3,3,size(x1n,3));

% switch different scalings
switch typ
    case {'rigid','uniform'}
        s=zeros(1,1,size(x1n,3));
    case 'component'
        s=zeros(3,1,size(x1n,3));
    case 'matrix'
        s=zeros(3,3,size(x1n,3));
end

% loop to compute transformations
for i=1:size(x1n,3)
    
    %figure; plot3(x2n(1,:,i),x2n(2,:,i),x2n(3,:,i),'.'); hold on; plot3(x1n(1,:,i),x1n(2,:,i),x1n(3,:,i),'.')
    
    % compute the covariance matrix
    C=x2n(:,:,i)*x1n(:,:,i)';
    
    % compute svd
    [U,S,V]=svd(C);
    
    % constrain proper
    D=eye(size(S));
    D(end)=det(U*V');
    
    % define rotation matrix
    R(:,:,i)=U*D*V';
    
    % rotate coordinates
    x1r(:,:,i)=R(:,:,i)*x1n(:,:,i);
    
    %figure; plot3(x2n(1,:,i),x2n(2,:,i),x2n(3,:,i),'.'); hold on; plot3(x1r(1,:,i),x1r(2,:,i),x1r(3,:,i),'.')
    %figure; plot3(x2n(1,:),x2n(2,:),x2n(3,:),'.'); hold on; plot3(x1r(1,:),x1r(2,:,i),x1r(3,:),'.')
    
    % switch different scalings
    switch typ
        case 'rigid' 
            
            % rigid body motion: s=1
            s(:,:,i)=1;
            
            % scale coordinates
            x1s(:,:,i)=x1r(:,:,i);
            
        case 'uniform'
            s=(reshape(x1r(:,:,i)',[],1)'*reshape(x2n(:,:,i)',[],1))...
                /(reshape(x1r(:,:,i)',[],1)'*reshape(x1r(:,:,i)',[],1));
            
            % scale coordinates
            x1s(:,:,i)=s(:,:,i)*x1r(:,:,i);
            
        case 'component' 
            % compute component LSQ solution
            s(:,:,i)=blkdiag(x1r(1,:,i),x1r(2,:,i),x1r(3,:,i))*blkdiag(x1r(1,:,i)',x1r(2,:,i)',x1r(3,:,i)')...
                \(blkdiag(x1r(1,:,i),x1r(2,:,i),x1r(3,:,i))*reshape(x2n(:,:,i)',[],1));
            
            % scale coordinates
            x1s(:,:,i)=diag(s(:,:,i))*x1r(:,:,i);
            
        case 'matrix'
            % compute LSQ to coefficients
            s(:,:,i)=reshape(...
                blkdiag(x1r(:,:,i),x1r(:,:,i),x1r(:,:,i))*blkdiag(x1r(:,:,i)',x1r(:,:,i)',x1r(:,:,i)')...
                \(blkdiag(x1r(:,:,i),x1r(:,:,i),x1r(:,:,i))*reshape(x2n(:,:,i)',[],1)),...
                3,3)';
            
            % scale coordinates
            x1s(:,:,i)=s(:,:,i)*x1r(:,:,i);
            
    end
    
    %figure; plot3(x2n(1,:,i),x2n(2,:,i),x2n(3,:,i),'.'); hold on; plot3(x1s(1,:,i),x1s(2,:,i),x1s(3,:,i),'.')
    %figure; plot3(x2n(1,:),x2n(2,:),x2n(3,:),'.'); hold on; plot3(x1s(1,:),x1s(2,:),x1s(3,:),'.')
    
end

% Compute Procrustes dissimilarity distance
d=[mean(sqrt(sum((x2n-x1n).^2,1)),2) % translational dissimilarity
	mean(sqrt(sum((x2n-x1r).^2,1)),2) % rotational dissimilarity
	mean(sqrt(sum((x2n-x1s).^2,1)),2)]; % scaling dissimilarity

%figure; plot(d(1,:),'.')
%hold on; plot(d(2,:),'.')
%plot(d(3,:),'.'); hold off

end

% covariance matrix x1 x2
% C=[sum(x2n(1,:,:)*x1n(1,:,:),2) sum(x2n(1,:,:)*x1n(2,:,:),2) sum(x2n(1,:,:)*x1n(3,:,:),2)
%     sum(x2n(2,:,:)*x1n(1,:,:),2) sum(x2n(2,:,:)*x1n(2,:,:),2) sum(x2n(2,:,:)*x1n(3,:,:),2)
%     sum(x2n(3,:,:)*x1n(1,:,:),2) sum(x2n(3,:,:)*x1n(2,:,:),2) sum(x2n(3,:,:)*x1n(3,:,:),2)];

% define rotation matrix
%  R=[U(1,:,:)*D(1,:,:)     *V';

%     % decompose rq sequence
%     [S(:,:,i),R(:,:,i)]=rq(A);