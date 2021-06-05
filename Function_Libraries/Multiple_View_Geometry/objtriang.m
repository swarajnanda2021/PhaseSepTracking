function [X,k,r] = objtriang( x, R, t ) %,H
%objtriang Linear least squares solution for optimal object triangulation
%   this code minimize the point line distance in object space by and OLS
%   min_{X^c,k^c}sum_c||x^c*k^c-P^c*[X;1]||_2^2 for N multiple matched 
%   coordinates
%   
%   Input for each of C views:
%       x dewarped image coordinates at reference plane x=K\x_im
%           (3xN matrix )xC cell or 3xCxN matrix
%       R rotation matrices (3x3)xC cell
%       t translation vectors (3x1)xC cell
%   
%   Output for each paired coordinate set:
%       X object coordinates 3*N matrix
%       k depth of fields for each view C*N scalings
%       r object residuals for each triangulation
%       supressed H sensitivity from hessian matrices (3+C)x(3+C)xN (TO BE CHECKED)
%   
%   remark: nan values will be ignored in triangulation
%   
%   note: I dont like cell input that much, neither preformatted input,
%       can't better make an indexed but unstructured data input ?

% correct input
if nargin<3 % then R is actually -> P=[R t]
    t=cell(size(R)); % split the translation vector
    for i=1:length(R)
        t{i}=R{i}(:,4);
        R{i}=R{i}(:,1:3);
    end
end
if iscell(x)
    x=permute(cat(3,x{:}),[1 3 2]);
end
if size(x,1)==2
    x=inhc2homc(x);
end

% data size
C=size(x,2); % number of cameras
N=size(x,3); % number of objects

% mark valid data
val=~isnan(x); % valid data
val=(sum(val,2)>1).*val; % remove any single coordinates

% assemble least squares problem
t=sparse(repmat(-cat(1,t{:}),N,1));
R=kron(speye(N),sparse(-cat(1,R{:}))); % assembly camera translation vectors
x=reshape(x,3,1,N*C);
[r,c,b]=size(x);
i=repmat((1:r)',1,c)+(r*(reshape(1:b,1,1,[])-1));
j=repmat((1:c),r,1)+(c*(reshape(1:b,1,1,[])-1));
x=sparse(i(:),j(:),x(:),r*b,c*b); % max(i(:)),max(j(:)));

% block nan values
v=reshape(val,3,1,N*C);
[r,c,b]=size(v);
i=repmat((1:r)',1,c)+(r*(reshape(1:b,1,1,[])-1));
j=repmat((1:c),r,1)+(c*(reshape(1:b,1,1,[])-1));
v=sparse(i(:),j(:),v(:),r*b,c*b); % max(i(:)),max(j(:)));
vk=sum(v,1)>0; % valid k in solution vector
vd=sum(v,2)'>0; % valid data in assembly
vX=vd*R~=0; % valid X in solution vector

% solve linear equations
A=cat(2,R(vd,vX),x(vd,vk));
B=A'*A;
b=A'*t(vd);
y=-B\b; % solution vector

% write valid solution vector
sol=nan*ones(3*N+C*N,1); % intiate nans
sol(cat(2,vX,vk))=y;

% write solution vectors
X=full(reshape(sol(1:3*N),3,[])); % object coordinates
k=full(reshape(sol(3*N+1:end),C,[])); % depth of fields

% compute residuals
if nargout>=3
    d=sqrt(nansum(reshape((t+cat(2,R,x)*sol).^2,3,C,[]),1)); % residual distances
    r=reshape(nansum(d,2)./nansum(nanmean(val,1),2),1,N); % average distances
end

% compute sensitivity
% if nargout>=4
%     H=sparse(length(vX)+length(vk),length(vX)+length(vk)); % sensitivity ellipse
%     H(cat(2,vX,vk),cat(2,vX,vk))=B; % slow
%     num=reshape([reshape(1:3*N,3,N);reshape(3*N+(1:C*N),C,[])],1,[]);
%     H=H(num,num); % slow
%     H=reshape(full(H(kron(speye(N),sparse(ones(3+C)))==1)),3+C,3+C,N); % slow
% end

end

