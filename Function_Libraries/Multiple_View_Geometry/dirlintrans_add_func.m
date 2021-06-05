function [ H ] = dirlintrans_add_func( x1 , x2 , C)
%hom_graph compute camera homography of a set of points on a plane
%   compact linear independent comp.

% if not, convert to cell
if ~iscell(x1)
    x1=num2cell(x1,[1 2]);
end
if ~iscell(x2)
    x2=num2cell(x2,[1 2]);
end

% data normalization
T1=cell(size(x1));
for i=1:length(x1)
    [T1{i},x1{i}]=datnorm(x1{i});
end
T2=cell(size(x2));
for i=1:length(x2)
    [T2{i},x2{i}]=datnorm(x2{i});
end

% compose dlt problem to solve multiple problems simulataneously 
A=cell(length(x1),length(x2));
for i=1:length(x1)
    for j=1:length(x2)
        % formulate coefficient matrix
        A{i,j}=sparse(...
            [zeros(size(x2{j}')) repmat(-x2{j}(3,:)',1,size(x2{j},1)).*x1{i}' repmat(x2{j}(2,:)',1,size(x2{j},1)).*x1{i}'
            repmat(x2{j}(3,:)',1,size(x2{j},1)).*x1{i}' zeros(size(x2{j}')) repmat(-x2{j}(1,:)',1,size(x2{j},1)).*x1{i}'
            repmat(-x2{j}(2,:)',1,size(x2{j},1)).*x1{i}' repmat(x2{j}(1,:)',1,size(x2{j},1)).*x1{i}' zeros(size(x2{j}'))]);
        
    end
end
A=cat(2,A{:});

% solve any constraints with prescribed 
if nargin==3 % basic constrain Cx=0
    % add rows if needed
    if size(C,1)<size(C,2)
        C=[C
            zeros(size(C,2)-size(C,1),size(C,2))];
    end
    
    % compute svd
    [~,~,V] = svds(C,1,'smallest');
    
    A=A*C;
    
% elseif nargin>=3
end

% solve the total l.s. problem on alg. dist, DLT solution
[~,~,V] = svds(A,1,'smallest');

% solve h
h=V(:,end);

% convert data to separate H
H=cell(length(x1),length(x2));
h=reshape(h,length(x1),length(x2),[]);
for i=1:length(x1)
    for j=1:length(x2)
        % select data
        hdat=squeeze(h(i,j,:));
        
        % in matrix form
        H{i,j}=reshape(hdat,size(x2{j},1),size(x1{i},1))';
        
        % data denormalization
        H{i,j}=T2{j}\H{i,j}*T1{i};
        
        % fix sign scale 
        H{i,j}=H{i,j}/sign(H{i,j}(end,end));
    end
end

end