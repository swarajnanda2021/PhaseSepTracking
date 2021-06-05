function [ W ] = absconhom( H )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% formulate coefficient matrix (dyadic structure, improve automation)
B=cell(size(H));
for k=1:length(H)
    A=cell(2);
    for i=1:2
        for j=1:2
            A{i,j}=[H{k}(1,i)*H{k}(1,j)...
                H{k}(1,i)*H{k}(2,j)+H{k}(2,i)*H{k}(1,j)...
                H{k}(2,i)*H{k}(2,j)...
                H{k}(3,i)*H{k}(1,j)+H{k}(1,i)*H{k}(3,j)...
                H{k}(3,i)*H{k}(2,j)+H{k}(2,i)*H{k}(3,j)...
                H{k}(3,i)*H{k}(3,j)];
        end
    end
    B{k}=[A{1,2}%+A{2,1}
        A{1,1}-A{2,2}]; % first term is based on zero addition, no prob due scale
end
B=cell2mat(B');

% solve l.s. problem on alg. dist, DLT solution
[~,~,V] = svd(B);

% solve h
w=V(:,end);

% in matrix form
W= [w(1) w(2) w(4)
    w(2) w(3) w(5)
    w(4) w(5) w(6)];

W=W/W(end,end); % make sure that W is pos definite here using a scale argument

end

