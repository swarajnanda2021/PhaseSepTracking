function c = projquad(q,Pmat)
%projquad Project world quadric
%   
%   Input
%       Q World coordinate
%       Pmat projection matrix
%   
%   Output
%       c projected conic
%   

% matrix form
Q=qvec2qmat(q);

% project
C=zeros(3,3,size(q,2));
C(:,:,:)=cell2mat(cellfun(@(x)inv(Pmat*(x\Pmat'))...
    ,num2cell(Q,[1 2]),'UniformOutput',false));

% write
c=cmat2cvec(C);

end

