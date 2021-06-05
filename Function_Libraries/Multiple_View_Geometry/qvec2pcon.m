function [ p , C ] = qvec2pcon( q , c )
%qvec2pcon quadric to plane conic described by origin

% quadric matrix
Q=qvec2qmat(q);

% replicate c
c=repmat(c,1,1,size(Q,3));

% find plane
p=cell2mat(cellfun(@(Q,c)Q*c,...
    num2cell(Q,[1 2]),num2cell(c,[1 2]),'UniformOutput',false));

% find coordinate system
[ M ] = pleq2pmat( p );

% find conic
C=cell2mat(cellfun(@(M,Q)cmat2cvec(M'*Q*M),...
    num2cell(M,[1 2]),num2cell(Q,[1 2]),'UniformOutput',false));

% squeeze
p=squeeze(p);
C=squeeze(C);

end

