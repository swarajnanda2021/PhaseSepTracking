function voightvec = symmat2voightvec(M)
%symmat2voightvec Symmetric 3x3 Matrix to Voight vector notation
%   
%   Input
%       M Symmetric (2x2) 3x3xN Matrix
%
%   Output
%       voightvec Voightvector notation 3xN

% correct input
if size(M,1)==2
    M=[M
        zeros(1,size(M,2))]; % zeropad 3
end
if size(M,2)==2
    M=[M zeros(size(M,1),1)]; % zeropad 3
end

% assemble voight vector
voightvec=[M(1,1,:)
        M(2,2,:)
        M(3,3,:)
        2*M(2,3,:)
        2*M(1,3,:)
        2*M(1,2,:)];


end

