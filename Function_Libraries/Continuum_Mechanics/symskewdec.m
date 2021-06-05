function [S,A] = symskewdec(M)
%symskewdec Symmetric [Anti/Skew]-symmetric Decomposition of a matrix
%   
%   Input
%       M Matrix, e.g. velocity gradient, data like M(i,j,:)
%   
%   Output
%       S Symmetric Part
%       A [Anti/Skew]-symmetric Part

% transpose
MT=permute(M,[2 1 3]); % use permute for multidimensional data input, see above

% symmetic part
S=(M+MT)/2;

% skew symetric part
A=(M-MT)/2;

end

