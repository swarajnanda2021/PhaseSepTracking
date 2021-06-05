function [ E ] = fmat2emat( F, K1,K2 )
%fmat2emat Fundamental matrix to essential matrix
%   input
%       Fundamental matrix
%       Camera matrices

E=K2'*F*K1;

end