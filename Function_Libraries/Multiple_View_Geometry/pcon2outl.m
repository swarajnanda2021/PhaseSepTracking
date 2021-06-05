function [ X ] = pcon2outl( p, C, x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% generate points on cirkel if not as input
if nargin<3
    x=genpntcir(100);
end

% change dimensions
C=cvec2cmat(C);
p=reshape(p,4,1,[]);

% find coordinate system
[ M ] = pleq2pmat( p );

% evaluate ellipse
Xe=cell2mat(cellfun(@(C)cvec2pnts(cmat2cvec(C),x),...
    num2cell(C,[1 2]),'UniformOutput',false));

% invert back coordinate system
X=cell2mat(cellfun(@(M,X)homc2inhc(M*inhc2homc(X)),...
    num2cell(M,[1 2]),num2cell(Xe,[1 2]),'UniformOutput',false));

end

