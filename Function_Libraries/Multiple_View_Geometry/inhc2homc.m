function [ x ] = inhc2homc( x , k )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    k=1;
end

siz=size(x);

x=reshape(x,siz(1),[]);

x=[ x
    k*ones(1,size(x,2)) ];

siz(1)=siz(1)+1;

x=reshape(x,siz); % preserves shape / dimensionality input

end

