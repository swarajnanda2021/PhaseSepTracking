function [ x ] = homc2inhc( x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

siz=size(x);

x=reshape(x,siz(1),[]);

x=x(1:end-1,:)./x(end,:);

siz(1)=siz(1)-1;

x=reshape(x,siz); % preserves shape / dimensionality input

end