function [ d ] = reprojdist( x,X, P )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xp=P*X;

xp=xp(1:2,:)./xp(3,:);
x=x(1:2,:)./x(3,:);

d=sqrt(sum((x-xp).^2,1));

end