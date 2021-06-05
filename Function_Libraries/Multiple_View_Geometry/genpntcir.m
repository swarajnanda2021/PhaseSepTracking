function [ x ] = genpntcir( N )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

theta=linspace(0,2*pi-2*pi/N,N);%pi/N,2*pi-pi/N
rho=ones(size(theta));
[x,y]=pol2cart(theta,rho);

x=[x;y];

end

