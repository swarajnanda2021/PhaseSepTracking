function [ T,xn ] = datnorm( x )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dim=size(x,1); % dimensionality

t=-mean(x,2);
t(end)=1;

s=sqrt(dim-1)./std(x,[],2); % unit circle / sphere
s(end)=1;

Tt=eye(dim);

Tt(:,end)=t;

Ts=diag(s);

T=Ts*Tt;

xn=T*x;

end