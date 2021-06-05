function [ x,y,X,Y ] = croproigrid( x,y,X,Y,roi )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% image grid at reference
x=x(roi(1):roi(3),roi(2):roi(4));
y=y(roi(1):roi(3),roi(2):roi(4));

% deformed grid
X=X(roi(1):roi(3),roi(2):roi(4));
Y=Y(roi(1):roi(3),roi(2):roi(4));

end

