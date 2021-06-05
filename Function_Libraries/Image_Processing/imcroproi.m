function [ I ] = imcroproi( I,roi )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

I=I(roi(1):roi(3),roi(2):roi(4),:);

end

