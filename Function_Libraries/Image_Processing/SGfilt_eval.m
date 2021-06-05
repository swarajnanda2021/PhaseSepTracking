function [ val ] = SGfilt_eval( coef, der )
%eval_coef_SG Evaluate coefficients Savintzky Golay filter and derivatives
%   
%   $Author: Koen Muller$
%   
% Input,
%   coef: coefficients from filt_coef_SG
%   der: desired derivative e.g. [1 0 0]->'dx' | [1 2 0]->'dxddy'
%
% Output,
%   val: nummeric evaluation

% 0) get size coefficients
% siz=size(coef);

% 1) get derivative info
w=prod(factorial(der));
c=coef{der(1)+1,der(2)+1,der(3)+1};

% 2) compute derivative
val=w*c;

end