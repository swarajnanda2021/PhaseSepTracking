function [ Cfit  ] = def_cfit
%def_cfit Define curve fit
%   
%   Input: processing setting
%
%	Output:
%       Cfit.mat defined curve as symbolic matlab function

%% get globals
global ctrl folder date cal 

%% define curve
% symbolic expression path
syms t
ot=(0:ctrl.cord)'; % order coef
a=reshape(sym('a',[2*(ctrl.cord+1) 1]),[],2);

% define curve
Cfit=a'*(t.^ot); % i.e.centered x by mapping + all h.o.t.

% export matlab function
Cfit=matlabFunction(Cfit,'Optimize',false);

%% save parametrized curve
save([folder date cal vsl 'Cfit.mat'],'Cfit')% save mapping

end