function [ Pmap ] = def_pmap
%def_pmap Define projection mapping pinhole camera
%   
%   Input: processing setting
%
%	Output:
%       Pmap.mat projection mapping as symbolic matlab function

%% get globals
global folder date cal 

%% Define mapping
% hom. positions
X=sym('X',[4,1]); 

% camera projection matrix
P=sym('P',[3 4]); 

% projected coordinates
xp=P*X;

% projection mapping
Pmap=xp(1:2)./xp(3);

% export matlab function
Pmap=matlabFunction(Pmap,'Optimize',false);

%% save mapping function
save([folder date cal vsl 'Pmap.mat'],'Pmap')% save mapping

end