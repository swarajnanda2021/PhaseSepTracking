function [ Pmap ] = def_pmap
%def_pmap Define projective mapping
%   this is a general mapping using the projection matrix [x/k y/k -]

% get globals
global ctrl folder date cal

% hom. projective positions
X=sym('X',[4,1]); 

% camera projection matrix
P=sym('P',[3 4]); 

% projected coordinates
xp=P*X;

% projectio mapping
Pmap=xp(1:2)./xp(3);

% export matlab function
Pmap=matlabFunction(Pmap,'Optimize',false);

save([folder date cal vsl 'Pmap.mat'],'Pmap')

end