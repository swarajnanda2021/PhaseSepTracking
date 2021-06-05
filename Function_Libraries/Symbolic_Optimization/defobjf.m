function [ objF, z ] = defobjf( argF , eqcF , typ)
%defobjf Define objective function
%   this function defines the lagragian and can therefore handdle equality
%   constraints

% input
if nargin<2
    eqcF=[]; % empty constraints
end
if nargin<3
    typ='comb'; % favorite
end

% convert to symbolic expressions
argF=sym(argF);
eqcF=sym(eqcF);

% define lagrange multipliers
z=sym('z',size(eqcF)); % append at end using z small

% type of function
switch typ
    case 'lagr' % consistent
        conF=dot(z,eqcF); % can be non-positive // non convex
    case 'merit' % convex penalty function 
        conF=dot(eqcF,eqcF); % scaling problem // insignificant
    case 'comb' % convex and consistent - self invented - to be checked
        conF=dot(z+eqcF,z+eqcF);
end

% define objective function (obj) + (lag) + (pen)    ?? dot(z+eqcF,z+eqcF);%+
objF=dot(argF,argF)+conF; 

end