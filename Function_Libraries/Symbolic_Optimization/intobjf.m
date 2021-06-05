function [ iObjF ] = intobjf( ObjF , ivars, ibounds )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% convert to symbolic
Obj=sym(ObjF);

% fix input variables
vars=sort(symvar(Obj));
if nargin<2
    ivars=vars; % integrate w.r.t. all variables
end

% integrate
if nargin<3
    iObj=Obj;
    for k=1:length(ivars)
        iObj=int(iObj,ivars(k));
    end
else
    iObj=Obj;
    for k=1:length(ivars)
        iObj=int(iObj,ivars(k),ibounds(k,1),ibounds(k,2));
    end
    % define new set of vars
    vars=sort(symvar(iObj)); % n vars may change w.r.t. objf
end

% variables
vars=sym2cell(vars);

% write gradient
iObjF=matlabFunction(iObj,'Optimize',false,'Vars',vars);

end

