function [ obj ] = subfunv( obj,old,new )
%subobjv Substitute objective variables
%   in order of obj func(arg(var))

obj=sym(obj);
old=sym(old);
new=sym(new);

obj=subs(obj,old,new);

vars=sort(symvar(obj));

obj=matlabFunction(obj,'Optimize',false,'Vars',vars);

end