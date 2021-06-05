function [ val ] = evalsfun( fun , varn )
% slow

vars=symvar(fun);

val=double(subs(fun,vars,varn));

end