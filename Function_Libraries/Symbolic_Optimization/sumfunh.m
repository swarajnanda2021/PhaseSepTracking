function [ fun ] = sumfunh( fun1,fun2 )
%UNTITLED Summary of this function goes here
%   Sum function handles

% symbolic
fun1=sym(fun1);
fun2=sym(fun2);

fun=fun1+fun2;


fun=matlabFunction(fun,'Optimize',false);

end

