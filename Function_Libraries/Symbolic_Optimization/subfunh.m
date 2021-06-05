function [ fun ] = subfunh( fun1, fun2)
%UNTITLED2 Summary of this function goes here
%   substract function handless

% symbolic
fun1=sym(fun1);
fun2=sym(fun2);

fun=fun1-fun2;


fun=matlabFunction(fun,'Optimize',false);

end

