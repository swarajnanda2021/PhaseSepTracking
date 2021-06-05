function [ fun ] = divfunh( fun1, fun2 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


% symbolic
fun1=sym(fun1);
fun2=sym(fun2);

fun=fun1./fun2;


fun=matlabFunction(fun,'Optimize',false);

end