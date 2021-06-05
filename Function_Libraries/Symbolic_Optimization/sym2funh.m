function [ funh ] = sym2funh( func , vars)
%try-out, faster to evaluate this function but slower on execution than
%    matlabFunction

% function variables
if nargin<2
    vars=sort(symvar(func));
end

% write anonymous argument
args = sprintf('%s,', vars);
args(end)=[];

% write anonymous expression
expr=[];
for i=1:size(func,1)
%     disp(i)
    for j=1:size(func,2)
%         disp(j)
        expr = sprintf('%s%s,',expr, func(i,j));
    end
    expr(end)=';';
end
expr(end)=[];

% write output
% tic
funh = eval(sprintf('@(%s)[%s]', args, expr));
% toc

end

