function [ M ] = pleq2pmat( p )
%pleq2plnc plane equation to 4x3 plane parametrization matrix

% reshape p
p=reshape(p,4,1,[]);

% construct matrix M (non-unique) 
M=cell2mat(cellfun(@(x)[-[x(2)/x(1) x(3)/x(1) x(4)/x(1)]' eye(3)]',...
    num2cell(p,[1 2]),'UniformOutput',false));

end

