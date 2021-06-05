function x = projcoord(X,Pmat)
%projcoord Project world coordinate
%   
%   Input
%       X World coordinate
%       Pmat projection matrix
%   
%   Output
%       x projected coordinate
%   

% correct input
if size(X,1)==3
    X=inhc2homc(X);
end

% transform reference
X=Pmat*X;

% project and write 
x=homc2inhc(X);

end

