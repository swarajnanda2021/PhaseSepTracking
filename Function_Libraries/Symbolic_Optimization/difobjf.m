function [ dObjF , ddObjF ] = difobjf( ObjF , dvars )
%difobjf Differentiate the objective function
%   exports the local gradient and hessian of the roblem to be asseembled

% convert to symbolic
ObjF=sym(ObjF);

% fix input variables w.r.t. consistency reg. objective func.
vars=sort(symvar(ObjF));
if nargin<2
    dvars=vars; % integrate w.r.t. all variables
end

if length(ObjF)==1 % scalar
    % gradient
    dObjF=gradient(ObjF,dvars);
    
    if nargout==2
        % hessian
        ddObjF=triu(hessian(ObjF,dvars)); % only store upper triangle (speed)
        
    end
else % vectorial
    % jacobian
    dObjF=jacobian(ObjF,dvars);
    
    if nargout==2
        % second der def
        if size(dObjF,1)==size(dObjF,2)
            ddObjF=det(dObjF);
            
        elseif length(dvars)==1
            ddObjF=jacobian(dObjF,dvars);
            
        end
        
    end
end

end