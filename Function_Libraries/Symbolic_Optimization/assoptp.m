function [ GlbF, gGlbF, hGlbF ] = assoptp( x , dat , wObjF , ObjF , gObjF , hObjF )
%assoptp Assemble optimization problem,
%
%   Input
%       symbolic function
%
%   Output
%       numeric evaluation
%
%   Note: Could also estimate a tau and amp(?) from hessian connectivity.

% corrct input by nargin for empty test
if nargin<5
    gObjF=[];
end
if nargin <6
    hObjF=[];
end

% numbering purpose
V=logical(dat(:,1));
num=dat(:,2);

% initiate
if ~isempty(ObjF)
    GlbF=0;
else
    GlbF=[];
end
if ~isempty(gObjF)
    gGlbF=zeros(length(x),1);
else
    gGlbF=[];
end
if ~isempty(hObjF)
    hGlbF=cell(size(dat,2)-2,1); % loop for indexed sparse array
else
    hGlbF=[];
end
for i=3:size(dat,2)
    % indexing
    ind=dat(V,i);
    
    % evaluate local objective
    inp=zeros(size(num));
    inp(num)=[dat(~V,i)
        x(dat(V,i))];
    inp=num2cell(inp'); % (num)
    
    % compute weight
    if ~isempty(wObjF)
        wobj=wObjF(inp{:});
    else
        wobj=1;
    end
    
    % objective value
    if ~isempty(ObjF)
        obj=ObjF(inp{:});
        GlbF=GlbF+obj*wobj;
    end
    
    % gradient vector
    if ~isempty(gObjF)
        gobj=gObjF(inp{:});
        gGlbF(ind)=gGlbF(dat(V,i))+gobj*wobj;
    end
    
    % hessian matrix
    if ~isempty(hObjF)
        hobj=(hObjF(inp{:}));
        [I,J,K]=find(hobj);
        hGlbF{i-2}=[ind(I),ind(J),K*wobj];
    end
%             if mod(i,10000)==0
%                 disp(num2str(i))
%             end
end
if ~isempty(hObjF) % if more nargin in future ?
    hGlbF=cell2mat(hGlbF);
    hGlbF=sparse(hGlbF(:,1),hGlbF(:,2),hGlbF(:,3),length(x),length(x)); % accumulate values
    hGlbF=hGlbF+tril(hGlbF',-1); % distribute the upper triangle
end

end