function [ Dmap ] = def_dmap
%def_dmap Define distortion mapping
%   can be any mapping, but here the programs are written for generic
%   polynomial one. (optics mappings can be complex)

% get globals
global ctrl folder date cal

% define distortion mapping
x=sym('x',[2 1]);
[oy,ox]=meshgrid(0:ctrl.mord,0:ctrl.mord);
T=(ox+oy)<=ctrl.mord & (ox+oy)>1;
ox=ox(T); 
oy=oy(T);
T=reshape(T,[],1);

% coefficients b for selected terms
b=reshape(sym('b',[6+2*nnz(T) 1]),2,[]);

% define mapping
Dmap=b(:,1)+b(:,2:3)*x+b(:,4:end)*((x(1).^ox).*(x(2).^oy)); % i.e. first + h.o.t.

% export matlab function
Dmap=matlabFunction(Dmap,'Optimize',false);

save([folder date cal vsl 'Dmap.mat'],'Dmap')

end