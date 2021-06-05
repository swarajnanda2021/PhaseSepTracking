function [ cwin, filt, pset ] = genellipse(wfun , fres , rad , ang , asp )
%genellipse generate ellipse

% window size
cwin=[-ceil(sqrt(2)*max(rad)*max(asp)) ceil(sqrt(2)*max(rad)*max(asp))];

% symbolic functional to generate level-set
x=sym('x',[2 1]);
thet=sym('thet');
eps=sym('eps');
R=[cos(thet) -sin(thet)
    sin(thet) cos(thet)];
D=[1 0
    0 eps];
x=R\D\x;
Lset=matlabFunction(sqrt(dot(x,x)),'Optimize',false);

% meshgrid for level-set
[y,x]=meshgrid(cwin(1):cwin(2),cwin(1):cwin(2));

% generate parameter set without degeneracy
[rho,thet,eps]=meshgrid(rad,ang,asp); % generate all combinations
rho=reshape(rho,[],1);
thet=reshape(thet,[],1);
eps=reshape(eps,[],1);
sel=eps(:)==1;
pset=unique([eps(~sel) thet(~sel) rho(~sel)
    eps(sel) 0*thet(sel) rho(sel)],'rows'); % list them

% initiate
filt.shap=zeros(size(y,1),size(y,2),size(pset,1));
% filt.edge=zeros(size(y,1),size(y,2),size(pset,1));
% filt.norm=zeros(size(y,1),size(y,2),size(pset,1));
figure
for i=1:size(pset,1)
    % radius
    r=reshape(Lset(pset(i,1),pset(i,2),x(:),y(:)),size(x));
    
    % convolution filters
    switch wfun
        case 'peak'
%             filt.shap=1/2*erfc(1/(sqrt(2)*fres)*(r-pset(i,3)));
%             filt.edge=exp(-1/(2*fres^2)*(r-pset(i,3)).^2);
%             filt.norm=double(r<sqrt(2)*pset(i,3));
        case 'edge'%;%erfc(1/(sqrt(2)*fres)*(r-pset(i,3)))
            filt.shap(:,:,i)=(2*double(r<pset(i,3))-double(r<sqrt(2)*pset(i,3)))/nnz(r<pset(i,3));
%             filt.edge(:,:,i)=exp(-(r-pset(i,3)).^2/(2*fres^2));
%             filt.norm(:,:,i)=double(r<sqrt(2)*pset(i,3));
    end
    
    % plotting
%     subplot(1,2,1)
    surf(x',y',filt.shap(:,:,i)')
    view(2)
    axis tight
    axis square
    shading flat
    drawnow
    
%     subplot(1,2,2)
%     surf(x',y',filt.edge(:,:,i)')
%     view(2)
%     axis tight
%     axis square
%     shading flat
%     drawnow
    
%     subplot(1,3,3)
%     surf(x',y',filt.norm(:,:,i)')
%     view(2)
%     axis tight
%     axis square
%     shading flat
%     drawnow
    
    pause
    
end

% close all

end