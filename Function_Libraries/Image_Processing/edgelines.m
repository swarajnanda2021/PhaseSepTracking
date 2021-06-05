function [lab,li,yi] = edgelines( imd, fobj ) % , typ, foi )
%edgelines Edgelines output and marks the edge data in an image in terms of
%   connected regions and corresponding point and line data
%   
%   Input:
%       imd Image data
%       fobj Filter object
%       
%   Output:
%       labE: labelled edges
%       li: pixel edgeline
%       yi: closest pixel point on line
%       
%   Remarks: ...
%       typ type of method (gradient | hessian.. | autocorrolation.. to be implemented)-tracing

% define filter, center filter, order, and window
switch class(fobj)
    case 'double'
        win=fobj;
        fobj=strel('arbitrary',ones(win));
    case 'strel'
        win=size(fobj.Neighborhood);
end
if length(win)<3
    win(3)=1;
end
ord=[1 1 win(3)>1]; % filter order
hwin=floor(win/2); % half flank filter window

% correct input
% if nargin<3
%     typ='undirected';
% end
if nargin<4
    foi=ceil(size(imd,3)/2); % define standard frame of interest
end

%-- #) savintsky golay image filtering

% define an image grid
[y,x]=meshgrid(1:size(imd,2),1:size(imd,1));

% find polynomial coefficients: based on S-G filter
coef=SGfilt_coef(imd,fobj,ord,1); % increasing cross order: shape deformation..

% attenuation coefficient
c0=coef{1,1,1}(:,:,foi);

% image gradient coefficient
cx=coef{2,1,1}(:,:,foi);
cy=coef{1,2,1}(:,:,foi);

% directional gradient
dg=sqrt(cx.^2 + cy.^2);

%figure; surf(x',y',dg'); axis tight; axis equal; view(2); shading flat

% point line distance
[pld,X]=pntlindist(inhc2homc(0*[x(:) y(:)]'),[cx(:) cy(:) c0(:)]','vector');
Y=reshape(X(2,:),size(imd,1),size(imd,2),size(imd,3));
X=reshape(X(1,:),size(imd,1),size(imd,2),size(imd,3));
pld=reshape(pld,size(imd,1),size(imd,2),size(imd,3));

%figure; surf(x',y',X'); axis tight; axis equal; view(2); shading flat
%figure; surf(x',y',pld'); axis tight; axis equal; view(2); shading flat

%-- #) Segmentation

% subpixel point line distance, for sharp line definition
% switch typ
%     case 'positive'
%         W= pld >=  0 & pld <= 1 ;  %
%     case 'negative'
%         W= pld >= -1 & pld <= 0 ;  %
%     case 'directed'
%         W= pld >= -1 & pld <= 1 ;  % half flank window /(norm(hwin,2)/sqrt(2)/2)
%     case 'undirected'
        W= pld >= -1 & pld <= 1 ;  % half flank window /(norm(hwin,2)/sqrt(2)/2)
% end

%figure; surf(x',y',double(W)'); axis tight; axis equal; view(2); shading flat

% part of prominent edge information, to remove weak edges
C=imsegscl(dg,0,1,4); % there is no definite setting

%figure; surf(x',y',double(C)'); axis tight; axis equal; view(2); shading flat

% feasible edge segmentation
F=W & C; % feasible set & msk

%figure; surf(x',y',double(F)'); axis tight; axis equal; view(2); shading flat

% find pixel indexing feasible segmentation
[i,j]=find(F);

% get linear indexing
k=sub2ind(size(F),i,j);

% define unique labeling feasible edges
labF=double(F);
labF(F)=1:length(k);

%figure; surf(x',y',labF'); view(2); axis tight; axis equal; shading flat

%-- #) List segmented edges

% transform lines to global grid
c0=c0-cx.*x-cy.*y;

% transform points on lines to global grid
X=x+X;
Y=y+Y;

% pixel coordinates edges
yi=cat(2,X(F),Y(F))'; % pixel location vector

%figure; plot(yi(1,:),yi(2,:),'.'); axis tight; axis equal; 

% % gradient vector corresponding to edges
% gi=cat(2,cx(F),cy(F))'; % image gradient vector

% tangent lines corresponding to edges
li=cat(2,cx(F),cy(F),c0(F))'; % homogenous line vector

%-- #) Define pixel adjacency 

% define neighbors within image window
[fi,fj]=find(fobj.Neighborhood); 
w=cat(2,fi-(hwin(1)+1),fj-(hwin(2)+1))';

% compute indexing neighboring pixels 
in=reshape(i-w(1,:),[],1); % i direction
jn=reshape(j-w(2,:),[],1); % j direction

% define valid neighbors within feasible segmentation
val=ismember([in jn],[i j],'rows');
in=in(val);
jn=jn(val);
kn=sub2ind(size(F),in,jn);

% get labeling feasible edges in valid connectivity
[lab,~]=find(reshape(val,[],size(w,2)));

% get labeling of the connected and corresponding neighboring edges
labn=labF(kn);

% compute angle and angle segmentation
alp=360./nnz(bwmorph(fobj.Neighborhood,'remove'));
aij=real(acosd(dot(li(1:2,lab)./vecnorm(li(1:2,lab),2,1),li(1:2,labn)./vecnorm(li(1:2,labn),2,1),1))); % tangent line direction pntlindist(inhc2homc(w),li,'matrix')' dtl=dtl(val);

%figure; plot(aij,'.')

% define pixel adjacency 
Adj=sparse(labn,lab,aij<alp,length(k),length(k));

Adj=Adj'*Adj>=(norm(hwin)/sqrt(2));
%figure; spy(Adj)

%-- #) Find connected components using the graph corrsponding to the adjacency

% define graph from edges
G=graph(Adj);

% find connected componenet from graph
lab=conncomp(G);

% labE=zeros(size(F));
% labE(F)=B;%double(C)*(1:size(C,2))';
% 
% %figure; surf(x',y',labE'); view(2); axis tight; colormap lines; axis equal; shading flat;

% %-- #) write output
% 
% % 
% !map=;
% !li=;
% !etc..

end

% % vote candidates - where is my most probable feature?
% [Vmap,~,~,I,J]=histcounts2(x(F)+X(F),y(F)+Y(F),1/2+(0:size(imd,1)),1/2+(0:size(imd,2))); % keep discrete, do not smooth
% 
% %figure; surf(x',y',Vmap'); axis tight; axis equal; view(2); shading flat
% 
% % linear indexing F 
% indF=find(F);

% w=[ -1  0  1 1 1 0 -1 -1
%     -1 -1 -1 0 1 1  1  0 ]; % 8 connected for now

% % compute pointline distance gradient direction
% dgl=reshape(pntlindist(inhc2homc(w),rotz(90)*inhc2homc(gi,0),'matrix')',[],1); % gradient line (centered at pixel) direction
% dgl=dgl(val);
% Dgl=sparse(labn,lab,abs(dgl)<=1,length(k),length(k));
% 
% %figure; plot(dgl,'.')
% %figure; spy(Dgl)

% aij=real(acosd(dot(li(1:2,lab)./vecnorm(li(1:2,lab),2,1),li(1:2,labn)./vecnorm(li(1:2,labn),2,1),1))); % tangent line direction pntlindist(inhc2homc(w),li,'matrix')' dtl=dtl(val);
% % aji=real(acosd(dot(-li(1:2,lab)./vecnorm(li(1:2,lab),2,1),li(1:2,labn)./vecnorm(li(1:2,labn),2,1),1))); % tangent line direction pntlindist(inhc2homc(w),li,'matrix')' dtl=dtl(val);
% % dji=reshape(pntlindist(inhc2homc(yi(:,lab)),li(:,labn),'vector')'; % tangent line direction pntlindist(inhc2homc(w),li,'matrix')' dtl=dtl(val);
% % dtl=dij-dji;
% Dtl=sparse(labn,lab,aij<5,length(k),length(k)); %  | aji<15 abs(dtl)<=sqrt(2)
% 
% %figure; plot(aij,'.')
% %figure; spy(Dtl)
% 
% % define pixel adjacency 
% Adj=Dtl & Dtl';%'*Dtl >=5;%( ( ( Dtl & Dtl' ) | ( Dgl & Dgl' ) ) - ( ( Dtl & Dgl' ) | ( Dtl' & Dgl ) ) ) > 0.5; % connect gradient direction or
% 
% %figure; spy(Adj)


% compute pointline distance tangent line neighboring pixels
% dij=reshape(pntlindist(inhc2homc(yi(:,labn)),li(:,lab),'vector')',[],1); % tangent line direction pntlindist(inhc2homc(w),li,'matrix')' dtl=dtl(val);
% dji=reshape(pntlindist(inhc2homc(yi(:,lab)),li(:,labn),'vector')',[],1); % tangent line direction pntlindist(inhc2homc(w),li,'matrix')' dtl=dtl(val);
% % dtl=dij-dji;
% Dtl=sparse(labn,lab,dij<=0.5 & dji<=0.5,length(k),length(k)); % abs(dtl)<=sqrt(2)

% % self adjacency
% Adj=Adj | speye(size(Adj));


% %-- #) Branch pixel neighbors
% 
% B=speye(size(Adj));
% 
% B=double(Adj+speye(size(Adj)))*B>0.5;
% 
% %figure; spy(B)
% 
% [~,dum]=max(double(B>0.5)*spdiags(fliplr(1:size(B,2))',0,size(B,1),size(B,2)),[],2);
% dum=unique(dum);
% 
% %figure; plot(dum,'.')
% 
% C=B(:,dum);
% 
% %figure; spy(C)
% %figure; plot(sum(C,1),'.')

% bwmorph(,'remove') bwmorph(,'dilate')
% lobj=strel('disk',8,0);


% % directional gradient magnitude edges
% dgi=dg(F);
% 
% %figure; scatter(yi(1,:),yi(2,:),[],dgi,'.'); axis tight; axis equal; 
% 
% % gradient directed point line distance edges
% pldi=pld(F);
% 
% %figure; scatter(yi(1,:),yi(2,:),[],pldi,'.'); axis tight; axis equal; 
