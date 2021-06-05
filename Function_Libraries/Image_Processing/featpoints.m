function [ XX, UU, CC, II, WW ] = featpoints( imd, fobj , typ , foi ) % , msk
%featpoints Obtain [moving] feature points from an image
%   input,
%       imd: image data
%       fobj: centered window filtering object (double | strel)
%       typ: find levset minima | maxima | saddles | real ellipses @ conic
%   output,
%       midpoint location
%       velocity vector
%       conic shape
%       levelset intensity attenuation
%
% this function does not use intensity segmentation but function properties
% instead, instensity values segmentation can be used outside this function


std_param = 0.5;
iter_param = 7;


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
ord=2*[1 1 win(3)>1]; % filter order
hwin=floor(win/2); % half flank filter window

% correct input
if nargin<4
    foi=ceil(size(imd,3)/2); % define standard frame of interest
end
% if nargin<5
%     msk=ones(size(imd,1),size(imd,2)); % define standard frame of interest
% end

% define an image grid
[y,x]=meshgrid(1:size(imd,2),1:size(imd,1));

% find polynomial coefficients: based on S-G filter
coef=SGfilt_coef(imd,fobj,ord,1); % increasing cross order: shape deformation..

% attenuation coefficient
c0=coef{1,1,1}(:,:,foi);
if hwin(3)>0
    ct=coef{1,1,2}(:,:,foi);
    ctt=coef{1,1,3}(:,:,foi);
else
    ct=zeros(size(imd,1),size(imd,2));
    ctt=zeros(size(imd,1),size(imd,2));
end

% image gradient coefficient
cx=coef{2,1,1}(:,:,foi);
cy=coef{1,2,1}(:,:,foi);

% image hessian coefficients (sym. mat)
cxx=coef{3,1,1}(:,:,foi);
cxy=coef{2,2,1}(:,:,foi);
cyy=coef{1,3,1}(:,:,foi);

% shift coefficients
if hwin(3)>0
    cxt=coef{2,1,2}(:,:,foi);
    cyt=coef{1,2,2}(:,:,foi);
else
    cxt=zeros(size(imd,1),size(imd,2));
    cyt=zeros(size(imd,1),size(imd,2));
end

% compute determinant and trace hessian
trH=2*cxx+2*cyy; % trace
detH=4*cxx.*cyy-cxy.^2; % determinant
% detC=(cxx.*cyy-cxy.^2/4).*c0+cxy.*cx.*cy/4-cxx.*(cy.^2)/4-cyy.*(cx.^2)/4;

% local feature points by x=-H\g0, invH=1/detH*adjH
X=-(2*cyy.*cx-cxy.*cy)./detH;
Y=-(-cxy.*cx+2*cxx.*cy)./detH;

% compute feature point shift by d=-H\g0t
U=-(2*cyy.*cxt-cxy.*cyt)./detH;
V=-(-cxy.*cxt+2*cxx.*cyt)./detH;

% compute attenuation peak
Ip=c0+cx.*X+cy.*Y+cxx.*X.^2+cxy.*X.*Y+cyy.*Y.^2;
% It=ct-cx.*U-cy.*V;
% Itt=ctt+cxx.*U.^2+cxy.*U.*V+cyy.*V.^2;

% feature point subwindow
W=( sqrt((X/(hwin(1)+1)).^2+(Y/(hwin(2)+1)).^2) <= 1 ); % within window 1/2 [px] tolerance @ sub-pix

% displacement subwindow
D=( sqrt((U*hwin(3)/(2*hwin(1)+2)).^2+(V*hwin(3)/(2*hwin(2)+2)).^2) <= 1 );  % disp out of window + 1/2 [px] tolerance @ sub-pix

% mask for real minima maxima saddles and ellipses
switch typ
    case 'minima'
        
        % analytic mimima
        C=imsegscl(detH,0,std_param,iter_param) & imsegscl(trH,0,std_param,iter_param) & imsegscl(-Ip,0,std_param,iter_param); % 0,2,8
        
    case 'maxima'
        
        % analytic maxima
        C=imsegscl(detH,0,std_param,iter_param) & imsegscl(-trH,0,std_param,iter_param) & imsegscl(Ip,0,std_param,iter_param); % there is no definite setting in imsegscl that make all other obsolete
%         C=detH>0 & -trH>0 & Ip>0; % there is no definite setting in imsegscl that make all other obsolete
        
    case 'saddles'
        
        % analytic saddle
        C=imsegscl(-detH,0,std_param,iter_param); % actually sigma =1 sets quick convergence , allows to limit iter=4
        
    case 'ellipses' 
        
        % analytic real ellipse
        C=    ( imsegscl(detH,0,std_param,iter_param) & imsegscl(trH,0,std_param,iter_param) & imsegscl(-Ip,0,std_param,iter_param) ) ... 
            | ( imsegscl(detH,0,std_param,iter_param) & imsegscl(-trH,0,std_param,iter_param) & imsegscl(Ip,0,std_param,iter_param) ); % see note above @ C
        
end

% transform conics
c0=c0-cx.*x-cy.*y+cxx.*x.^2+cxy.*x.*y+cyy.*y.^2;
cx=cx-2*cxx.*x-cxy.*y;
cy=cy-cxy.*x-2*cyy.*y;

% feasible domain within window displacement and classification
F=W & D & C; % feasible set & msk

% linear indexing F 
indF=find(F);

% vote candidates - where is my most probable feature?
[Vmap,~,~,I,J]=histcounts2(x(F)+X(F),y(F)+Y(F),1/2+(0:size(imd,1)),1/2+(0:size(imd,2))); % keep discrete, do not smooth

%figure; surf(x',y',Vmap'); shading flat; view(2); axis tight; axis equal

% % average voting on size and image fit region for ellipse
% Vmap = imrelscl( Vmap , size(fobj.Neighborhood) ); % does not improve things
% [Vs,Vl] = imrelscl( Vmap , size(fobj.Neighborhood) );
% Vmap=Vs-Vl;

% remove zero indexing (outside bounds)
rem=I==0 | J==0; I(rem)=[]; J(rem)=[]; F(indF(rem))=0; indF(rem)=[];
indV=sub2ind(size(Vmap),I,J); % actually size(F) but note that equal to size(Vmap)

% maximize voting likelyhood - most probably midpoint
Vdom=(Vmap)>0 & F; % segment repeatibility vs remove noisy ident  -Bmap
Vmax=imdilate(Vmap,fobj)==Vmap & Vdom; % seed
V=Vdom;
while nnz(V-Vmax)>0
    V=Vmax;
    Vmax=bwmorph(V,'dilate',1) & Vdom; % 8 neigh groth inside i.e. clean edg?
end

% connected labels connected regions
Lmap=bwlabeln(V,8);
indL=find(V);
L=Lmap(indL);

% index mapping
dumV=sparse(indV,1:length(indV),ones(size(indV)),numel(Vmap),length(indV)); % automaticly weights occurence
dumL=sparse(indL,1:length(indL),ones(size(indL)),numel(Lmap),length(indL)); % numel(Lmap)==numel(Vmap)
FL=double(dumV'*dumL>0.5); % bsxfun(@eq,sparse(indV),sparse(indL)'); % optimize..d
%figure;spy(FL)
F=zeros(size(F));
F(indF)=FL*L; % map to labels

% retrieve coefficient after averaging cell2mat(struct2cell(regionprops(F,Vlev,'MeanIntensity'))')./
W=cell2mat(struct2cell(regionprops(F,'Area'))');
c0=cell2mat(struct2cell(regionprops(F,c0,'MeanIntensity'))');
cx=cell2mat(struct2cell(regionprops(F,cx,'MeanIntensity'))');
cy=cell2mat(struct2cell(regionprops(F,cy,'MeanIntensity'))');
cxx=cell2mat(struct2cell(regionprops(F,cxx,'MeanIntensity'))');
cxy=cell2mat(struct2cell(regionprops(F,cxy,'MeanIntensity'))');
cyy=cell2mat(struct2cell(regionprops(F,cyy,'MeanIntensity'))');
cxt=cell2mat(struct2cell(regionprops(F,cxt,'MeanIntensity'))');
cyt=cell2mat(struct2cell(regionprops(F,cyt,'MeanIntensity'))');
ct=cell2mat(struct2cell(regionprops(F,ct,'MeanIntensity'))');
ctt=cell2mat(struct2cell(regionprops(F,ctt,'MeanIntensity'))'); % regionprops is slowest part here

% compute determinant and trace hessian
trH=2*cxx+2*cyy; % trace
detH=4*cxx.*cyy-cxy.^2; % determinant

% local feature points by x=-H\g0, invH=1/detH*adjH
X=-(2*cyy.*cx-cxy.*cy)./detH;
Y=-(-cxy.*cx+2*cxx.*cy)./detH;

% compute feature point shift by d=-H\g0t
U=-(2*cyy.*cxt-cxy.*cyt)./detH;
V=-(-cxy.*cxt+2*cxx.*cyt)./detH;

% compute attenuation peak -  here @ fluctuation of the int levels (not variations in c0)
Ip=c0+cx.*X+cy.*Y+cxx.*X.^2+cxy.*X.*Y+cyy.*Y.^2;
It=ct;%+cx.*U+cy.*V;
Itt=ctt;%-cxx.*U.^2-cxy.*U.*V-cyy.*V.^2;

% time varying peak value over frames
T=( (1:size(imd,3)) - foi ).^( (0 : 2)' );
if ~isempty(Ip)
    segI=[Ip It Itt]*T;
else
    segI=Ip;
end

% segment final data for real minima maxima saddles and ellipses
switch typ
    case 'minima'
        
        % analytic mimima
        wrt=detH>0 & trH>0 & max(segI,[],2)<0; %Ip
        
    case 'maxima'
        
        % analytic maxima
        wrt=detH>0 & trH<0 & min(segI,[],2)>0;
        
    case 'saddles'
        
        % analytic saddles
        wrt=detH<0;
        
    case 'ellipses' 
        
        % analytic ellipses
        wrt=  ( detH>0 & trH>0 & max(segI,[],2)<0 ) ... 
            | ( detH>0 & trH<0 & min(segI,[],2)>0 );
        
end

% write conics
len=length(c0(wrt));
CC=zeros(len,6);
CC(:,1:6)=[cxx(wrt),cxy(wrt),cyy(wrt),cx(wrt),cy(wrt),c0(wrt)];

% write feature points
XX=zeros(len,2);
XX(:,1:2)=[X(wrt),Y(wrt)];

% write displacement
UU=zeros(len,2);
UU(:,1:2)=[U(wrt),V(wrt)];

% write intensity information
II=zeros(len,3);
II(:,1:3)=[Ip(wrt) It(wrt) Itt(wrt)];

% weight from repeated similar identifycation
WW=zeros(len,1);
WW(:,1)=W(wrt);

end