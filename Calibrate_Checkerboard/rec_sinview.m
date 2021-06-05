function rec_sinview
%rec_sinview Rectify single view to obtain the distortion mapping and 
%   distortion correct the camera-plane by minimizing the point-line 
%   distance. 
%   
%   Input:
%       Processed image data
%       filesel.dis to indicate selected distortion files
%   
%   Output:
%       DWblin.dat Dewarped baseline
%       pld.dat Point-(on curve)-line distance from minimization
%       Dpar.mat Distortion parameters in .mat file

%% get globals
global folder date cal rec ctrl prop Cfit Dmap filesel

%% load data
IMcfit=cell(2,length(filesel));
IMblin=cell(2,length(filesel));
for k=1:length(filesel) % new data set indicator
    if ~isempty(filesel(k).dis)
        n=unique(filesel(k).dis(:,1));
        c=unique(filesel(k).dis(:,2));
        
        imcfit=fload([folder date cal vsl filesel(k).name '_IMcfit.dat']);
        imcfit=imcfit(:,ismember(imcfit(1,:),n) & ismember(imcfit(2,:),c));
        
        IMcfit{2,k}=imcfit;
        IMcfit{1,k}=k*ones(1,size(IMcfit{2,k},2));
        
        imblin=fload([folder date cal vsl filesel(k).name '_IMblin.dat']);
        imblin=imblin(:,ismember(imblin(1,:),n) & ismember(imblin(2,:),c));
        
        IMblin{2,k}=imblin;
        IMblin{1,k}=k*ones(1,size(IMblin{2,k},2));
    end
end
IMcfit=cell2mat(IMcfit);
IMblin=cell2mat(IMblin);

%% initiate nummeric integration
dt=1/10; % numeric integration stepsize
t=linspace(0,1,1/dt); % and steps; later we use sum*dt=mean in 0 to 1

%% Define objective function
% substitute path in mapping
x=sym(Dmap);
b=sort(symvar(x));
b=b(1:end-2);

% subs
x=subs(x,sym('x',[2 1]),sym(Cfit));

% line
l=sym('c',[2 1]); % inh norm def.
delx=subs(x,sym('t'),t(end))-subs(x,sym('t'),t(1)); % simple def

% def global opt problem
ObjF=(x'*l+1)^2; % point line distance (!) 
wObjF=1/dot(delx,delx)/dot(l,l); % non dimensionalize by length on curve

% gradient and hessian iObjF
vars=sort(symvar(ObjF));

% matlab function file
objf=matlabFunction(ObjF,'file',[folder date cal vsl 'mfunc' vsl 'objf.m' ],'Vars',vars);%,'Optimize',false);
wobjf=matlabFunction(wObjF,'file',[folder date cal vsl 'mfunc' vsl 'wobjf.m' ],'Vars',vars);%,'Optimize',false);

%% Initiate parameters

% no distortion map
switch ctrl.dmap
    case 'interface-mod'
        b0=[0
            -1/2*prop.res(1)
            1
            -1/2*prop.res(2)
            0
            0
            (1-prop.nref^2)/((prop.flens/sqrt(prop.areachip/prod(prop.res(1:2))))^2)]; %zeros(length(nu),1);%[0 0 0 0 0]';%
        H0=reshape([1;b0(1:2);0;b0(3:6);1],3,3)';
    case 'interface-cor'
        b0=zeros(length(b),1);
        b0(1:7)=[0
            -1/2*prop.res(1)
            1
            -1/2*prop.res(2)
            0
            0
            (1-prop.nref^2)/((prop.flens/sqrt(prop.areachip/prod(prop.res(1:2))))^2)];% initial distortion can be scaled / ignored
        H0=reshape([1;b0(1:2);0;1;b0(4:6);1],3,3)';
    otherwise
        b0=zeros(length(b),1);
        H0=eye(3);
end

%% Find initial best estimate line to curve from baseline
disp('Fit base lines')

% numbering knowns
nk=[1:length(vars)-3,length(vars)]';

% numbering unknowns
nu=(1:length(vars))';
nu(nk)=[];

% take gradient
[gObjF,hObjF]=difobjf(ObjF,vars(nu)); 
gobjf=matlabFunction(gObjF,'file',[folder date cal vsl 'mfunc' vsl 'gobjf.m' ],'Vars',vars);%,'Optimize',false);
hobjf=matlabFunction(hObjF,'file',[folder date cal vsl 'mfunc' vsl 'hobjf.m' ],'Optimize',false,'Vars',vars);

% solution
pld=zeros(8,size(IMcfit,2));
pld(1:4,:)=IMcfit(1:4,:); % pointline distance stats
DWblin=zeros(7,size(IMcfit,2));
DWblin(1:4,:)=IMcfit(1:4,:); % dewarped baseline
for c=1:prop.res(4)
    % get curve data
    C=IMcfit(3,:)==c;
    a=IMcfit(5:end,C); % coef data
    
    % get baseline data
    l=homc2inhc(inv(H0)'*(IMblin(5:end,C)./IMblin(end,C))); %
    
    % initial sol vec.
    x0=reshape(l,[],1);
    
    % numering assembly
    asn=reshape(repmat(reshape((1:numel(l)),2,[]),length(t),1),2,[]); % assembly numbering
    
    % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
    dat=[zeros(size(nk)) nk [reshape(repmat(a,length(t),1),size(a,1),[]);repmat(b0,1,length(t)*size(a,2));repmat(t,1,size(a,2))]
        ones(size(nu)) nu asn];
    
    % solve global optimization at low precision improving intial estimate baselines
    [xsol,~] = optobjf( wobjf, objf, gobjf, hobjf , x0 , dat , ctrl.optl ,'newton');% , 'disp', 1, 2);
    
    % write lines
    l=reshape(xsol,2,[]);
    
    % initiated point line distance
    inp=num2cell([reshape(repmat(a,length(t),1),size(a,1),[]);...
        repmat(b0,1,size(a,2)*length(t));...
        reshape(repmat(l,length(t),1),size(l,1),[]);...
        repmat(t,1,size(a,2))],2);
    pld(5,C)=mean(reshape(sqrt(objf(inp{:})),length(t),[]),1); % distance
    pld(6,C)=mean(reshape(sqrt(objf(inp{:}).*wobjf(inp{:})),length(t),[]),1); % normalized distance
    
    % write line data
    DWblin(5:end,C)=inhc2homc(reshape(xsol(1:end),2,[])); % coef data
    
end % c

%% solve distortion mappings and plot results

disp('Distortion correct views')

% numbering knowns
nk=[1:2*(ctrl.cord+1),length(vars)]';

% numbering unknowns
nu=(1:length(vars))';
nu(nk)=[];

tic
[gObjF,hObjF]=difobjf(ObjF,vars(nu));
gobjf=matlabFunction(gObjF,'file',[folder date cal vsl 'mfunc' vsl 'gobjf.m' ],'Vars',vars);%,'Optimize',false);
hobjf=matlabFunction(hObjF,'file',[folder date cal vsl 'mfunc' vsl 'hobjf.m' ],'Optimize',false,'Vars',vars);
toc

% solution
Dpar=cell(1,prop.res(4));
for c=1:prop.res(4)
    % get curve data
    C=IMcfit(3,:)==c;
    a=IMcfit(5:end,C); % coef data
    
    % get baseline data
    l=homc2inhc(DWblin(5:end,C));
    
    % initial sol vec.
    x0=[b0;reshape(l,[],1)];
    
    % numering assembly
    asn=[repmat((1:length(b0))',1,length(t)*size(a,2))
        reshape(repmat(reshape(...
        length(b0)+(1:numel(l)),2,[]),length(t),1),2,[])]; % assembly numbering
    
    % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
    dat=[zeros(size(nk)) nk [reshape(repmat(a,length(t),1),size(a,1),[]);repmat(t,1,size(a,2))]
        ones(size(nu)) nu asn];
    
    % solve global optimization
    [xsol,~] = optobjf( wobjf, objf, gobjf, hobjf , x0 , dat , ctrl.optl ,'newton');
    
    % write solution vector ( trash/forget redefined lines )
    bsol=xsol(1:length(b0));
    
    % get curve data
    DWblin(5:end,C)=inhc2homc(reshape(xsol(length(b0)+1:end),2,[])); % coef data
    
    % initial point line distance
    l=homc2inhc(DWblin(5:end,C));
    inp=num2cell([reshape(repmat(a,length(t),1),size(a,1),[]);...
        repmat(bsol,1,size(a,2)*length(t));...
        reshape(repmat(l,length(t),1),size(l,1),[]);...
        repmat(t,1,size(a,2))],2);
    pld(7,C)=mean(reshape(sqrt(objf(inp{:})),length(t),[]),1); % distance
    pld(8,C)=mean(reshape(sqrt(objf(inp{:}).*wobjf(inp{:})),length(t),[]),1); % normalized distance
    
    % some image
    imd=double(import_frames({folder date cal rec},prop.ext,2,c)); % board scatters
    
    % downgrade resolution
    [y,x]=meshgrid(1:3:prop.res(2),1:3:prop.res(1));
    imd=interp2(imd,y,x);
    
    % plot
    figure
    subplot(1,3,1)
    surf(x',y',log(imd'+1))
    view(2)
    shading flat
    axis tight
    axis equal
    colormap gray
    xlabel('x [px]','interpreter','latex')
    ylabel('y [px]','interpreter','latex')
    title(['Camera = ',num2str(c)],'interpreter','latex')
    drawnow
    
	% compute homography to reference plane
	xref=[1 prop.res(1) prop.res(1) 1
	      1 1           prop.res(2) prop.res(2)]; % line near horizon image
    inp=num2cell([repmat(bsol,1,size(xref,2)) ; xref],2);
    xcor=Dmap(inp{:}); % deformed
    H=dirlintrans(inhc2homc(xcor),inhc2homc(xref)); % rotation part
    DWblin(5:end,C)=H'\DWblin(5:end,C);
    
    % write solution vector
    Dpar{c}.map=bsol';
	Dpar{c}.H=H;
	
    % dewarp distortionsym('b',[1 length(Dpar{c})])
    inp=num2cell([repmat(Dpar{c}.map',1,size(x(:)',2)) ; x(:)';y(:)'],2);
    XX=Dmap(inp{:}); % deformed
    Y=reshape(XX(2,:),size(x));
    X=reshape(XX(1,:),size(x));
    
    % plot distortion mapping
    subplot(1,3,2)
    surf(X',Y',log(imd'+1))
    view(2)
    shading flat
    axis tight
    axis equal
    colormap gray
    xlabel('x [px]','interpreter','latex')
    ylabel('y [px]','interpreter','latex')
    title('Dewarped Image','interpreter','latex')
    drawnow
    
    % dewarp distortionsym('b',[1 length(Dpar{c})])
    XX=homc2inhc(Dpar{c}.H*inhc2homc(XX)); % deformed
    Y=reshape(XX(2,:),size(x));
    X=reshape(XX(1,:),size(x));
    
    % plot reference plane by image resolution
    subplot(1,3,3)
    surf(X',Y',log(imd'+1))
    view(2)
    shading flat
    axis tight
    axis equal
    colormap gray
    xlabel('x [px]','interpreter','latex')
    ylabel('y [px]','interpreter','latex')
    title('Dewarped Image Regularized Reference Plane','interpreter','latex')
    drawnow
    
end % c

%% residual data ksq and distortion
for f=unique(pld(1,:))
    % file indicator
    F1=DWblin(1,:)==f;
    F2=pld(1,:)==f;
    
    % write optimized lines
    fsave([folder date cal vsl filesel(f).name '_DWblin.dat'],DWblin(2:end,F1),'w'); % overwrite
        
    % write residual
    fsave([folder date cal vsl filesel(f).name '_pld.dat'],pld(2:end,F2),'w'); % overwrite
    
end

% Save distortion parameters
save([folder date cal vsl 'Dpar.mat'],'Dpar')% save mapping

end