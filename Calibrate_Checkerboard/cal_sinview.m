function cal_sinview
%cal_sinview Calibrate single view by optimizing the pixel errors for each
%   view.
%   
%   Input:
%       Processed image data
%       Solution to the distortion map
%       filesel.cal to indicate selected calibration files
%   
%   Output:
%       rpe.dat reprojection errors
%       hcb.dat Initial image homography for each checkerboard
%       tcb.dat optimized translation vectors for each checkerboard
%       qcb.dat optimized unit quaternion for each checkerboard
%       Kmat.mat Camera calibration matrix for each view

%% Get global variables
global folder date cal ctrl prop chkboard Dmap Pmap rec filesel

%% load data
IMcnod=cell(2,length(filesel));
for k=1:length(filesel) % new data set indicator
    if ~isempty(filesel(k).cal)
        n=unique(filesel(k).cal(:,1));
        c=unique(filesel(k).cal(:,2));
        
        imcnod=fload([folder date cal vsl filesel(k).name '_IMcnod.dat']);
        imcnod=imcnod(:,ismember(imcnod(1,:),n) & ismember(imcnod(2,:),c));
        
        IMcnod{2,k}=imcnod;
        IMcnod{1,k}=k*ones(1,size(IMcnod{2,k},2));
    end
end
IMcnod=cell2mat(IMcnod);

%% load previous
Dpar=importdata([folder date cal vsl 'Dpar.mat']);

%% create objective function fit homography
% image and checkerboard defined coordinates
x=sym('x',[2 1]);
X=sym('X',[2 1]);

% define argument objective function
ObjA=x-subs(sym(Pmap),sym('X',[4 1]),[X;0;1]);

% parametrize homography up to fixed scale
H=reshape(sym('h',[9,1]),3,3)';
H(end,end)=1;

% substitute homography
P=sym('P',[3 4]);
P(:,3)=[];
ObjA=subs(ObjA,P,H);

% define objective function
ObjF=defobjf(ObjA,[]);

% variabes
vars=sort(symvar(ObjF));

% write function
objf=matlabFunction(ObjF,'file',[folder date cal vsl 'mfunc' vsl 'objf.m'],'Vars',vars);%,'Optimize',false);

%% get homography between camera's and checkerboard
% numbering knowns
nk=[1,2,11,12]';

% numbering unknowns
nu=(1:length(vars))';
nu(nk)=[];

% find grad and hess w.r.t. unknowns
[gObjF,hObjF]=difobjf(ObjF,vars(nu));
gobjf=matlabFunction(gObjF,'file',[folder date cal vsl 'mfunc' vsl 'gobjf.m'],'Vars',vars);%,'Optimize',false);
hobjf=matlabFunction(hObjF,'file',[folder date cal vsl 'mfunc' vsl 'hobjf.m'],'Vars',vars);

% initiate Hcb
Hcb=zeros(12,size(IMcnod,2)/prod(ctrl.nnod)); % initiate data cell Blines
rpe=zeros(6,size(IMcnod,2));
rpe(1:4,:)=IMcnod(1:4,:);
for f=unique(IMcnod(1,:))
    % file indicator set
    F=IMcnod(1,:)==f;
    
    for n=unique(IMcnod(2,F))
        % frame indicator
        N=IMcnod(2,:)==n;
        
        for c=unique(IMcnod(3,F&N))
            % camera ind.
            C=IMcnod(3,:)==c;
            
            % camera data
            x=IMcnod(5:6,F&N&C);
            
            % dewarped coordinates
            inp=num2cell([repmat(Dpar{c}.map',1,prod(ctrl.nnod))
                x],2);
            x=Dmap(inp{:}); % distortion mapping
            x=homc2inhc(Dpar{c}.H*inhc2homc(x)); % back to image resolution
            
            %figure; plot(x(1,:),x(2,:),'.')
            
            % chk nodes
            X=chkboard.points; % [X Y - 1] 
            
            % hom coord.
            x=inhc2homc(x);
            X=inhc2homc(X);
            
            % homography checkerboard
            h=reshape(dirlintrans(X,x)',[],1);
            h(end)=[];
            
            % back to homogenous coordinates
            X=chkboard.points;
            x=homc2inhc(x);
            
            % initial reproj error
            inp=num2cell([X;repmat(h,1,size(X,2));x],2);
            rpe(5,F&N&C)=sqrt(objf(inp{:}));
            
            % numering assembly
            asn=repmat((1:length(nu))',1,size(x,2)); % assembly numbering
            
            % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
            dat=[zeros(size(nk)) nk [X;x]
                ones(size(nu)) nu asn];
            
            % initiate solution vector
            h0=h; % or do 1 ord coef a
            
            % solve global optimization (newton is here exact solution, but do generalize)
            hsol = optobjf( [], objf, gobjf, hobjf , h0 , dat ,  ctrl.optl ,'newton','no');
            
            % optimal reproj error
            inp=num2cell([X;repmat(hsol,1,size(X,2));x],2);
            rpe(6,F&N&C)=sqrt(objf(inp{:}));
            
            % write
            ind=(find(F&N&C,1,'first')-1)/prod(ctrl.nnod)+1;
            Hcb(:,ind)=[f
                n
                c
                hsol
                1]; % same ordering unstruct data
            
        end
    end
    
    % message
    disp(['fit chkb. fil. ',num2str(f)])
    
    % save homography
    F=Hcb(1,:)==f;
    fsave([folder date cal vsl filesel(f).name '_Hcb.dat'],Hcb(2:end,F),'w'); % overwrite

end

%% create objective function fit parametrized homography
% image and checkerboard defined coordinates
x=sym('x',[2 1]);
X=sym('X',[2 1]);

% define argument objective function
ObjA=x-subs(sym(Pmap),sym('X',[4 1]),[X;0;1]);

% parametrize homography up to fixed scale
k=sym('k',[5,1]);
K=[k(1) k(2) k(3)
    0 k(4) k(5)
    0 0 1];

q=sym('q',[4,1]);
t=sym('t',[3 1]);
R=qcom2rotm(q);
H=K*[R(:,1:2),t]; %

% substitute homography
P=sym('P',[3 4]);
P(:,3)=[];
ObjA=subs(ObjA,P,H);

% constraints
eqcF=dot(q,q)-1;

% define objective function
ObjF=defobjf(ObjA,eqcF);

% variables
vars=sort(symvar(sym(ObjF)));

% write function
objf=matlabFunction(ObjF,'file',[folder date cal vsl 'mfunc' vsl 'objf.m'],'Vars',vars);%,'Optimize',false);

%% Compute camera matrices by homographies
Kmat=cell(1,prop.res(4));
qcb=zeros(7,size(Hcb,2));
tcb=zeros(6,size(Hcb,2));
rpe(7:8,:)=zeros(2,size(rpe,2)); % 7:10
for c=unique(Hcb(3,:))
    % camera set indicators
    C1=Hcb(3,:)==c;
    C2=IMcnod(3,:)==c;
    
    % get data
    Hdat=cellfun(@(x)reshape(x,3,3)',num2cell(Hcb(4:end,C1),1),'UniformOutput',false);% the @isempty vs 'isempty' speed bug has been resolved it looks like
    
    % asolute conic by homographies reference object
    [ W ] = absconhom( Hdat );
    
    % camera matrix
    K=inv(chol(W));
    K=K/K(end); % normalize scale
    
    % get camera parameters
    k=kmat2cpar(K)';
    
    % numbering knowns
    nk=[1:7,15,16]';
    
    % numbering unknowns
    nu=(1:length(vars))';
    nu(nk)=[];
            
    % find grad and hess w.r.t. unknowns
    [gObjF,hObjF]=difobjf(ObjF,vars(nu));
    gobjf=matlabFunction(gObjF,'file',[folder date cal vsl 'mfunc' vsl 'gobjf.m'],'Optimize',false,'Vars',vars);%,'Optimize',false);
    hobjf=matlabFunction(hObjF,'file',[folder date cal vsl 'mfunc' vsl 'hobjf.m'],'Optimize',false,'Vars',vars);
    
    % find rigid body motion checkerboard, per file per board
    for f=unique(Hcb(1,C1))
        % file set indicators
        F1=Hcb(1,:)==f;
        F2=IMcnod(1,:)==f;
        
        for n=unique(Hcb(2,F1&C1))
            % frame set inidicators
            N1=Hcb(2,:)==n;
            N2=IMcnod(2,:)==n;
            
            % checkerboard homography
            H=reshape(Hcb(4:end,C1&F1&N1),3,3)';
            
            % homography best R-t decomposition
            Q=K\H;
            lam=0.5/norm(Q(:,1),2)+0.5/norm(Q(:,2),2);
            r1=lam*Q(:,1);
            r2=lam*Q(:,2);
            r3=cross(r1,r2);
            R=[r1 r2 r3];
            [R,~]=frotm(R); % incl stetch
            t=lam*Q(:,3);%S\
            
            % quaternion
            q=rotm2qcom(R);
            
            %optimize reprojection
            X=chkboard.points;
            x=(IMcnod(5:6,F2&N2&C2));
            
            % dewarped coordinates
            inp=num2cell([repmat(Dpar{c}.map',1,prod(ctrl.nnod))
                x],2);
            x=Dmap(inp{:}); % distortion mapping
            x=homc2inhc(Dpar{c}.H*inhc2homc(x)); % back to image resolution
            
            % initial reproj error
            inp=num2cell([X;repmat([k;q;t],1,size(X,2));x;zeros(1,size(X,2))],2);
            rpe(7,F2&N2&C2)=objf(inp{:});
            
            % numering assembly w.r.t solution vector
            asn=repmat((1:length(nu))',1,size(x,2)); % assembly numbering
            
            % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
            dat=[zeros(size(nk)) nk [X;repmat(k,1,size(x,2));x]
                ones(size(nu)) nu asn];
            
            % initiate solution vector
            x0=[q
                t
                0]; % or do 1 ord coef a
            
            % solve global optimization 
            xsol = optobjf( [], objf, gobjf, hobjf , x0 , dat , ctrl.optl,'newton','no');
            
            % write least data repres
            q=xsol(1:4)/norm(xsol(1:4),2);
            t=xsol(5:7);
            
            % final reproj error
%             inp=num2cell([X;repmat([k;q;t],1,size(X,2));x;zeros(1,size(X,2))],2); % not need mult, due corr quat
%             rsq(8,F2&N2&C2)=objf(inp{:});
            
            % write            
            qcb(:,C1&F1&N1)=[f
                n
                c
                q];
            tcb(:,C1&F1&N1)=[f
                n
                c
                t]; % same order as Hcb
        end
        
        % message
        disp(['pos. chkb. cam. ',num2str(c),' fil. ',num2str(f)])
        
    end
    
    % message
    disp('ref. cam. mat.')
        
    % numbering knowns
    nk=[1:2,8:17]';
    
    % numbering unknowns
    nu=(1:length(vars))';
    nu(nk)=[];
    
    % find grad and hess. w.r.t opt. rot m
    [gObjF,hObjF]=difobjf(ObjF,vars(nu));
    gobjf=matlabFunction(gObjF,'file',[folder date cal vsl 'mfunc' vsl 'gobjf.m' ],'Optimize',false,'Vars',vars);%,'Optimize',false);
    hobjf=matlabFunction(hObjF,'file',[folder date cal vsl 'mfunc' vsl 'hobjf.m' ],'Optimize',false,'Vars',vars);
        
    % optimize reprojection
    x=(IMcnod(5:6,C2));
    
    % dewarped coordinates
    inp=num2cell([repmat(Dpar{c}.map',1,size(x,2))
        x],2);
    x=Dmap(inp{:}); % distortion mapping
    x=homc2inhc(Dpar{c}.H*inhc2homc(x)); % back to image resolution
    
    % param board.
    X=repmat(chkboard.points,1,size(x,2)/prod(ctrl.nnod));
    
    % quaternion and transl
    q=cell2mat(cellfun(@(x)repmat(x,1,prod(ctrl.nnod)),num2cell(qcb(4:end,C1),1),'UniformOutput',false));
    t=cell2mat(cellfun(@(x)repmat(x,1,prod(ctrl.nnod)),num2cell(tcb(4:end,C1),1),'UniformOutput',false));
    
    % numering assembly w.r.t solution vector
    asn=repmat((1:length(nu))',1,size(x,2)); % assembly numbering
    
    % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
    dat=[zeros(size(nk)) nk [X;q;t;x;zeros(1,size(q,2))]
        ones(size(nu)) nu asn];
    
    % initiate solution vector
    k0=k; % or do 1 ord coef a
    
    % solve global optimization
    ksol = optobjf( [], objf, gobjf, hobjf , k0 , dat , ctrl.optl ,'newton');
    
    % simulateneous refinement
    k=ksol;
    
    % final reproj error
%     inp=num2cell([X;repmat(k,1,size(X,2));q;t;x;zeros(1,size(X,2))],2); % not need mult, due corr quat
%     rsq(9,C2)=objf(inp{:});
    
    % message
    disp('ref. cam. cal.')
     
    % numbering knowns
    nk=[1,2,15,16]';
    
    % numbering unknowns
    nu=(1:length(vars))';
    nu(nk)=[];
    
    % find grad and hess. w.r.t opt. rot m
    [gObjF,hObjF]=difobjf(ObjF,vars(nu));%can upgrade speed here, by only diff to need
    gobjf=matlabFunction(gObjF,'file',[folder date cal vsl 'mfunc' vsl 'gobjf.m' ],'Optimize',false,'Vars',vars);%,'Optimize',false);
    hobjf=matlabFunction(hObjF,'file',[folder date cal vsl 'mfunc' vsl 'hobjf.m' ],'Optimize',false,'Vars',vars);
    
    % quaternion and camera parameters
    q=qcb(4:end,C1);
    t=tcb(4:end,C1);
    
    %optimize reprojection
    x=(IMcnod(5:6,C2));
    
    % dewarped coordinates
    inp=num2cell([repmat(Dpar{c}.map',1,size(x,2))
        x],2);
    x=Dmap(inp{:}); % distortion mapping
    x=homc2inhc(Dpar{c}.H*inhc2homc(x)); % back to image resolution
    
    % param board.
    X=repmat(chkboard.points,1,size(x,2)/prod(ctrl.nnod));
        
    % initiate solution vector
    x0=[k
        reshape([q;t;zeros(1,size(q,2))],[],1)]; % or do 1 ord coef a
    
    % numering assembly w.r.t solution vector
    asn=[repmat((1:5)',1,size(x,2))
        reshape(...
        repmat(...
        reshape(6:length(x0),[],size(x,2)/prod(ctrl.nnod))...
        ,prod(ctrl.nnod),1)...
        ,[],size(x,2))]; % assembly numbering
    
    % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
    dat=[zeros(size(nk)) nk [X;x]
        ones(size(nu)) nu asn];
    
    % solve global optimizationclc
    xsol = optobjf( [], objf, gobjf, hobjf , x0 , dat , ctrl.optl ,'newton');
    
    % write
    k=xsol(1:5);
    xsol=reshape(xsol(6:end),8,[]);
    q=xsol(1:4,:)./sqrt(sum(xsol(1:4,:).^2,1));
    t=xsol(5:7,:);
    
    % write
    Kmat{c}=cpar2kmat(k(1),k(2),k(3),k(4),k(5));
    qcb(4:end,C1)=q;
    tcb(4:end,C1)=t;
    
    % final reproj error% quaternion and transl
    q=cell2mat(cellfun(@(x)repmat(x,1,prod(ctrl.nnod)),num2cell(qcb(4:end,C1),1),'UniformOutput',false));
    t=cell2mat(cellfun(@(x)repmat(x,1,prod(ctrl.nnod)),num2cell(tcb(4:end,C1),1),'UniformOutput',false));
    inp=num2cell([X;repmat(k,1,size(X,2));q;t;x;zeros(1,size(X,2))],2); % not need mult, due corr quat
    rpe(8,C2)=sqrt(objf(inp{:})); %10
    
end

%% save
for f=unique(qcb(1,:))
    % file indicator
    F1=qcb(1,:)==f;
    F2=rpe(1,:)==f;
    
    % rigid body motion
    fsave([folder date cal vsl filesel(f).name '_qcb.dat'],qcb(2:end,F1),'w'); % overwrite
    fsave([folder date cal vsl filesel(f).name '_tcb.dat'],tcb(2:end,F1),'w'); % overwrite
    fsave([folder date cal vsl filesel(f).name '_rpe.dat'],rpe(2:end,F2),'w'); % overwrite
    
end
save([folder date cal vsl 'Kmat.mat'],'Kmat')
  
%% Plot image-plane
% show dewarped and calibrated image K\Distcor..
for c=1:prop.res(4)
    % some figure
    figure
    
    % show result of calibrated image
    imd=log(double(import_frames({folder date cal rec},prop.ext,1,c))+1); % board scatters
    
    % downgrade resolution
    [y,x]=meshgrid(1:3:prop.res(2),1:3:prop.res(1));
    imd=interp2(imd,y,x); % downsample
    
    %dewarp distortion
    map=subfunv(Dmap,sym('b',[1 length(Dpar{c}.map)]),Dpar{c}.map);
    X=map(x(:)',y(:)'); % distortion
    X=homc2inhc(Dpar{c}.H*inhc2homc(X)); % image reference
    
    % rectify
    X=homc2inhc(Kmat{c}\inhc2homc(X));
    
    % meshgrid
    Y=reshape(X(2,:),size(x));
    X=reshape(X(1,:),size(x));
    
    % plot
    surf(X',Y',imd')
    view(2)
    shading flat
    axis tight
    axis equal
    colormap gray
    xlabel('x $[\rm{m}]$','interpreter','latex')
    ylabel('y $[\rm{m}]$','interpreter','latex')
    title(['Calibrated camera ',num2str(c)],'interpreter','latex')
    drawnow
    
end

%% Plot calibration results
for c=1:prop.res(4)
    % some figure
    figure
    
    % camera set indicators
    C=qcb(3,:)==c;
    
    % find rigid body motion checkerboard, per file per board
    for f=unique(qcb(1,C))
        % file set indicators
        F=qcb(1,:)==f;
        
        for n=unique(qcb(2,F&C))
            % frame set inidicators
            N=qcb(2,:)==n;
            
            % get points
            X=inhc2homc(chkboard.points,0);
            
            % get transformation
            R=qcom2rotm(qcb(4:end,F&N&C));
            t=tcb(4:end,F&N&C);
            
            % position board
            X=R*X+t;
            
            %plotting
            hold on
            plot3(X(1,:),X(2,:),X(3,:),'.')
            hold off
            
        end
    end
    
    % axis 
    axis tight
    axis equal
    view(3)
    grid minor
    xlabel('x $[\rm{m}]$','interpreter','latex')
    ylabel('y $[\rm{m}]$','interpreter','latex')
    zlabel('z $[\rm{m}]$','interpreter','latex')
    box on
    title(['Checkerboard positioning view ',num2str(c)],'interpreter','latex')
    drawnow
    
end

end