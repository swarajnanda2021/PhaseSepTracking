function cal_mulview
%cal_mulview Calibrate multiple view consistently by minimizing the spatial
%   error for all views together.
%   
%   Input:
%       Processed image data
%       Solution to the distortion map
%       Single view calibration
%       filesel.cal to indicate selected calibration files
%   
%   Output:
%       rpe.dat reprojection errors
%       tcb.dat refined translation vectors for each checkerboard
%       qcb.dat refined unit quaternion for each checkerboard
%       Kmat.mat Camera calibration matrix for each view
%       tvec.mat Camera pose translation vector
%       Rmat.mat Camera pose rotation matrix

%% Get global variables
global folder date cal ctrl chkboard Pmap Dmap filesel prop % prop rec

%% load data
IMcnod=cell(2,length(filesel));
qcb=cell(2,length(filesel));
tcb=cell(2,length(filesel));
rpe=cell(2,length(filesel));
for k=1:length(filesel) % new data set indicator
    if ~isempty(filesel(k).cal)
        n=unique(filesel(k).cal(:,1));
        c=unique(filesel(k).cal(:,2));
        
        imcnod=fload([folder date cal vsl filesel(k).name '_IMcnod.dat']);
        imcnod=imcnod(:,ismember(imcnod(1,:),n) & ismember(imcnod(2,:),c));
        
        IMcnod{1,k}=k*ones(1,size(imcnod,2));
        IMcnod{2,k}=imcnod;
        
        q=fload([folder date cal vsl filesel(k).name '_qcb.dat']);
        q=q(:,ismember(q(1,:),n) & ismember(q(2,:),c));
        
        qcb{1,k}=k*ones(1,size(q,2));
        qcb{2,k}=q;
        
        t=fload([folder date cal vsl filesel(k).name '_tcb.dat']);
        t=t(:,ismember(t(1,:),n) & ismember(t(2,:),c));
        
        tcb{1,k}=k*ones(1,size(t,2));
        tcb{2,k}=t;
        
        r=fload([folder date cal vsl filesel(k).name '_rpe.dat']);
        r=r(:,ismember(r(1,:),n) & ismember(r(2,:),c));
        
        rpe{1,k}=k*ones(1,size(r,2));
        rpe{2,k}=r;
    end
end
IMcnod=cell2mat(IMcnod);
qcb=cell2mat(qcb);
tcb=cell2mat(tcb);
rpe=cell2mat(rpe);

%% load previous
Dpar=importdata([folder date cal vsl 'Dpar.mat']);
Kmat=importdata([folder date cal vsl 'Kmat.mat']);

%% define lsqlin solver options
opt=optimoptions('lsqlin','Algorithm','interior-point'); % choose by doc (doc is two-fold)

%% find the relative rotation and translation between the views
c=unique(qcb(3,:));
X=inhc2homc(chkboard.points,0);
Rvw=cell(max(qcb(3,:)),max(qcb(3,:))); %relative
tvw=cell(max(qcb(3,:)),max(qcb(3,:))); %relative
for m=c
    for s=c
        
        % cam index
        M=qcb(3,:)==m;
        S=qcb(3,:)==s;
        
        % rigid body motion
        qm=qcb(4:end,M)';
        qs=qcb(4:end,S)';
        tm=tcb(4:end,M)';
        ts=tcb(4:end,S)';
        
        % overlap
        fm=qcb(1,M);
        fs=qcb(1,S);
        nm=tcb(2,M);
        ns=tcb(2,S);
        O=bsxfun(@eq,fm',fs) & bsxfun(@eq,nm',ns);
        
        % map data
        qs=O*qs;
        ts=O*ts;
        
        % squeeze data
        o=sum(O,2)==1;
        qm=qm(o,:)';
        qs=qs(o,:)';
        tm=tm(o,:)';
        ts=ts(o,:)';
        
        % rigid body motion
        Rm=cell2mat(cellfun(@(x)qcom2rotm(x),num2cell(qm,1),'UniformOutput',false)');
        tm=reshape(tm,[],1);
        Rs=cell2mat(cellfun(@(x)qcom2rotm(x),num2cell(qs,1),'UniformOutput',false)');
        ts=reshape(ts,[],1);
                
        % triangulation (not separate rot trans, this would be weighted rot trans)
        Xm=reshape(Rm*X+tm,3,[]);
        Xs=reshape(Rs*X+ts,3,[]);
        
        % use kabsch
        [ Rsm,tsm ] = kabsch( Xm , Xs ); % defined to act as tensor over ms
        
        % write
        Rvw{s,m}=Rsm;
        tvw{s,m}=tsm;
        
        % plot
%         Xp=Rsm*Xm+tsm;
%         figure; plot3(Xs(1,:),Xs(2,:),Xs(3,:),'.'); hold on;  plot3(Xp(1,:),Xp(2,:),Xp(3,:),'.'); axis equal
%         zlim([0 25])
%         xlim([-10 10])
%         ylim([-10 10])
        drawnow
    end 
end 

%% fix view 1 and define camera translation
% initiate minimization
c=unique(qcb(3,:));
A=cell(length(c),length(c));
b=cell(length(c),length(c));

% loop true permutation of cameras
for m=c
    for s=c(c~=m)
        
        A{m,s}=sparse(3,3*length(c));
        b{m,s}=sparse(3,1);
        
        Rsm=Rvw{s,m};
        tsm=tvw{s,m};
        
        A{m,s}(:,1+3*(m-1):3*m)=-Rsm;
        A{m,s}(:,1+3*(s-1):3*s)=eye(3);
                
        b{m,s}=tsm;
        
    end
end

% assemle minization
A=cat(1,A{:});
b=cat(1,b{:});

% define the base system by strict equlatiy constr.
C=speye(3,size(A,2));
d=sparse(3,1);

% solve problem
x=lsqlin(A,b,[],[],C,d,[],[],[],opt);

% write translations
tvec=cell(1,length(c));
for m=c
    tvec{m}=x(1+3*(m-1):3*m);
end

%% fix view 1 and define camera rotations
% initiate the linear minimization
c=unique(qcb(3,:));
A=cell(length(c),length(c));

% loop over cameras
for m=c
    
    for s=c(c~=m)
        
        A{m,s}=sparse(9,9*length(c));
        A{m,s}(:,1+9*(m-1):9*m)=speye(9);
    
        Rms=Rvw{m,s}; % to be checked
                
        A{m,s}(:,1+9*(s-1):9*s)=-kron(eye(3),Rms);
    end
    
end
A=cat(1,A{:});
b=sparse(size(A,1),1);

% define constraint that fix global frame to first view
C=speye(9,9*length(c));
d=reshape(speye(3),[],1);

% solve problem
x=lsqlin(A,b,[],[],C,d,[],[],[],opt); % impose something kabsch like?

% write translations
Rmat=cell(1,length(c));
for m=c
    Q=reshape(x(1+9*(m-1):9*m),3,3);
    Rmat{m}=frotm(Q);
end

%% averaged consistent triangulation by averaging proc
for f=unique(qcb(1,:))
    % file
    F=qcb(1,:)==f;
    
    for n=unique(qcb(2,F))
        % frame
        N=qcb(2,:)==n;
        
        % camera averaged input data, not using any quaternion rules for now
        R=zeros(3,3);
        t=zeros(3,1);
        c=unique(qcb(3,F&N));
        for m=c
            M=qcb(3,:)==m;
            R=R+Rmat{m}\qcom2rotm(qcb(4:end,F&N&M));
            t=t+Rmat{m}\(tcb(4:end,F&N&M)-tvec{m});
        end
        R=frotm(R/length(c));
        t=t/length(c);
        q=rotm2qcom(R);
        
        for m=c
            M=qcb(3,:)==m; % some nasty fix to normalization
            qcb(4:end,F&N&M)=rotm2qcom(Rmat{m}*qcom2rotm(q))/norm(rotm2qcom(Rmat{m}*qcom2rotm(q)),2);
            tcb(4:end,F&N&M)=Rmat{m}*t+tvec{m};
            
        end
        
    end
end

%% create objective function fit parametrized homography
% image and checkerboard defined coordinates
x=sym('x',[2 1]);
X=sym('X',[2 1]);

% checkerboard rotation and translation
q=sym('q',[4,1]);
t=sym('t',[3 1]);
r=qcom2rotm(q);
X=[r,t]*[X;0;1];

% parametrize homography up to fixed scale
Q=sym('Q',[4,1]);
R=qcom2rotm(Q);
T=sym('T',[3,1]);
K=sym('K',[5,1]);
K=[K(1) K(2) K(3)
    0 K(4) K(5)
    0 0 1];
P=K*[R T]; %#a9=1, scale is fixed: les unknown

% projected coordinate
xp=subs(subs(sym(Pmap),sym('X',[4 1]),[X;1]),sym('P',[3 4]),P);

% define argument objective function
ObjA=x-xp;

% local projection area
J=difobjf(xp,sym('X',[2 1]));
ObjW=1/(ctrl.tsiz^2*abs(J(1,1)*J(2,2)-J(1,2)*J(2,1))); % local area by jacobian

% constraints
eqcF=[dot(q,q)-1
    dot(Q,Q)-1];

% define objective function
ObjF=defobjf(ObjA,eqcF);

% variables
vars=sort(symvar(sym(ObjF)));

% write function
objf=matlabFunction(ObjF,'file',[folder date cal vsl 'mfunc' vsl 'objf.m'],'Vars',vars);%,'Optimize',false);
wobjf=matlabFunction(ObjW,'file',[folder date cal vsl 'mfunc' vsl 'wobjf.m'],'Vars',vars);%,'Optimize',false);

%% sequential geometric refiment step
% numbering knowns
nk=[1:5+4+3+2,length(vars)-3,length(vars)-2,length(vars)]';

% numbering unknowns
nu=(1:length(vars))';
nu(nk)=[];

% find grad and hess w.r.t. unknowns
[gObjF,hObjF]=difobjf(ObjF,vars(nu));
gobjf=matlabFunction(gObjF,'file',[folder date cal vsl 'mfunc' vsl 'gobjf.m' ],'Optimize',false,'Vars',vars);%,'Optimize',false);
hobjf=matlabFunction(hObjF,'file',[folder date cal vsl 'mfunc' vsl 'hobjf.m' ],'Optimize',false,'Vars',vars);

% init
rpe(9:10,:)=zeros(2,size(rpe,2)); %11:12
for f=unique(qcb(1,:))
    % file
    F1=qcb(1,:)==f;
    F2=IMcnod(1,:)==f;
    
    for n=unique(qcb(2,F1))
        % frame
        N1=qcb(2,:)==n;
        N2=IMcnod(2,:)==n;
                
        % checkerboard from first view dewarped to world 
        % qcb(3,N1&F1)
        c=min(qcb(3,N1&F1));
        C=qcb(3,:)==c;
        q=rotm2qcom(Rmat{c}\qcom2rotm(qcb(4:end,N1&F1&C)));
        t=Rmat{c}\(tcb(4:end,N1&F1&C)-tvec{c});
        
        % camera coordinates
        x=(IMcnod(5:6,F2&N2));
        
        %checkerboard
        X=repmat(chkboard.points,1,size(x,2)/prod(ctrl.nnod));
        
        % camera sorting
        c=IMcnod(3,F2&N2);
        
        % quaternion and camera parameters
        T=cell2mat(tvec(c));
        Q=cell2mat(cellfun(@(x)rotm2qcom(x),Rmat(c),'UniformOutput',false));
        K=cell2mat(cellfun(@(x)kmat2cpar(x)',Kmat(c),'UniformOutput',false));
        
        % dewarped coordinates
        b=cell2mat(cellfun(@(D)D.map',Dpar(c),'UniformOutput',false));
        inp=num2cell([b;x],2);
        x=Dmap(inp{:});% distortion mapping
        x=cell2mat(cellfun(@(D,x)homc2inhc(D.H*inhc2homc(x)),Dpar(c),num2cell(x,1),'UniformOutput',false)); % back to image resolution
        
        %initial error
        inp=num2cell([K;Q;T;X;repmat(q,1,size(X,2));repmat(t,1,size(X,2));x;zeros(2,size(x,2))],2);
        rpe(9,F2&N2)=sqrt(objf(inp{:})); % err
        rpe(10,F2&N2)=sqrt(objf(inp{:}).*wobjf(inp{:})); % norm err
        
        % numering assembly w.r.t solution vector
        asn=repmat((1:length(nu))',1,size(x,2)); % assembly numbering
        
        % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
        dat=[zeros(size(nk)) nk [K;Q;T;X;x;zeros(1,size(x,2))]
            ones(size(nu)) nu asn];
        
        % initiate solution vector
        x0=[q
            t
            0]; % or do 1 ord coef a
        
        % solve global optimization
        xsol = optobjf([], objf, gobjf, hobjf , x0 , dat , ctrl.optl ,'newton','no');%wobjf
%         xsol = optobjf( ObjF, gObjF, hObjF , x0 , dat , 10^(-3) ,'impdesc'); 
        
        % write
        q=xsol(1:4)/norm(xsol(1:4),2);
        t=xsol(5:7);
        for c=(qcb(3,F1&N1))
            C=qcb(3,:)==c;
            qcb(4:end,F1&N1&C)=rotm2qcom(Rmat{c}*qcom2rotm(q));
            tcb(4:end,F1&N1&C)=Rmat{c}*t+tvec{c};
            
        end
        
        % final error
%         inp=num2cell([K;Q;T;X;repmat(q,1,size(X,2));repmat(t,1,size(X,2));x;zeros(2,size(x,2))],2);
%         rsq(12,F2&N2)=objf(inp{:});
        
    end
    
    % message
    disp(['fit chkb. fil. ',num2str(f)])
    
end

%% update camera matrix and positioning
% numbering knowns
nk=(1+5+4+3:length(vars)-1)';

% numbering unknowns
nu=(1:length(vars))';
nu(nk)=[];

% find grad and hess. w.r.t opt. rot m
[gObjF,hObjF]=difobjf(ObjF,vars(nu));%can upgrade speed here, by only diff to need
gobjf=matlabFunction(gObjF,'file',[folder date cal vsl 'mfunc' vsl 'gobjf.m' ],'Optimize',false,'Vars',vars);%,'Optimize',false);
hobjf=matlabFunction(hObjF,'file',[folder date cal vsl 'mfunc' vsl 'hobjf.m' ],'Optimize',false,'Vars',vars);

% init
% rsq(11:12,:)=zeros(2,size(rsq,2));
for c=unique(qcb(3,:))
    % indexing
    C1=qcb(3,:)==c;
    C2=IMcnod(3,:)==c;
    
    % camera stuff
    K=kmat2cpar(Kmat{c})';
    Q=rotm2qcom(Rmat{c}); 
    T=tvec{c}; 
    
    % camera coordinates
    x=(IMcnod(5:6,C2));
    
    % dewarped coordinates
    inp=num2cell([repmat(Dpar{c}.map',1,size(x,2))
        x],2);
    x=Dmap(inp{:}); % distortion mapping
    x=homc2inhc(Dpar{c}.H*inhc2homc(x)); % back to image resolution
    
    % checkerboard
    X=repmat(chkboard.points,1,size(x,2)/prod(ctrl.nnod));
    
    % rgid body motion
    q=qcb(4:end,C1);
    t=tcb(4:end,C1);
    q=cellfun(@(x)rotm2qcom(Rmat{c}\qcom2rotm(x)),num2cell(q,1),'UniformOutput',false);
    qq=cell2mat(cellfun(@(x)repmat(x,1,prod(ctrl.nnod)),q,'UniformOutput',false));
    t=cellfun(@(x)Rmat{c}\(x-tvec{c}),num2cell(t,1),'UniformOutput',false);
    tt=cell2mat(cellfun(@(x)repmat(x,1,prod(ctrl.nnod)),t,'UniformOutput',false));
    
    % init error
%     inp=num2cell([repmat([K;Q;T],1,size(X,2));X;qq;tt;x;zeros(2,size(x,2))],2);
%     rsq(13,C2)=objf(inp{:});
    
    % numering assembly w.r.t solution vector
    asn=repmat((1:length(nu))',1,size(x,2)); % assembly numbering
    
    % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
    dat=[zeros(size(nk)) nk [X;qq;tt;x;zeros(1,size(x,2))]
        ones(size(nu)) nu asn];
    
    % initiate solution vector
    x0=[K
        Q
        T
        0]; % or do 1 ord coef a
    
    % solve global optimization
    xsol = optobjf([], objf, gobjf, hobjf , x0 , dat , ctrl.optl ,'newton');%wobjf
    
    % write
    K=xsol(1:5);
    Q=xsol(6:9)/norm(xsol(6:9),2);
    T=xsol(10:12);
    
    Kmat{c}=cpar2kmat(K(1),K(2),K(3),K(4),K(5));
    Rmat{c}=qcom2rotm(Q);
    tvec{c}=T;
    
    % final error
%     inp=num2cell([repmat([K;Q;T],1,size(X,2));X;qq;tt;x;zeros(2,size(x,2))],2);
%     rsq(14,C2)=objf(inp{:});
    
    % q and t also changed now!
    q=cellfun(@(x)rotm2qcom(Rmat{c}*qcom2rotm(x)),q,'UniformOutput',false); 
    t=cellfun(@(x)Rmat{c}*x+tvec{c},t,'UniformOutput',false);    
    qcb(4:end,C1)=cell2mat(q);
    tcb(4:end,C1)=cell2mat(t);
    
end

%% full problem to optimize excl distmap
% numbering knowns
nk=[13:14,length(vars)-3,length(vars)-2]';

% numbering unknowns
nu=(1:length(vars))';
nu(nk)=[];

% find grad and hess. w.r.t opt. rot m
[gObjF,hObjF]=difobjf(ObjF,vars(nu));
gobjf=matlabFunction(gObjF,'file',[folder date cal vsl 'mfunc' vsl 'gobjf.m' ],'Optimize',false,'Vars',vars);%,'Optimize',false);
hobjf=matlabFunction(hObjF,'file',[folder date cal vsl 'mfunc' vsl 'hobjf.m' ],'Optimize',false,'Vars',vars);

% init final reproj err
rpe(11:12,:)=zeros(2,size(rpe,2));

% get each board
[~,ia,ic]=unique(qcb(1:2,:)','rows','stable'); % unique board
q=num2cell(qcb(3:end,ia),1);
q=cell2mat(cellfun(@(x)rotm2qcom(Rmat{x(1)}\qcom2rotm(x(2:5))),q,'UniformOutput',false));
t=num2cell(tcb(3:end,ia),1);
t=cell2mat(cellfun(@(x)Rmat{x(1)}\(x(2:4)-tvec{x(1)}),t,'UniformOutput',false));
z=zeros(1,size(t,2));

% solution vect
K=cellfun(@(x)kmat2cpar(x)',Kmat,'UniformOutput',false);
Q=cellfun(@(x)rotm2qcom(x),Rmat,'UniformOutput',false);
T=tvec;
Z=num2cell(zeros(1,size(T,2)),1);

% rgid body motion
x0=[reshape(cell2mat([K;Q;T;Z]),[],1)
    reshape([q;t;z],[],1)];

% camera coordinates
x=(IMcnod(5:6,:));
c=(IMcnod(3,:));

% checkerboard
X=repmat(chkboard.points,1,size(x,2)/prod(ctrl.nnod));

% dewarped coordinates
b=cell2mat(cellfun(@(D)D.map',Dpar(c),'UniformOutput',false));
inp=num2cell([b;x],2);
x=Dmap(inp{:});% distortion mapping
x=cell2mat(cellfun(@(D,x)homc2inhc(D.H*inhc2homc(x)),Dpar(c),num2cell(x,1),'UniformOutput',false)); % back to image resolution

% init error
% inp=num2cell([cell2mat(cellfun(@(x)kmat2cpar(x)',Kmat(c),'UniformOutput',false))
%     cell2mat(cellfun(@(x)rotm2qcom(x),Rmat(c),'UniformOutput',false))
%     cell2mat(tvec(c));X
%     reshape(repmat(q(:,ic),prod(ctrl.nnod),1),4,[])
%     reshape(repmat(t(:,ic),prod(ctrl.nnod),1),3,[]);x;zeros(2,size(x,2))],2);
% rsq(15,:)=objf(inp{:});

% assembly numbering
ASN=num2cell(reshape(1:13*prop.res(4),[],prop.res(4)),1);
asn=reshape(13*prop.res(4)+1:length(x0),8,[]);
asn=asn(:,ic);
num=[cell2mat(cellfun(@(x)x(1:end-1),ASN(c),'UniformOutput',false))
    reshape(repmat(asn,prod(ctrl.nnod),1),8,[])
    cell2mat(cellfun(@(x)x(end),ASN(c),'UniformOutput',false))];

% data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
dat=[zeros(size(nk)) nk [X;x]
    ones(size(nu)) nu num];

% solve global optimization
xsol = optobjf([], objf, gobjf, hobjf , x0 , dat , ctrl.optl ,'newton');%wobjf
xsol = optobjf(wobjf, objf, gobjf, hobjf , xsol , dat , ctrl.optl ,'newton');

% write camera stuff
Xdat=reshape(xsol(1:13*prop.res(4)),13,[]);
for c=1:prop.res(4)
    Kmat{c}=cpar2kmat(Xdat(1,c),Xdat(2,c),Xdat(3,c),Xdat(4,c),Xdat(5,c));
    Rmat{c}=qcom2rotm(Xdat(6:9,c)/norm(Xdat(6:9,c),2));
    tvec{c}=Xdat(10:12,c);
end

% write triang in cam frm
xdat=reshape(xsol(13*prop.res(4)+1:end),8,[]);
c=qcb(3,:);
xdat=xdat(:,ic);
q=num2cell([c;xdat(1:4,:)./sqrt(sum(xdat(1:4,:).^2,1))],1);
q=cell2mat(cellfun(@(x)rotm2qcom(Rmat{x(1)}*qcom2rotm(x(2:5))),q,'UniformOutput',false));
t=num2cell([c;xdat(5:7,:)],1);
t=cell2mat(cellfun(@(x)Rmat{x(1)}*x(2:4)+tvec{x(1)},t,'UniformOutput',false));
qcb(4:end,:)=q;
tcb(4:end,:)=t;

% end error
c=(IMcnod(3,:));
inp=num2cell([cell2mat(cellfun(@(x)kmat2cpar(x)',Kmat(c),'UniformOutput',false))
    cell2mat(cellfun(@(x)rotm2qcom(x),Rmat(c),'UniformOutput',false))
    cell2mat(tvec(c));X
    reshape(repmat(xdat(1:4,:)./sqrt(sum(xdat(1:4,:).^2,1)),prod(ctrl.nnod),1),4,[])
    reshape(repmat(xdat(5:7,:),prod(ctrl.nnod),1),3,[]);x;zeros(2,size(x,2))],2);
rpe(11,:)=sqrt(objf(inp{:}));
rpe(12,:)=sqrt(objf(inp{:}).*wobjf(inp{:}));

%% define global coordinate system
% triangulate camera centers
xcen=cell(1,prop.res(4));
Xcen=cell(1,prop.res(4));
for c=1:prop.res(4)
    xcen{c}=-Rmat{c}\tvec{c};
    Xcen{c}=Rmat{c}\eye(3);
end

% define global coordinates system
A=cat(2,xcen{:});
tgl=mean(A,2);

% get global transform
B=cat(3,Xcen{:});
D=frotm(mean(B,3));
[Rgl,~,~]=svd(A-tgl); % unitary matrix aligning to point cloud
Rgl=Rgl.*sign(diag(Rgl\D)'); % correction to look same avg direction

% rotate to new frame
for c=1:prop.res(4)
    
    tvec{c}=tvec{c}+Rmat{c}*tgl;
    Rmat{c}=Rmat{c}*Rgl;
    
    
end

% interchange x=x z=y y=k
Q=[1 0 0
    0 0 1
    0 1 0];
% Rgl=Q*Rgl;

% flip acis to new frame
Pmat=cell(1,prop.res(4));
for c=1:prop.res(4)
    
    Rmat{c}=Rmat{c}*Q;
    
    xcen{c}=-Rmat{c}\tvec{c};
    
    Pmat{c}=krtm2pmat(Kmat{c},Rmat{c},tvec{c});
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
save([folder date cal vsl 'Rmat.mat'],'Rmat')
save([folder date cal vsl 'tvec.mat'],'tvec')
% save([folder date cal '\Pmat.mat'],'Pmat') % dont need
% save([folder date cal '\xcen.mat'],'xcen')

%% Plot calibration results
% checkerboard points
X=inhc2homc(chkboard.points,0);

% some figure
figure
hold on

% global coordinates
quiver3(0,0,0,1,0,0,0,'k.')
quiver3(0,0,0,0,1,0,0,'k.')
quiver3(0,0,0,0,0,1,0,'k.')

% camera coordinates
for c=1:length(xcen)
    text(xcen{c}(1),xcen{c}(2),xcen{c}(3),num2str(c))
    Xcen=Rmat{c}*eye(3);
    quiver3(xcen{c}(1),xcen{c}(2),xcen{c}(3),Xcen(1,1),Xcen(1,2),Xcen(1,3),0,'r.')
    quiver3(xcen{c}(1),xcen{c}(2),xcen{c}(3),Xcen(2,1),Xcen(2,2),Xcen(2,3),0,'g.')
    quiver3(xcen{c}(1),xcen{c}(2),xcen{c}(3),Xcen(3,1),Xcen(3,2),Xcen(3,3),0,'b.')
    
end

q=num2cell(qcb(3:end,:),1);
R=cell2mat(cellfun(@(x)Rmat{x(1)}\qcom2rotm(x(2:5)),q,'UniformOutput',false)');
t=num2cell(tcb(3:end,:),1);
t=cell2mat(cellfun(@(x)Rmat{x(1)}\(x(2:4)-tvec{x(1)}),t,'UniformOutput',false)');
X=reshape(R*X+t,3,[]);%5,[]);

% plot checkerboards
plot3(X(1,:),X(2,:),X(3,:),'.')
 
% axis
axis tight
axis equal
grid minor
box on
xlabel('x $[\rm{m}]$','interpreter','latex')
ylabel('y $[\rm{m}]$','interpreter','latex')
zlabel('z $[\rm{m}]$','interpreter','latex')
view(3)
title('Multiple view calibration','interpreter','latex')
drawnow

end