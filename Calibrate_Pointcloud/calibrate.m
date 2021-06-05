function [Dpar,Kmat,Rmat,tvec,camprop]=calibrate
%calibrate Perform the calibration from the porcessed .xml file
%   In sequential order: fit projective mapping from world to camera plane
%   correct distortion by polynomial mapping to user define order.
%
%   Notes,
%       - Input: caliibration data processed into binary from .xml file
%       Davis
%       - Output,
%           Kmat (upper triangular intrinsic camera matrix)
%           Rmat (extrinsic rotationmatrix w.r.t. world points)
%           tvec ( extrinsictranslation vector
%           xcen ( (virtual-)camera centers)
%
% NOTE that this calibration goes in the wrong order calibration dewarping...
%   RIGHT would be to first calibrate the distortion by colinear points
%   and then perform the camera calibration
%   when affine imaging this is ok, the distortion does not vary much over
%   depth of field and therefore can be corrected for afterwards

%% Get global variables
global folder date cal Pmap Dmap ctrl

%% Load calibration marks
Cmark=fload([folder date cal vsl 'Cmark.dat']);

%% Define objective function
% world coordinates
X=sym('X',[3 1]);

% image coordinates
x=sym('x',[2 1]);

% parametrized projection matrix
p=sym('p',[11,1]);
P=[p(1) p(2) p(3) p(4)
    p(5) p(6) p(7) p(8)
    p(9) p(10) p(11) 1];

% substitute parametrization
xp=subs(sym(Pmap),sym('P',[3 4]),P);
xp=subfunv(xp,sym('X',[4 1]),[X;1]);
% xd=subfunv(Dmap,sym('b',[1 6]),[0 0 1 0 0 1]); % only include dist coef.
xd=subs(sym(Dmap),sym('x',[2 1]),x);
bvar=sort(symvar(sym(Dmap)));
bvar(end-1:end)=[];

% objective argument
ObjA=xd-xp;

% define objective function
ObjF=defobjf(ObjA,[]);

% variabes
vars=sort(symvar(sym(ObjF)));

objf=matlabFunction(ObjF,'Optimize',false,'Vars',vars);

%% Compute best projective fit
% numbering knowns
nk=[1:3+length(bvar),15+length(bvar):length(vars)]';

% numbering unknowns
nu=(1:length(vars))';
nu(nk)=[];

% find grad and hess w.r.t. unknowns
[gObjF,hObjF]=difobjf(ObjF,vars(nu));

gobjf=matlabFunction(gObjF,'Optimize',false,'Vars',vars);

hobjf=matlabFunction(hObjF,'Optimize',false,'Vars',vars);

%% initiate projective matrix and camera mapping
Dpar=cell(1,max(Cmark(1,:)));
Pmat=cell(1,max(Cmark(1,:)));
Kmat=cell(1,max(Cmark(1,:)));
Rmat=cell(1,max(Cmark(1,:)));
tvec=cell(1,max(Cmark(1,:)));
for c=unique(Cmark(1,:))
    % display
    display(['projective mapping camera ',num2str(c)])
    
    % define camera index
    C=Cmark(1,:)==c ;% & Cmark(6,:)==c;
        
    % get image data in homogeneous format
    x=inhc2homc(Cmark(2:3,C));
    X=inhc2homc(Cmark(4:6,C));
    
    % initiate projection matrix
    P=dirlintrans(X,x); % direct
    P=P/P(end,end); % fix scale
    p=reshape(P',[],1); % get entries
    p(end)=[]; % remove fixed scale
    
    % convert to homogenous coordinates for regression
    x=homc2inhc(x);
    X=homc2inhc(X);
    
    % initiate polynomial distortion correction
    b=[0 0 1 0 0 1 zeros(1,length(nk)-6-size(X,1)-size(x,1))]';
    Dpar{c}=b';
    b=repmat(b,1,size(X,2));
    
    % solve objective function
    asn=repmat((1:length(nu))',1,size(X,2)); % assembly numbering
    
    % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
    dat=[zeros(size(nk)) nk [X;b;x]
        ones(size(nu)) nu asn];
    
    % initiate solution vector
    x0=p; % or do 1 ord coef a
    
    % solution
    xsol = optobjf( [], objf, gobjf, hobjf , x0 , dat , ctrl.optl ,'newton');
    
    % Pmat
    P=reshape([xsol;1],4,3)';
    
    % write output
    Pmat{c}=P;
    
    % decompose pmat
    [ K,R,t  ] = pmat2krtm( P );
    Kmat{c}=K;
    Rmat{c}=R;
    tvec{c}=t;
end

%% Correct for distortion
% numbering knowns
nk=[1:3,4+length(bvar):length(vars)]';

% numbering unknowns
nu=(1:length(vars))';
nu(nk)=[];

% find grad and hess w.r.t. unknowns
[gObjF,hObjF]=difobjf(ObjF,vars(nu));

gobjf=matlabFunction(gObjF,'Optimize',false,'Vars',vars);

hobjf=matlabFunction(hObjF,'Optimize',false,'Vars',vars);

%% initiate camera mapping
camprop=struct([]);
for c=unique(Cmark(1,:))
    % display
    display(['distortion correction camera ',num2str(c)])
    tic
    
    % define camera indexing
    C=Cmark(1,:)==c ;% & Cmark(6,:)==c;
        
    % get image data in homogeneous format
    x=Cmark(2:3,C);
    X=Cmark(4:6,C);
    
    % projective matrix
    p=repmat(reshape(Pmat{c}',[],1),1,size(X,2));
    p(end,:)=[];
    
    % solve objective function
    asn=repmat((1:length(nu))',1,size(X,2)); % assembly numbering
    
    % data assembly [unknown set indicator | numering w.r.t. input | [data -- assembly] ]
    dat=[zeros(size(nk)) nk [X;p;x]
        ones(size(nu)) nu asn];
    
    % initiate solution vector
    b=Dpar{c}';
    x0=b;
    
    % solution
    [xsol,res] = optobjf( [], objf, gobjf, hobjf , x0 , dat , ctrl.optl ,'newton');
    
    % write output
    Dpar{c}=xsol';
    
    % write a simple reproj error
    camprop(c).reprjavg=sqrt(res/(size(dat,2)-2));
    camprop(c).reprjstd=sqrt(res/(size(dat,2)-2));
end

%% save
save([folder date cal vsl 'Dpar.mat'],'Dpar')
save([folder date cal vsl 'Kmat.mat'],'Kmat')
save([folder date cal vsl 'Rmat.mat'],'Rmat')
save([folder date cal vsl 'tvec.mat'],'tvec')
save([folder date cal vsl 'camprop.mat'],'camprop')

end