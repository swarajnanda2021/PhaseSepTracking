function postproc_cotranslatingframe
%postproc_cotranslateframe Process frame to center of gravity points cloud
%   This is a dynamic reference frame mapping to the translating center
%   of gravity.
%
%   This function maps data to the center of gravity of the group
%
%   This is the simplest lagrangian reference to the group.

%% get globals
global folder date rec prop post grps

%% Get data

% trajectory data
TrajIndex=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Time=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Time.dat']);
Quadric=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Quadric.dat']);
Pos=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
Vel=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
Acc=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accelaration.dat']);

% group data
GroupIndex=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupIndex.dat']); % Group Index
SegIndex=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'SegmentationIndex.dat']); % Group Index
CenGrav=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'CenterGravity.dat']); % Center Rotation

%% create path for groupframe if not there
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'GroupFrame'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'GroupFrame'])
end

%% make a symbolic problem for moving reference frame

% define curve translation
ot=(0:grps.tord)';
t=sym('t');
T=reshape(sym('b',[1 3*length(ot)]),3,[])*t.^ot;

% Time derivatives
dT=diff(T,t);
ddT=diff(dT,t);

% variabes
vars=sort(symvar(T));
vars=vars(1:end-1);

% substitute t=0
T=subs(T,t,0);
dT=subs(dT,t,0);
ddT=subs(ddT,t,0);

% make anonymous functions
tfun=matlabFunction(T,'Vars',vars); % ,'Optimize',false);
dtfun=matlabFunction(dT,'Vars',vars); % ,'Optimize',false);
ddtfun=matlabFunction(ddT,'Vars',vars); % ,'Optimize',false);

%% create path for processing co translate frame if not there
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'GroupFrame' vsl 'CoTranslatingFrame'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'GroupFrame' vsl 'CoTranslatingFrame'])
end

%% initiate group frame
grpfrm_index=zeros(2,0); % track indexing
grpfrm_time=zeros(1,0); % track quadrics
grpfrm_quad=zeros(10,0); % track quadrics
grpfrm_pos=zeros(3,0); % track positions
grpfrm_vel=zeros(3,0); % track velocity stil reference
grpfrm_relvel=zeros(3,0); % track velocity relative to moving reference
grpfrm_acc=zeros(3,0); % track accelaration
grpfrm_relacc=zeros(3,0); % relative in moving frame

%% Map data to moving (co-translating) reference

disp('post proc')
tic

% loop over time frames in the data
for n=post.tproc(1):post.tproc(2)
    
    % Data Segmentation
    G=GroupIndex(2,:)>=n-floor(grps.tres/2) & GroupIndex(2,:)<=n+floor(grps.tres/2);
    
    % group index
    gind=GroupIndex(:,G);
    
    % Center Rotation
    cengrav=CenGrav(:,G);
    
    % fit spline
    coef=polytraj([gind ; cengrav],grps.tord);
    
    % curve data for moving reference
    curve=num2cell(reshape(coef(3:5,:),[],1)',1);
    
    % evaluate moving frame of reference
    tvec=tfun(curve{:}); % translation vector
    dtvec=dtfun(curve{:})*prop.fps; % time derivative translation vector
    ddtvec=ddtfun(curve{:})*prop.fps*prop.fps; % second time derivative translation vector
    
    % homogeneous rigid body motion
    H=inv([eye(3) tvec ; zeros(1,3) 1]);
    
    % timespan
    N=TrajIndex(2,:)==n & ismember(TrajIndex(1,:),SegIndex(2,:)) ;
    
    % get data in frame span trajectory
    tind=TrajIndex(:, N); % indexing data
    time=Time(:, N); % quadric data
    quad=Quadric(:, N); % quadric data
    pos=Pos(:, N); % position data
    vel=Vel(:, N); % velocity data
    acc=Acc(:, N); % accelaration data
    
    %figure; plot3(pos(1,:),pos(2,:),pos(3,:),'.'); axis equal; view(2)
    
    % decomposed relative velocity in frame of reference
    relvel= vel - dtvec ;% get velocity difference
    
    % decomposed relative accelaration in frame of reference
    relacc= acc - ddtvec ; % difference in accelaration
    
    % transform quadric to new reference
    quad=quadtrans(quad,H); % quadric
    
    % transform data to new reference // this could be done first decompos dRmat and ddRmat if you like
    pos= pos - tvec ; % position in moving reference
    
    %figure; plot3(pos(1,:),pos(2,:),pos(3,:),'.'); axis equal; view(2)
    
    %%% write variables beneath
    
    % trajectory frame indexing
    grpfrm_index=cat(2,grpfrm_index,tind);
    
    % trajectory frame indexing
    grpfrm_time=cat(2,grpfrm_time,time);
    
    % perform coordinate transformation quadric
    grpfrm_quad=cat(2,grpfrm_quad,quad); % to be checked
    
    % perform coordinate transformation position
    grpfrm_pos=cat(2,grpfrm_pos,pos);
    
    % perform coordinate transformation velocity
    grpfrm_vel=cat(2,grpfrm_vel,vel);
    
    % relative velocity in moving reference
    grpfrm_relvel=cat(2,grpfrm_relvel,relvel);
    
    % perform coordinate transformation accelaration
    grpfrm_acc=cat(2,grpfrm_acc,acc);
    
    % relative accelaration
    grpfrm_relacc=cat(2,grpfrm_relacc,relacc);
    
    % message
    disp(['processed cotranslating groupframe at frame ',num2str(n)])
    
end

%% save
% Indexing
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'CoTranslatingFrame' vsl 'Index.dat']...
    ,grpfrm_index,'w');

% Time
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'CoTranslatingFrame' vsl 'Time.dat']...
    ,grpfrm_time,'w');

% Quadric
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'CoTranslatingFrame' vsl 'Quadric.dat']...
    ,grpfrm_quad,'w');

% Position
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'CoTranslatingFrame' vsl 'Position.dat']...
    ,grpfrm_pos,'w');

% Velocity
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'CoTranslatingFrame' vsl 'Velocity.dat']...
    ,grpfrm_vel,'w');

% Relative Velocity
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'CoTranslatingFrame' vsl 'RelativeVelocity.dat']...
    ,grpfrm_relvel,'w');

% Accelaration
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'CoTranslatingFrame' vsl 'Accelaration.dat']...
    ,grpfrm_acc,'w');

% Relative Accelaration
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'CoTranslatingFrame' vsl 'RelativeAccelaration.dat']...
    ,grpfrm_relacc,'w');

end

