function postproc_rotationcenterframe
%postproc_rotationcenterframe Process frame to center of rotation
%   This is a nondynamic reference frame mapping to the instateneous center
%   of rotation.
%
%   This function maps data to the center of rotation of the group
%
%   This is the simplest instateneous [eulerian] reference to the group.
%
%   Todo: Material accelaration fluctuation

%% get globals
global folder date rec post

%% Get data

% trajectory data
TrajIndex=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Time=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Time.dat']);
Quadric=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Quadric.dat']);
Pos=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
Vel=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
% Acc=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accelaration.dat']);

% group data
GroupIndex=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupIndex.dat']); % Group Index
SegIndex=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'SegmentationIndex.dat']); % Group Index
CenRot=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'CenterRotation.dat']); % Center Rotation
GroupRotation=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupRotation.dat']); % Center Rotation
GroupStrainRate=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupStrainRate.dat']); % Center Rotation
% GroupAccelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupAccelaration.dat']); % Center Rotation

%% create path for groupframe if not there
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'GroupFrame'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'GroupFrame'])
end

%% create path for processing co translate frame if not there
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'GroupFrame' vsl 'RotationCenterFrame'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'GroupFrame' vsl 'RotationCenterFrame'])
end

%% initiate group frame
grpfrm_index=zeros(2,0); % track indexing
grpfrm_time=zeros(1,0); % time
grpfrm_quad=zeros(10,0); % track quadrics
grpfrm_pos=zeros(3,0); % track positions
grpfrm_vel=zeros(3,0); % track velocity 
grpfrm_velflucrot=zeros(3,0); % velocity fluctuation to group rotation
grpfrm_velflucdisp=zeros(3,0); % velocity fluctuation to displacement group
% grpfrm_acc=zeros(3,0); % track accelaration
% grpfrm_accfluc=zeros(3,0); % accelaration fluctuation to [material] accelaration

%% Map data to moving (co-translating) reference

disp('post proc')
tic

% loop over time frames in the data
for n=post.tproc(1):post.tproc(2)
    
    % Data Segmentation
    G=GroupIndex(2,:)==n;
    
    % Center Rotation
    cenrot=CenRot(:,G);
    
    % Rotation motion
    W=GroupRotation(:,G);
    W=[  0    -W(3)  W(2)
        W(3) 0     -W(1)
        -W(2) W(1) 0   ];
    
    % Strain motion
    F=GroupStrainRate(:,G);
    F=[ F(1)   F(6)/2 F(5)/2
        F(6)/2 F(2)   F(4)/2
        F(5)/2 F(4)/2 F(3)  ];
    
    % homogeneous rigid body motion
    H=inv([eye(3) cenrot ; zeros(1,3) 1]);
    
    % timespan
    N=TrajIndex(2,:)==n & ismember(TrajIndex(1,:),SegIndex(2,:)) ;
    
    % get data in frame span trajectory
    tind=TrajIndex(:, N); % indexing data
    time=Time(:, N); % indexing data
    quad=Quadric(:, N); % quadric data
    pos=Pos(:, N); % position data
    vel=Vel(:, N); % velocity data
%     acc=Acc(:, N); % accelaration data
    
    %figure; plot3(pos(1,:),pos(2,:),pos(3,:),'.'); axis equal; view(2)
    
    % transform data to new reference // this could be done first decompos dRmat and ddRmat if you like
    pos= pos - cenrot ; % position in moving reference
    
    % rotation fluctuation
    velflucrot=vel-W*pos;
    
    % strain motion fluctuation
    velflucdisp=vel-(W+F)*pos;
    
%     % decomposed relative accelaration in frame of reference
%     relacc= acc - ddtvec ; % difference in accelaration
    
    % transform quadric to new reference
    quad=quadtrans(quad,H); % quadric
    
    %figure; plot3(pos(1,:),pos(2,:),pos(3,:),'.'); axis equal; view(2)
    
    %%% write variables beneath
    
    % trajectory frame indexing
    grpfrm_index=cat(2,grpfrm_index,tind);
    
    % write time
    grpfrm_time=cat(2,grpfrm_time,time);
    
    % perform coordinate transformation quadric
    grpfrm_quad=cat(2,grpfrm_quad,quad); % to be checked
    
    % perform coordinate transformation position
    grpfrm_pos=cat(2,grpfrm_pos,pos);
    
    % perform coordinate transformation velocity
    grpfrm_vel=cat(2,grpfrm_vel,vel);
    
    % write velocity rotation fluctuation
    grpfrm_velflucrot=cat(2,grpfrm_velflucrot,velflucrot); % velocity fluctuation to group rotation
    
    % write velocity full displacement field fluctuation
    grpfrm_velflucdisp=cat(2,grpfrm_velflucdisp,velflucdisp); % velocity fluctuation to displacement group

%     % perform coordinate transformation accelaration
%     grpfrm_acc=cat(2,grpfrm_acc,acc);
%     
%     % relative accelaration
%     grpfrm_relacc=cat(2,grpfrm_relacc,relacc);
    
    % message
    disp(['Processed center rotation groupframe at frame ',num2str(n)])
    
end

%% save
% Indexing
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'RotationCenterFrame' vsl 'Index.dat']...
    ,grpfrm_index,'w');

% Time
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'RotationCenterFrame' vsl 'Time.dat']...
    ,grpfrm_time,'w');

% Quadric
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'RotationCenterFrame' vsl 'Quadric.dat']...
    ,grpfrm_quad,'w');

% Position
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'RotationCenterFrame' vsl 'Position.dat']...
    ,grpfrm_pos,'w');

% Velocity
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'RotationCenterFrame' vsl 'Velocity.dat']...
    ,grpfrm_vel,'w');

% Velocity fluctuation rotation
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'RotationCenterFrame' vsl 'VelocityFluctuationRotation.dat']...
    ,grpfrm_velflucrot,'w');

% Velocity fluctuation displacement
fsave([folder date rec vsl ...
    'Postproc3DTracking' vsl 'GroupFrame' vsl 'RotationCenterFrame' vsl 'VelocityFluctuationDisplacement.dat']...
    ,grpfrm_velflucdisp,'w');

% % Accelaration
% fsave([folder date rec vsl ...
%     'Postproc3DTracking' vsl 'GroupFrame' vsl 'RotationCenterFrame' vsl 'Accelaration.dat']...
%     ,grpfrm_acc,'w');
% 
% % Relative Accelaration
% fsave([folder date rec vsl ...
%     'Postproc3DTracking' vsl 'GroupFrame' vsl 'RotationCenterFrame' vsl 'RelativeAccelaration.dat']...
%     ,grpfrm_relacc,'w');

end

