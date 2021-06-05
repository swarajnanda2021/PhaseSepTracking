function postproc_trajectframe
%postproc_trajectframe Process trajectory in frame w.r.t. xy-projected 
%   velocity heading and imposed z-gravity vector and copy all data in that 
%   frame.

global folder date rec prop post traj

%% Create memory maps
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Curve=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Curve.dat']);
Quadric=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Quadric.dat']);
Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
Accelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accelaration.dat']);

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame'])
end

%% compute function for transformation
% define curve
ot=(0:traj.tord)';
t=sym('t');
y=reshape(sym('b',[1 3*length(ot)]),3,[])*t.^ot;

% define base vectors from curve
tdif=diff(y,t);
zdir=sym([0 0 1]');
ydir=cross(tdif,zdir)/norm(cross(tdif,zdir),2);
xdir=cross(zdir,ydir);

% define transformations to curve reference
R=[xdir ydir zdir];
T=y;

dR=diff(R,t);
dT=diff(T,t);

ddR=diff(dR,t);
ddT=diff(dT,t);

% substitute t=0
R=subs(R,t,0);
T=subs(T,t,0);

dR=subs(dR,t,0);
dT=subs(dT,t,0);

ddR=subs(ddR,t,0);
ddT=subs(ddT,t,0);

% variabes
vars=sort(symvar(y));
vars=vars(1:end-1);

% make anonymous functions
Rfun=matlabFunction(R,'Vars',vars); % ,'Optimize',false);
tfun=matlabFunction(T,'Vars',vars); % ,'Optimize',false);

dRfun=matlabFunction(dR,'Vars',vars); % ,'Optimize',false);
dtfun=matlabFunction(dT,'Vars',vars); % ,'Optimize',false);

ddRfun=matlabFunction(ddR,'Vars',vars); % ,'Optimize',false);
ddtfun=matlabFunction(ddT,'Vars',vars); % ,'Optimize',false);

%% post processing
disp('post proc')
tic

% loop over the matched and indexed Plinks
for l=unique(Index(1,:))
    % get trajectory of interest
    L=Index(1,:)==l ;
    
    %%% vectorize beneath for speed %%%
    
    % initiate trajectory frame
    tfrm_index=zeros(2,0);
    tfrm_quad=zeros(10,0);
    tfrm_pos=zeros(3,0);
    tfrm_vel=zeros(3,0);
    tfrm_delvel=zeros(3,0);
    tfrm_rotvel=zeros(3,0);
    tfrm_relvel=zeros(3,0);
    tfrm_acc=zeros(3,0);
    tfrm_rotacc=zeros(3,0);
    tfrm_motacc=zeros(3,0);
    tfrm_delacc=zeros(3,0);
    tfrm_relacc=zeros(3,0);
    
    % loop timespan
    for n=unique(Index(2,L))
        
        % timespan
        N=ismember(Index(2,:),n);
        
        % get coefficient to local curve in right order
        curve=num2cell(reshape(reshape(Curve(:,L & N),[],3)',[],1),2);
        
        % evaluate moving frame of reference
        Rmat=Rfun(curve{:}); % rotation matrix
        tvec=tfun(curve{:}); % translation vector
        dRmat=dRfun(curve{:})*prop.fps; % time derivative rotation matrix
        dtvec=dtfun(curve{:})*prop.fps; % time derivative translation vector
        ddRmat=ddRfun(curve{:})*prop.fps*prop.fps; % second time derivative rotation matrix
        ddtvec=ddtfun(curve{:})*prop.fps*prop.fps; % second time derivative translation vector
        
        % homogeneous ridig body motion
        H=inv([Rmat tvec ; zeros(1,3) 1]);
        
        % get data in frame span trajectory
        index=Index(:, N); % indexing data
        quad=Quadric(:, N); % quadric data
        pos=Position(:, N); % position data
        vel=Velocity(:, N); % velocity data
        acc=Accelaration(:, N); % accelaration data
        
        % decomposed relative velocity in frame of reference
        delvel=( Rmat'*( vel - dtvec ) );% get velocity difference
        rotvel=( dRmat'*( pos - tvec ) );% get solid body rotational velocity induced by motion reference
        relvel=( rotvel + delvel ); % get total relative velocity to in reference
        
        % decomposed relative accelaration in frame of reference
        rotacc=ddRmat'*( pos - tvec ); % solid body rotation accelaration
        motacc=dRmat'*( vel - dtvec ); % rotation motion accelaration
        delacc=Rmat'*( acc - ddtvec ); % difference in accelaration
        relacc=( rotacc + motacc + delacc ); % get total relative accelaration 
        
        % transform quadric to new reference
        quad=quadtrans(quad,H); % quadric
        
        % transform data to new reference // this could be done first decompos dRmat and ddRmat if you like
        pos=( Rmat'*( pos - tvec ) ); % position in moving reference
        vel=( Rmat'*( vel ) ); % true velocity in moving reference
        acc=( Rmat'*( acc ) ); % true accelaration in moving reference
        
        %%% write variables beneath
        
        % trajectory frame indexing
        tfrm_index=cat(2,tfrm_index,index);
        
        % perform coordinate transformation quadric
        tfrm_quad=cat(2,tfrm_quad,quad); % to be checked
        
        % perform coordinate transformation position
        tfrm_pos=cat(2,tfrm_pos,pos);
        
        % perform coordinate transformation velocity
        tfrm_vel=cat(2,tfrm_vel,vel);
        
        % velocity difference to moving frame
        tfrm_delvel=cat(2,tfrm_delvel,delvel);
        
        % rotation velocity to moving frame
        tfrm_rotvel=cat(2,tfrm_rotvel,rotvel);
        
        % relative velocity in moving reference
        tfrm_relvel=cat(2,tfrm_relvel,relvel);
        
        % perform coordinate transformation accelaration
        tfrm_acc=cat(2,tfrm_acc,acc);
        
        % perform coordinate transformation accelaration
        tfrm_rotacc=cat(2,tfrm_rotacc,rotacc);
        
        % perform coordinate transformation accelaration
        tfrm_motacc=cat(2,tfrm_motacc,motacc);
        
        % perform coordinate transformation accelaration
        tfrm_delacc=cat(2,tfrm_delacc,delacc);
        
        % relative accelaration
        tfrm_relacc=cat(2,tfrm_relacc,relacc);
        
    end
    
    % create path for processing trajectory reference
    if exist([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame'],'dir')~=7
        mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame'])
    end
    
    % Indexing
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'Index.dat']...
        ,tfrm_index,'w');
    
    % Quadric
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'Quadric.dat']...
        ,tfrm_quad,'w');
    
    % Position
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'Position.dat']...
        ,tfrm_pos,'w');
    
    % Velocity
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'Velocity.dat']...
        ,tfrm_vel,'w');
    
    % VelocityDifference
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'VelocityDifference.dat']...
        ,tfrm_delvel,'w');
    
    % Rotation Velocity
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'RotationVelocity.dat']...
        ,tfrm_rotvel,'w');
    
    % Relative Velocity
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'RelativeVelocity.dat']...
        ,tfrm_relvel,'w');
    
    % Accelaration
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'Accelaration.dat']...
        ,tfrm_acc,'w');
    
    % Accelaration
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'RotationAccelaration.dat']...
        ,tfrm_rotacc,'w');
    
    % Accelaration
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'MotionAccelaration.dat']...
        ,tfrm_motacc,'w');
    
    % Accelaration
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'AccelarationDifference.dat']...
        ,tfrm_delacc,'w');
    
    % Relative Accelaration
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'TrajectFrame' vsl 'RelativeAccelaration.dat']...
        ,tfrm_relacc,'w');
    
    % message
    disp(['post processed traject frame trajectory ',num2str(l)])
    
end

toc

end

