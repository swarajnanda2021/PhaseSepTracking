function postproc_visualfield
%postproc_visualfield Summary of this function goes here
%   Detailed explanation goes here

global folder date rec post traj

%% loop over existing folders
trajfolder=dir([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl '*.']); % '\2017*.']);

tic

for f=3:length(trajfolder)
    %% get data
    Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'Index.dat']);
    
    Quadric=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'Quadric.dat']);
    Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'Position.dat']);
    
    Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'Velocity.dat']);
    RotationVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'RotationVelocity.dat']);
    VelocityDifference=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'VelocityDifference.dat']);
    RelativeVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'RelativeVelocity.dat']);
    
    Accelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'Accelaration.dat']);
    RotationAccelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'RotationAccelaration.dat']);
    MotionAccelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'MotionAccelaration.dat']);
    AccelarationDifference=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'AccelarationDifference.dat']);
    RelativeAccelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(f).name vsl 'TangentFrame' vsl 'RelativeAccelaration.dat']);
    
    %% processing
    
    % indexing of trajectory
    l=str2double(trajfolder(f).name(12:end));%trajfolder(f).name(end));
    L=Index(1,:)==l;
    
    %%% convert index l to nan placeholders %%%
    Quadric(:,L)=nan;
    Position(:,L)=nan;
    
    Velocity(:,L)=nan;
    RotationVelocity(:,L)=nan;
    VelocityDifference(:,L)=nan;
    RelativeVelocity(:,L)=nan;
    
    Accelaration(:,L)=nan;
    RotationAccelaration(:,L)=nan;
    MotionAccelaration(:,L)=nan;
    AccelarationDifference(:,L)=nan;
    RelativeAccelaration(:,L)=nan;
    
    % initiate visual field
    visf_index=zeros(2,0);
    visf_pln=zeros(4,0);
    visf_con=zeros(6,0);
    visf_pos=zeros(3,0);
    visf_vel=zeros(3,0);
    visf_delvel=zeros(3,0);
    visf_rotvel=zeros(3,0);
    visf_relvel=zeros(3,0);
    visf_acc=zeros(3,0);
    visf_rotacc=zeros(3,0);
    visf_motacc=zeros(3,0);
    visf_delacc=zeros(3,0);
    visf_relacc=zeros(3,0);
    
    % loop timespan
    for n=unique(Index(2,:))
        
        % timespan
        N=ismember(Index(2,:),n);
        
        % indexing
        index=Index(:,N );
        
        % get quadric
        quad=Quadric(:,N );
        
        % get position
        pos=Position(:,N );
        
        % get velocity
        vel=Velocity(:,N );
        rotvel=RotationVelocity(:,N );
        delvel=VelocityDifference(:,N );
        relvel=RelativeVelocity(:,N );
        
        % get accelaration
        acc=Accelaration(:,N );
        rotacc=RotationAccelaration(:,N );
        motacc=MotionAccelaration(:,N );
        delacc=AccelarationDifference(:,N );
        relacc=RelativeAccelaration(:,N );
        
        % transform position data to spherical coordinates
        [az,el,rho]=cart2sph(pos(1,:),-pos(2,:),pos(3,:));
        pos=[az
            el
            rho];
        
        % rule out nan behaviour
        az(isnan(az))=0;
        el(isnan(el))=0;
        
        % get rotation matrices by basevector transpose
        A=cellfun(@(az,el)azelaxes(az,el),...
            num2cell(reshape(az,1,1,[]),[1 2]),...
            num2cell(reshape(el,1,1,[]),[1 2]),'UniformOutput',false);
        Rmat=cell2mat(cellfun(@(A)[A(:,2) A(:,3) A(:,1)]',A,'UniformOutput',false));
        
        % velocity decomposition in spherical coordinates
        vel=cell2mat(cellfun(@(R,x)R*x,squeeze(num2cell(Rmat,[1 2]))',...
            num2cell(vel,1),'UniformOutput',false));
        rotvel=cell2mat(cellfun(@(R,x)R*x,squeeze(num2cell(Rmat,[1 2]))',...
            num2cell(rotvel,1),'UniformOutput',false));
        delvel=cell2mat(cellfun(@(R,x)R*x,squeeze(num2cell(Rmat,[1 2]))',...
            num2cell(delvel,1),'UniformOutput',false));
        relvel=cell2mat(cellfun(@(R,x)R*x,squeeze(num2cell(Rmat,[1 2]))',...
            num2cell(relvel,1),'UniformOutput',false));
        
        % accelaration decomposition in spherical coordinates
        acc=cell2mat(cellfun(@(R,x)R*x,squeeze(num2cell(Rmat,[1 2]))',...
            num2cell(acc,1),'UniformOutput',false));
        rotacc=cell2mat(cellfun(@(R,x)R*x,squeeze(num2cell(Rmat,[1 2]))',...
            num2cell(rotacc,1),'UniformOutput',false));
        motacc=cell2mat(cellfun(@(R,x)R*x,squeeze(num2cell(Rmat,[1 2]))',...
            num2cell(motacc,1),'UniformOutput',false));
        delacc=cell2mat(cellfun(@(R,x)R*x,squeeze(num2cell(Rmat,[1 2]))',...
            num2cell(delacc,1),'UniformOutput',false));
        relacc=cell2mat(cellfun(@(R,x)R*x,squeeze(num2cell(Rmat,[1 2]))',...
            num2cell(relacc,1),'UniformOutput',false));
        
        % compute ellipsiod contour
        [pln,con]=qvec2pcon(quad,[0 0 0 1]');
        
        %%% write variables beneath
        
        % trajectory frame indexing
        visf_index=cat(2,visf_index,index);
        
        % perform coordinate transformation position
        visf_pos=cat(2,visf_pos,pos);
        
        % perform coordinate transformation velocity
        visf_vel=cat(2,visf_vel,vel);
        
        % velocity difference to moving frame
        visf_delvel=cat(2,visf_delvel,delvel);
        
        % rotation velocity to moving frame
        visf_rotvel=cat(2,visf_rotvel,rotvel);
        
        % relative velocity in moving reference
        visf_relvel=cat(2,visf_relvel,relvel);
        
        % perform coordinate transformation accelaration
        visf_acc=cat(2,visf_acc,acc);
        
        % perform coordinate transformation accelaration
        visf_rotacc=cat(2,visf_rotacc,rotacc);
        
        % perform coordinate transformation accelaration
        visf_motacc=cat(2,visf_motacc,motacc);
        
        % perform coordinate transformation accelaration
        visf_delacc=cat(2,visf_delacc,delacc);
        
        % relative accelaration
        visf_relacc=cat(2,visf_relacc,relacc);
        
        % plane
        visf_pln=cat(2,visf_pln,pln);
        
        % conic
        visf_con=cat(2,visf_con,con);
        
    end
    
    % create path for processing trajectory reference
    if exist([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField'],'dir')~=7
        mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField'])
    end
    
    % Indexing
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'Index.dat']...
        ,visf_index,'w');
    
    % Position
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'Position.dat']...
        ,visf_pos,'w');
    
    % Velocity
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'Velocity.dat']...
        ,visf_vel,'w');
    
    % RotationVelocity
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'RotationVelocity.dat']...
        ,visf_rotvel,'w');
    
    % VelocityDifference
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'VelocityDifference.dat']...
        ,visf_delvel,'w');
    
    % RelativeVelocity
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'RelativeVelocity.dat']...
        ,visf_relvel,'w');
    
    % Accelation
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'Accelaration.dat']...
        ,visf_acc,'w');
    
    % MotionAccelaration
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'MotionAccelaration.dat']...
        ,visf_motacc,'w');
    
    % AccelarationDifference
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'AccelarationDifference.dat']...
        ,visf_delacc,'w');
    
    % RelativeAccelaration
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'RelativeAccelaration.dat']...
        ,visf_relacc,'w');
    
    % Plane
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'Plane.dat']...
        ,visf_pln,'w');
    
    % Conic
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'TrajectoryFrame' vsl 'Trajectory_' num2str(l) vsl 'VisualField' vsl 'Conic.dat']...
        ,visf_con,'w');
    
    % message
    disp(['post processed visual field trajectory ',num2str(l)])
    
end

toc

end

