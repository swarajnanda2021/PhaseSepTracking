function postproc_axisymcoord
%postproc_axisymcoord Summary of this function goes here
%   Detailed explanation goes here

global folder date rec prop traj

%% loop over existing folders
trajfolder=dir([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl '*.']); % '\2017*.']);

tic

for f=3:length(trajfolder)
    %% get data
    Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'Index.dat']);
    
%     Quadric=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'Quadric.dat']);
    Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'Position.dat']);
    
    Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'Velocity.dat']);
    RotationVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'RotationVelocity.dat']);
    VelocityDifference=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'VelocityDifference.dat']);
    RelativeVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'RelativeVelocity.dat']);
    
    Accelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'Accelaration.dat']);
    RotationAccelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'RotationAccelaration.dat']);
    MotionAccelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'MotionAccelaration.dat']);
    AccelarationDifference=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'AccelarationDifference.dat']);
    RelativeAccelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Tangentframe' vsl 'RelativeAccelaration.dat']);
    
    %% processing
    
    % indexing of trajectory
    l=str2double(trajfolder(f).name(end));
    L=Index(1,:)==l;
    
    %%% convert index l to nan placeholders %%%
%     Quadric(:,L)=nan;
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
    
    % initiate axi symetric coordinates
    axic_index=zeros(2,0);
    axic_pos=zeros(3,0);
    axic_vel=zeros(3,0);
    axic_delvel=zeros(3,0);
    axic_rotvel=zeros(3,0);
    axic_relvel=zeros(3,0);
    axic_acc=zeros(3,0);
    axic_rotacc=zeros(3,0);
    axic_motacc=zeros(3,0);
    axic_delacc=zeros(3,0);
    axic_relacc=zeros(3,0);
    
    % loop timespan
    for n=unique(Index(2,:))
        
        % timespan
        N=ismember(Index(2,:),n);
        
        % indexing
        index=Index(:,N );
        
        % get quadric
%         quad=Quadric(:,N);
        
        % get position
        pos=Position(:,N );
        
        % get velocity
        vel=Velocity(:,N );
        rotvel=RotationVelocity(:,N );
        delvel=VelocityDifference(:,N );
        relvel=RelativeVelocity(:,N);
        
        % get accelaration
        acc=Accelaration(:,N );
        rotacc=RotationAccelaration(:,N );
        motacc=MotionAccelaration(:,N );
        delacc=AccelarationDifference(:,N );
        relacc=RelativeAccelaration(:,N );
        
        % transform position data to spherical coordinates
        [az,el,rho]=cart2sph(pos(3,:),pos(2,:),pos(1,:));
        pos=[rho
            pi-el
            az];
        
        % rule out nan behaviour
        az(isnan(az))=0;
        el(isnan(el))=0;
        
        % get rotation matrices by basevector transpose
        A=cellfun(@(az,el)azelaxes(az,el),...
            num2cell(reshape(az,1,1,[]),[1 2]),...
            num2cell(reshape(el,1,1,[]),[1 2]),'UniformOutput',false);
        Rmat=cell2mat(cellfun(@(A)[A(:,1) -A(:,3) A(:,2)]',A,'UniformOutput',false));
        
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
        
        %%% write variables beneath
        
        % trajectory frame indexing
        axic_index=cat(2,axic_index,index);
        
        % perform coordinate transformation position
        axic_pos=cat(2,axic_pos,pos);
        
        % perform coordinate transformation velocity
        axic_vel=cat(2,axic_vel,vel);
        
        % velocity difference to moving frame
        axic_delvel=cat(2,axic_delvel,delvel);
        
        % rotation velocity to moving frame
        axic_rotvel=cat(2,axic_rotvel,rotvel);
        
        % relative velocity in moving reference
        axic_relvel=cat(2,axic_relvel,relvel);
        
        % perform coordinate transformation accelaration
        axic_acc=cat(2,axic_acc,acc);
        
        % perform coordinate transformation accelaration
        axic_rotacc=cat(2,axic_rotacc,rotacc);
        
        % perform coordinate transformation accelaration
        axic_motacc=cat(2,axic_motacc,motacc);
        
        % perform coordinate transformation accelaration
        axic_delacc=cat(2,axic_delacc,delacc);
        
        % relative accelaration
        axic_relacc=cat(2,axic_relacc,relacc);
        
    end
    
    % create path for processing trajectory reference
    if exist([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord'],'dir')~=7
        mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord'])
    end
    
    % Indexing
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord' vsl 'Index.dat']...
        ,axic_index,'w');
    
    % Position
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord' vsl 'Position.dat']...
        ,axic_pos,'w');
    
    % Velocity
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord' vsl 'Velocity.dat']...
        ,axic_vel,'w');
    
    % RotationVelocity
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord' vsl 'RotationVelocity.dat']...
        ,axic_rotvel,'w');
    
    % VelocityDifference
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord' vsl 'VelocityDifference.dat']...
        ,axic_delvel,'w');
    
    % RelativeVelocity
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord' vsl 'RelativeVelocity.dat']...
        ,axic_relvel,'w');
    
    % Accelation
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord' vsl 'Accelaration.dat']...
        ,axic_acc,'w');
    
    % MotionAccelaration
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord' vsl 'MotionAccelaration.dat']...
        ,axic_motacc,'w');
    
    % AccelarationDifference
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord' vsl 'AccelarationDifference.dat']...
        ,axic_delacc,'w');
    
    % RelativeAccelaration
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'AxiSymCoord' vsl 'RelativeAccelaration.dat']...
        ,axic_relacc,'w');
    
    % message
    disp(['post processed axi-symetric coordinates trajectory ',num2str(l)])
    
end

toc

end

