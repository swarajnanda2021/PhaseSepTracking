function postproc_ranking
%postproc_ranking Summary of this function goes here
%   Detailed explanation goes here

global folder date rec post traj

%% loop over existing folders
trajfolder=dir([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl '*.']); % '\2017*.']);

tic

for f=3:length(trajfolder)
    %% get data
    Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Visualfield' vsl 'Index.dat']);
    
    Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Visualfield' vsl 'Position.dat']);
    
    Plane=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Visualfield' vsl 'Plane.dat']);
    Conic=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Visualfield' vsl 'Conic.dat']);
    
%     Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Visualfield' vsl 'Velocity.dat']);
%     RotationVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Visualfield' vsl 'RotationVelocity.dat']);
    VelocityDifference=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Visualfield' vsl 'VelocityDifference.dat']);
%     RelativeVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl trajfolder(f).name vsl 'Visualfield' vsl 'RelativeVelocity.dat']);
    
    % initiate axi symetric coordinates
    rank_rad=zeros(1,0);
    rank_cor=zeros(1,0);
    rank_ovlp=zeros(1,0);
    proj_area=zeros(2,0);
    
    % indexing of trajectory
    l=str2double(trajfolder(f).name(end));
    
    % loop timespan
    for n=unique(Index(2,:))
        
        % timespan
        N=ismember(Index(2,:),n);
        
        % get position
        pos=Position(3,N );
        
        % get velocity
        delvel=sqrt(sum(VelocityDifference(:,N ).^2,1));
        
        % get plane
        pln=Plane(:,N);
        
        % get conic
        con=Conic(:,N);
        
        % define radial rank
        [~,rad]=sort(pos); % by sorting radial position
        
        % define corrolation rank (i.e. spatio temporal rank)
        [~,cor]=sort(delvel); 
        if 0
        % compute outline
        [ outl ] = pcon2outl( pln, con );
        
        % convert to spherical coordinates
        [az,el,rho]=cart2sph(outl(1,:,:),-outl(2,:,:),outl(3,:,:));
        
        %%% vectorize beneath %%%
        end
        % visual overlap 
        area=zeros(size(rad)); % projected area
        ovlp=ones(size(rad)); % initiate all first rank
        wght=ones(size(rad)); % initiate all wght one
        if 0
        for i=rad
            if ~isnan(az(:,:,i))
            [~,area(i)]=boundary(az(:,:,i)',el(:,:,i)'); % projected area
            
            rnk=ovlp(i); % select current rank
            
            j=rad(rad~=i); % other obj.
            
            [I,~,K]=find(squeeze(sum(inpolygon(az(:,:,j),el(:,:,j),az(:,:,i),el(:,:,i)),2))); % adjacency
            
            if ~isempty(I)
                ovlp(I)=rnk+1; % raise rank
                wght(I)=wght(I).*(1-K'/size(rho,2)); % poor but working estimate
            end
            end
        end
        
        %%% write metrics beneath
        end
        % trajectory frame indexing
        rank_rad=cat(2,rank_rad,rad);
        rank_cor=cat(2,rank_cor,cor);
        rank_ovlp=cat(2,rank_ovlp,ovlp);
        proj_area=cat(2,proj_area,[area ; wght]);
        
    end
    
    % create path for processing trajectory reference
    if exist([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'Ranking'],'dir')~=7
        mkdir([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'Ranking'])
    end
    
    % Indexing
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'Ranking' vsl 'Index.dat']...
        ,Index,'w');
    
    % Radial Rank
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'Ranking' vsl 'RadialRank.dat']...
        ,rank_rad,'w');
    
    % Corrolation Rank
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'Ranking' vsl 'CorrolationRank.dat']...
        ,rank_cor,'w');
    
    % Rank Visual Overlap
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'Ranking' vsl 'VisualOverlapRank.dat']...
        ,rank_ovlp,'w');
    
    % Projected Area
    fsave([folder date rec vsl ...
        'Postproc3DTracking' vsl 'Trajectoryframe' vsl 'Trajectory_' num2str(l) vsl 'Ranking' vsl 'ProjectedArea.dat']...
        ,proj_area,'w');
    
    % message
    disp(['post processed ranking trajectory ',num2str(l)])
    
end

toc

end

