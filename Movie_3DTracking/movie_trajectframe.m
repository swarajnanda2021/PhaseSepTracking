function movie_trajectframe(name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% get globals
global folder date rec prop mset

%% initiate figure
h=figure(1);
h.Color=[1 1 1];
h.Position=[50 50 800 600];
cmap=colormap('parula');

%% loop over existing folders
trajfolder=dir([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl '*.']); % '\2017*.']);

for i=3:length(trajfolder)
    %% get data
    % indexing
    Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl 'Index.dat']);
    
    % quadric reconstruction
    Quadric=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl 'Quadric.dat']);
    
    % position midpiont
    Position=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl 'Position.dat']);
    
    % Spherical Coordinates
    SphereCoord=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl 'Position.dat']);
    
    % velocity
    Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl 'Velocity.dat']);
%     RotationVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl 'RotationVelocity.dat']);
    VelocityDifference=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl 'VelocityDifference.dat']);
%     RelativeVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl 'RelativeVelocity.dat']);
    
    %% directory
    % create path for processing trajectory reference
    if exist([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame'],'dir')~=7
        mkdir([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame'])
    end
    
    %% start movie
    
    % figure or movie
    switch mset.ext
        case {'avi' 'mp4'} % movie
            
            % video object
            if strcmp(mset.ext,'avi')
                vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl name '.' mset.ext],'Motion JPEG AVI');
            elseif strcmp(mset.ext,'mp4')
                vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl name '.' mset.ext],'MPEG-4');
            end
            vidObj.Quality = mset.qual;
            vidObj.FrameRate = prop.fps;
            open(vidObj);
            
        case {'-dbmp' '-depsc'} % figure
            
            % folder
            if exist([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl name],'dir')~=7
                mkdir([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl name])
            end
            
    end

    %% render movie of 3d tracks
    % mesh sphere
    Xsph=genpntsph(100);
    
    % indexing of trajectory
    l=str2double(trajfolder(i).name(end));
    
    % intiate plotting
    trs_rad=1;
    trs_vel=2;
    
    % loop frames
    for n=unique(Index(2,:))
        % time set
        N=Index(2,:)==n ;
        
        % get indexing
        index=Index(:,N);
        
        % get index focal trajectory
        L=index(1,:)==l;
        
        % get quadrics
        quad=Quadric(:,N);
        
        % get position
        pos=Position(:,N);
        
        % get velocity
        vel=Velocity(:,N);
        delvel=VelocityDifference(:,N);
%         rotvel=RotationVelocity(:,N);
%         relvel=RelativeVelocity(:,N);
        
        % accelaration
        
        % radial position
        rad=SphereCoord(3,N);
        
        % magnitude delvel
        delvel=sqrt(sum(delvel.^2,1));
                
        % segmentation within radius
        seg=rad<trs_rad & delvel<trs_vel;
        
        % plotting
        switch name
            case 'EllipsoidReconstruction'
                % mesh a quadric
                Xshp=qvec2pnts( quad , Xsph);
                
                % repmat radial position
                rad=ceil(64*reshape(repmat(rad,size(Xsph,2),1),1,size(Xsph,2),[])/trs_rad);
                
                % plotting
                Xpnts=reshape(Xshp(:,:,~seg & ~L),3,[]);
                plot3(Xpnts(1,:),Xpnts(2,:),Xpnts(3,:),'k.');
                hold on
                Xpnts=reshape(Xshp(:,:,seg & ~L),3,[]);
                Rad=reshape(rad(:,:,seg & ~L),1,[]);
                scatter3(Xpnts(1,:),Xpnts(2,:),Xpnts(3,:),[],cmap(Rad,:),'.');
                hold off
                
                % local axis
                hold on
                Xpnts=reshape(Xshp(:,:,L),3,[]);
                plot3(Xpnts(1,:),Xpnts(2,:),Xpnts(3,:),'r.');
                quiver3(0,0,0,0.5,0,0,0,'k')
                quiver3(0,0,0,0,0.25,0,0,'k')
                quiver3(0,0,0,0,0,0.25,0,'k')
                text(0.5,0,0,'$\hat{x}$','interpreter','latex')
                text(0,0.25,0,'$\hat{y}$','interpreter','latex')
                text(0,0,0.25,'$\hat{z}$','interpreter','latex')
                hold off
                
                % label title and drawing
                axis tight
                axis equal
                xlim([-1 1])
                ylim([-1 1])
                zlim([-1 1])
                set(gca, 'YDir','reverse')
                view([45 22.5])
                grid minor
                camproj('perspective')
                
                % title and axis
                title(['frame ',num2str(n)],'interpreter','latex','FontSize',14)
                xlabel('$x \;[\rm{m}]$','interpreter','latex','FontSize',14)
                ylabel('$y \;[\rm{m}]$','interpreter','latex','FontSize',14)
                zlabel('$z \;[\rm{m}]$','interpreter','latex','FontSize',14)
                cb=colorbar;
                cb.Ticks=linspace(0,trs_rad,5);
                ylabel(cb,'$d \ [\rm{m}]$','interpreter','latex','FontSize',14)
                
            case 'MidpointTracking'
                % velocity coding
                cvel=ceil(64*delvel/trs_vel);
                
                % plotting
                plot3(pos(1,~seg & ~L),pos(2,~seg & ~L),pos(3,~seg & ~L),'k.');
                hold on
                scatter3(pos(1,seg & ~L),pos(2,seg & ~L),pos(3,seg & ~L),[],...
                    cmap(cvel(seg & ~L),:),'.');
                quiver3(pos(1,seg & ~L),pos(2,seg & ~L),pos(3,seg & ~L)...
                    ,vel(1,seg & ~L),vel(2,seg & ~L),vel(3,seg & ~L),0,'g','LineWidth',1);
                hold off
                
                % local axis
                hold on
                plot3(pos(1,L),pos(2,L),pos(3,L),'r.')
                quiver3(0,0,0,0.5,0,0,0,'k')
                quiver3(0,0,0,0,0.25,0,0,'k')
                quiver3(0,0,0,0,0,0.25,0,'k')
                quiver3(pos(1,L),pos(2,L),pos(3,L),vel(1,L),vel(2,L),vel(3,L),0,'r','LineWidth',1);
                text(0.5,0,0,'$\hat{x}$','interpreter','latex')
                text(0,0.25,0,'$\hat{y}$','interpreter','latex')
                text(0,0,0.25,'$\hat{z}$','interpreter','latex')
                hold off
                
                % label title and drawing
                axis tight
                axis equal
                xlim([-1 1])
                ylim([-1 1])
                zlim([-1 1])
                set(gca, 'YDir','reverse')
                view([45 22.5])
                grid minor
                camproj('perspective')
                
                % title and axis
                title(['frame ',num2str(n)],'interpreter','latex','FontSize',14)
                xlabel('$x \;[\rm{m}]$','interpreter','latex','FontSize',14)
                ylabel('$y \;[\rm{m}]$','interpreter','latex','FontSize',14)
                zlabel('$z \;[\rm{m}]$','interpreter','latex','FontSize',14)
                cb=colorbar;
                cb.Ticks=linspace(0,trs_vel,5);
                ylabel(cb,'$\Delta V \ [\rm{m / s}]$','interpreter','latex','FontSize',14)
                
        end
        
        % drawnow
        drawnow
        
        % write
        switch mset.ext
            case {'avi' 'mp4'} % movie
                
                % getframe for movie
                frm=getframe(h);
                
                % write video
                writeVideo(vidObj, frm);%
                
            case {'-dbmp' '-depsc'} % figure
                
                % print
                print(h,[folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'TrajectFrame' vsl name vsl 'frame_' num2str(n)],mset.rend,mset.qual,mset.ext)
                
        end
        
        % display message
        disp(['processed movie traject frame ' trajfolder(i).name ' ' name ' frame ',num2str(n)])
        
    end % n
    
    %% close video
    close(vidObj)
    
end

end

