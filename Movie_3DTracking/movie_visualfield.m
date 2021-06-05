function movie_visualfield(name)
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
    Index=fload([folder date rec  vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField' vsl 'Index.dat']);
    
    % position midpiont
    Position=fload([folder date rec  vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField' vsl 'Position.dat']);
    
    % velocity
    Velocity=fload([folder date rec  vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField' vsl 'Velocity.dat']);
%  !   RotationVelocity=fload([folder date rec '\postproc\trajectoryframe\',trajfolder(i).name,'\tangentframe\RotationVelocity.dat']);
    VelocityDifference=fload([folder date rec  vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField' vsl 'VelocityDifference.dat']);
%  !   RelativeVelocity=fload([folder date rec '\postproc\trajectoryframe\',trajfolder(i).name,'\tangentframe\RelativeVelocity.dat']);
    
    % visual cone
    Plane=fload([folder date rec  vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField' vsl 'Plane.dat']);
    Conic=fload([folder date rec  vsl 'Postproc3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField' vsl 'Conic.dat']);
    
    % radial rank
%  !   RadialRank=fload([folder date rec '\Postproc\TrajectoryFrame\',trajfolder(i).name,'\ranking\RadialRank.dat']);
%   !  CorrolationRank=fload([folder date rec '\Postproc\TrajectoryFrame\',trajfolder(i).name,'\ranking\CorrolationRank.dat']);
    
    %% directory
    % create path for processing trajectory reference
    if exist([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField'],'dir')~=7
        mkdir([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField'])
    end
    
    %% start movie
    
    % figure or movie
    switch mset.ext
        case {'avi' 'mp4'} % movie
            
            % video object
            if strcmp(mset.ext,'avi')
                vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField' vsl name '.' mset.ext],'Motion JPEG AVI');
            elseif strcmp(mset.ext,'mp4')
                vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField' vsl name '.' mset.ext],'MPEG-4');
            end
            vidObj.Quality = mset.qual;
            vidObj.FrameRate = prop.fps;
            open(vidObj);
            
        case {'-dbmp' '-depsc'} % figure
            
            % folder
            if exist([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField' vsl name],'dir')~=7
                mkdir([folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField' vsl name])
            end
            
    end

    %% render movie of 3d tracks
    
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
        
        % get plane conic
        pln=Plane(:,N);
        con=Conic(:,N);
        
        % get position
        pos=Position(:,N);
        
        % get velocity
%         vel=Velocity(:,N);
        delvel=VelocityDifference(:,N);
%         rotvel=RotationVelocity(:,N);
%         relvel=RelativeVelocity(:,N);
        
        % accelaration
        
        % magnitude delvel
        delvel=sqrt(sum(delvel.^2,1));
        
        % segmentation within radius
        seg=pos(3,:)<trs_rad & delvel<trs_vel;
        
        % plotting
        switch name
            case 'ProjectedContours'
                % outline
                Xout = pcon2outl( pln, con );
                
                % compute projections
                [az,el,rho]=cart2sph(Xout(1,:,:),Xout(2,:,:),Xout(3,:,:));
                
                % plot grid
                plot([-pi pi],[0 0],'k--')
                hold on
                plot([pi/2 pi/2],[-pi/2 pi/2],'k--')
                plot([0 0],[-pi/2 pi/2],'k--')
                plot([-pi/2 -pi/2],[-pi/2 pi/2],'k--')
                hold off
                
                % plot
                hold on
                plot(reshape(az(:,:, ~L),1,[]),...
                    reshape(el(:,:, ~L),1,[]),'.','Color',[0.5 0.5 0.5])
                
                scatter(reshape(az(:,:, seg & ~L),1,[]),...
                    reshape(el(:,:, seg & ~L),1,[]),[],...
                    reshape(rho(:,:, seg & ~L),1,[]),'.')
                hold off
                                
                % label title and drawing
                box on
                axis equal
                ax=gca;
                ax.XLim=[-pi pi];
                ax.YLim=[-pi/2 pi/2];
                ax.XTickLabel={'-\pi' '-3\pi/4' '-\pi/2' '-\pi/4' '0' '\pi/4' '\pi/2' '3\pi/4' '\pi'};
                ax.XTick=[-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi];
                ax.YTickLabel={'-\pi/2' '-\pi/4' '0' '\pi/4' '\pi/2'};
                ax.YTick=[-pi/2 -pi/4 0 pi/4 pi/2];
                grid minor
                colorbar
                
                % title and axis
                title(['frame ',num2str(n)],'interpreter','latex')
                xlabel('$\theta \;[\rm{rad}]$','interpreter','latex')
                ylabel('$\varphi \;[\rm{rad}]$','interpreter','latex')
                
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
                print(h,[folder date rec vsl 'Movie3DTracking' vsl 'TrajectoryFrame' vsl trajfolder(i).name vsl 'VisualField' vsl name vsl 'frame_' num2str(n)],mset.rend,mset.qual,mset.ext)
                
        end
        
        % display message
        disp(['processed movie visual field ' trajfolder(i).name ' ' name ' frame ',num2str(n)])
        
    end % n
    
    %% close video
    close(vidObj)
    
end

end

