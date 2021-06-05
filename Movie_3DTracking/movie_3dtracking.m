function movie_3dtracking(dat, slc, vpnt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% get globals
global folder date rec prop eulr mset % plotting

%% get data

% general trajectories
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Time=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Time.dat']);
Pos=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
load([folder date rec vsl 'Postproc3DTracking' vsl 'TrackQualityMetric.mat'],'t_spurious');
Quadric = fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Quadric.dat']);

prop.fps = 10;
% load case sensitive data
switch dat
%     case 'FilteredNoise'
%         
%         % Tracking data
%         Ndata=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Ndata.dat']);
%         
%     case 'ReprojectionError'
%         
%         % Tracking data
%         Tdata=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Tdata.dat']);
%         
%     case 'SweptReprojectionError'
%         
%         Tracking data
%         Tdata=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Tdata.dat']);
        
    case 'Velocity'
        
        % Velocity
        Vel=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
        
    case 'Accelaration'
        
        % Accelaration
        Acc=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accelaration.dat']);
        
%     case 'FluctuationTime'
%         
%         % Accelaration
%         Vel=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);
%         
%         % Accelaration
%         Acc=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Accelaration.dat']);
        
end

%% create path for processing function storage
if exist([folder date rec vsl 'Movie3DTracking' vsl '3DTracking'],'dir')~=7
    mkdir([folder date rec vsl 'Movie3DTracking' vsl '3DTracking'])
end

%% start movie

% name
name=[vpnt dat slc];

% figure or movie
switch mset.ext
    case {'avi' 'mp4'} % movie
        
        % video object
        if strcmp(mset.ext,'avi')
            vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl '3DTracking' vsl name '.' mset.ext],'Uncompressed AVI');
        elseif strcmp(mset.ext,'mp4')
            vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl '3DTracking' vsl name '.' mset.ext],'MPEG-4');
        end
%         vidObj.Quality = mset.qual;
        vidObj.FrameRate = prop.fps;
        open(vidObj);
        
    case {'-dbmp' '-depsc'} % figure
        
        % folder
        if exist([folder date rec vsl 'Movie3DTracking' vsl '3DTracking' vsl name],'dir')~=7
            mkdir([folder date rec vsl 'Movie3DTracking' vsl '3DTracking' vsl name])
        end
        
end

%% Initiate figure
h=figure(1);
h.Color=[1 1 1];
h.Position=[50 50 800 600];%1.5
colormap parula

%% render movie of 3d tracks
for n=200%mset.tspan(1):mset.tspan(2) % loop frames
    
    % time set
    N=Index(2,:)>=n-floor(mset.tlen/2) & Index(2,:)<=n+floor(mset.tlen/2) ;
    
     % Long tracks
    [~,test,~] = unique(Index(1,:),'rows');
    l = unique(Index(1,:));
    L=histcounts(Index(1,:),[l-1/2,max(l)+1/2]);
    l=l(L>=mset.lmin);
    L=ismember(Index(1,:),l);

                
    % remove ghost tracks inside the cavity
    inval = ~ismember(Index(1,:),unique(t_spurious));
    % position
    pos=Pos(:,N&L&inval);
    
    % Slice
    switch slc
        case 'Volume'
            seg=true(1,size(pos,2));
        case 'XYsliceZp1'
            seg=pos(3,:)>=(1-eulr.sres) & pos(3,:)<=(1+eulr.sres);
        case 'XYsliceZ0'
            seg=pos(3,:)>=(0-eulr.sres) & pos(3,:)<=(0+eulr.sres);
        case 'XYsliceZm1'
            seg=pos(3,:)>=(-1-eulr.sres) & pos(3,:)<=(-1+eulr.sres);
        case 'XYsliceZm2'
            seg=pos(3,:)>=(-2-eulr.sres) & pos(3,:)<=(-2+eulr.sres);
        case 'XYsliceZm3'
            seg=pos(3,:)>=(-3-eulr.sres) & pos(3,:)<=(-3+eulr.sres);
        case 'XYsliceZm4'
            seg=pos(3,:)>=(-4-eulr.sres) & pos(3,:)<=(-4+eulr.sres);
    end
    
    % case view
    switch dat
        
%         case 'FilteredNoise'
%             
%             % velocity
%             rpe=Ndata(19,N); % incl scaling dt=1
%             
%             % segment
%             pos=pos(:,seg);
%             rpe=rpe(:,seg);
%             
%             % check data
%             if nnz(seg)>0
%                 
%                 % color
%                 col=rpe;
%                 
%                 % plot trajectories
%                 scatter3(pos(1,:),pos(2,:),pos(3,:),[],col,'.')
%                 
%                 % colorbar
%                 c=colorbar;
%                 ylabel(c,[dat ' $ \rm{ [-] }$'],'interpreter','latex','FontSize',14)
%                 caxis([0 1])
%                 
%             end
%               
%         case 'ReprojectionError'
%             
%             % velocity
%             rpe=Tdata(19,N); % incl scaling dt=1
%             
%             % segment
%             pos=pos(:,seg);
%             rpe=rpe(:,seg);
%             
%             % check data
%             if nnz(seg)>0
%                 
%                 % color
%                 col=rpe;
%                 
%                 % plot trajectories
%                 scatter3(pos(1,:),pos(2,:),pos(3,:),[],col,'.')
%                 
%                 % colorbar
%                 c=colorbar;
%                 ylabel(c,[dat ' $ \rm{ [-] }$'],'interpreter','latex','FontSize',14)
%                 caxis([0 1])
%                 
%             end
            
        case 'Velocity'
            
            % velocity
            vel=Vel(:,N); % incl scaling dt=1
            
            % segment
            pos=pos(:,seg);
            vel=vel(:,seg);
            
            % check data
            if nnz(seg)>0
                
                if strcmp(mset.ellipsoid ,'true')==0
                    % color
                    col=sqrt(sum(vel.^2,1));

                    % plot trajectories
                    scatter3(pos(1,:),pos(2,:),pos(3,:),[],col,'.')

                    % colorbar
                    c=colorbar;
                    ylabel(c,[dat ' $ \rm{ [',mset.units,' / s] }$'],'interpreter','latex','FontSize',14)
                    caxis([0 (nanmean(col)+6*nanstd(col))])
                else
                    col=sqrt(sum(vel.^2,1));
                   
                    [ X_c, X_p, X_ang ] = qvec2eshp( Quadric(:,N&L&inval) );
                    hold on
%                     scatter3(X_c(1,:),X_c(2,:),X_c(3,:),'.') 
                    
                    for i=1:length(X_c(1,:))
                        [Y_el,X_el,Z_el] = ellipsoid(X_c(1,i),X_c(2,i),X_c(3,i),0.5*X_p(1,i),0.5*X_p(2,i),0.5*X_p(3,i),19);
                        % If you knew the angle the axis was rotated by, you could multiply x, y, and z by the rotation matrix.
                        rotm = eula2rotm(X_ang(:,i));
                        points = [X_el(:)-X_c(1,i), Y_el(:)-X_c(2,i), Z_el(:)-X_c(3,i)];
                        points_new = points * rotm + X_c(:,i)';

                        
                        surf(reshape(points_new(:,1),size(X_el)),reshape(points_new(:,2),size(X_el)),...
                            reshape(points_new(:,3),size(X_el)),ones(size((X_el))).*col(i),'edgecolor','none');
                        drawnow
                        
                    end
                    hold off
                    c=colorbar;
                    ylabel(c,[dat ' $ \rm{ [',mset.units,' / s] }$'],'interpreter','latex','FontSize',14)
                    caxis([0 (nanmean(col)+6*nanstd(col))])

                end
            end
            
        case 'Accelaration'
            
            % velocity
            acc=Acc(:,N); % incl scaling dt=1
            
            % segment
            pos=pos(:,seg);
            acc=acc(:,seg);
            
            % check data
            if nnz(seg)>0
                
                % color
                col=sqrt(sum(acc.^2,1));
                
                % plot trajectories
                scatter3(pos(1,:),pos(2,:),pos(3,:),[],col,'.')
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s^2] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
%         case 'FluctuationTime'
%             
%             % velocity
%             fluctime=sqrt(sum(Vel(:,N).^2,1))...
%                 ./sqrt(sum(Acc(:,N).^2,1))/prop.fps; % incl scaling dt=1
%             
%             % segment
%             pos=pos(:,seg);
%             fluctime=fluctime(:,seg);
%             
%             % check data
%             if nnz(seg)>0
%                 
%                 % color
%                 col=fluctime;
%                 
%                 % plot trajectories
%                 scatter3(pos(1,:),pos(2,:),pos(3,:),[],col,'.')
%                 
%                 % colorbar
%                 c=colorbar;
%                 ylabel(c,[dat ' $ \rm{ [s] }$'],'interpreter','latex','FontSize',14)
%                 caxis([0 2])
%                 
%             end
            
    end
    
    % axis
    view(2)
    axis tight
    axis equal
    
    caxis([4000 6000])
    
    
    
    load([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(n) '.mat'],'Vox_int','X_vox', 'Y_vox', 'Z_vox','Centroid_cav')
    [fo,vo] = isosurface(Y_vox,X_vox,Z_vox,permute(Vox_int,[1 2 3]));
    [fe,ve,ce] = isocaps(Y_vox,X_vox,Z_vox,permute(Vox_int,[1 2 3]));
    
    p1 = patch('Faces', fo, 'Vertices', vo);       % draw the outside of the volume
    p1.FaceColor = 'blue';
    p1.EdgeColor = 'none';

    p2 = patch('Faces', fe, 'Vertices', ve, ...    % draw the end caps of the volume
       'FaceVertexCData', ce);
    p2.FaceColor = 'interp';
    p2.EdgeColor = 'none';

    camlight(40,40)                                % create two lights 
    camlight(-20,-20)
%     lighting gouraud
    set(gca,'Color','k')
    
    
    % Get viewpoint
    switch vpnt
        case 'Isometric'
            view([0.1 0.1 1])
        case 'TopView'
            view([0 90])
        case 'FrontView'
            view([0 0])
        case 'SideView'
            view([90 0])
        case 'CameraMotion'
            view([0 90*(2-erfc( 5*(n-(range(mset.tspan)+1)/2)/((range(mset.tspan)+1)/2)))/2])
    end
    h.WindowState = 'maximized';
    % label title and drawing
%     xlim(mset.dom(1,:))
%     ylim(mset.dom(2,:))

    xlim([mean(Centroid_cav(2,:))-15 mean(Centroid_cav(2,:))+15])
    ylim([mean(Centroid_cav(1,:))-15 mean(Centroid_cav(1,:))+15])

    zlim(mset.dom(3,:))
    grid minor
    box on
    camproj('perspective')
    
    % title
    title([ vpnt ' ' dat ' ' slc ' Frame ' num2str(n) ],'interpreter','latex')
    xlabel(['$X \;[\rm{' mset.units '}]$'],'interpreter','latex','FontSize',14)
    ylabel(['$Y \;[\rm{' mset.units '}]$'],'interpreter','latex','FontSize',14)
    zlabel(['$Z \;[\rm{' mset.units '}]$'],'interpreter','latex','FontSize',14)
    
    % print
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
            print(h,[folder date rec vsl 'Movie3DTracking' vsl '3DTracking' vsl name vsl 'frame_' num2str(n)],mset.rend,mset.qual,mset.ext)
            
    end
    
    % display message
    disp(['Processed ' vpnt ' ' dat ' ' slc ' Frame ' num2str(n)])
    
    pause
    clf
    
    
end % n

%% close video
close(vidObj)

end

