function movie_groupmotion(dat, slc, vpnt)
%movie_groupmotion Make movie of the group motion
%
%   Input
%       dat Particular data to plot
%           - CenterGravityBoundingBox
%           - RotationAxisTrajectories
%           - FitVelocityField
%           - GroupVelocity
%           - GroupRotation
%           - GroupStrainRate
%       slc Type of plot
%           - Volume
%           - Slice_[#]
%       vpnt Viewpoint
%           - CameraMotion
%           - TopView
%           - FrontView
%           - SideView
%
%   Output
%       Generate the movie in subfolder
%
%   TODO: Make slc a name value pair, e.g. why is obvious see inline code

%% get globals
global folder date rec prop traj eulr mset

%% load grid
vecgrid=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'vecgrid.mat']);

%% get data

% trajectories
TrajIndex=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Index.dat']);
Time=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Time.dat']);
Pos=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Position.dat']);
Vel=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Trajectories' vsl 'Velocity.dat']);

% segmentation
GroupIndex=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupIndex.dat']); % Group Index
SegIndex=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'SegmentationIndex.dat']); % Group Segmentation

% load case sensitive data
switch dat
    case 'CenterGravityBoundingBox'
        
        % Center of gravity
        GravCen=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'CenterGravity.dat']); % Center of Gravity
        OrienAxis=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'PrincipleAxisOrientation.dat']); % Principle Axis
        SizeAxis=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'PrincipleAxisSize.dat']); % Principle Size
        
        error('movie_groupmotion.m: CenterGravityBoundingBox under construction')
        
    case 'RotationAxisTrajectories'
        
        % Center of rotation
        CenRot=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'CenterRotation.dat']); % Center Rotation
        RotAxis=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'RotationAxis.dat']); % Rotation Axis
        
    case 'FitVelocityField'
        
        % Displacement field
        GravCen=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'CenterGravity.dat']); % Center of Gravity
        DispField=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupDisplacementField.dat']); % Displacement Field Group
        
    case 'GroupVelocity' % Drift
        
        % Velocity drift
        GroupVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupVelocity.dat']); % Displacement Field Group
        
    case 'GroupRotation'
        
        % Drift and rotation
        CenRot=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'CenterRotation.dat']); % Center Rotation
        GroupRotation=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupRotation.dat']); % Displacement Field Group
        
    case 'GroupStrainRate'
        
        % Drift and strain
        CenRot=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'CenterRotation.dat']); % Center Rotation
        GroupStrainRate=fload([folder date rec vsl 'Postproc3DTracking' vsl 'Groups' vsl 'GroupStrainRate.dat']); % Displacement Field Group
        
end

%% create path for processing function storage
if exist([folder date rec vsl 'Movie3DTracking' vsl 'GroupMotion'],'dir')~=7
    mkdir([folder date rec vsl 'Movie3DTracking' vsl 'GroupMotion'])
end

%% start movie
% name
name=[vpnt dat slc];
        
% figure or movie
switch mset.ext
    case {'avi' 'mp4'} % movie
        
        % video object
        if strcmp(mset.ext,'avi')
            vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'GroupMotion' vsl name '.' mset.ext],'Motion JPEG AVI');
        elseif strcmp(mset.ext,'mp4')
            vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'GroupMotion' vsl name '.' mset.ext],'MPEG-4');
        end
        vidObj.Quality = mset.qual;
        vidObj.FrameRate = prop.fps;
        open(vidObj);
        
    case {'-dbmp' '-depsc'} % figure
        
        % folder
        if exist([folder date rec vsl 'Movie3DTracking' vsl 'GroupMotion' vsl name],'dir')~=7
            mkdir([folder date rec vsl 'Movie3DTracking' vsl 'GroupMotion' vsl name])
        end
        
end

%% Initiate figure
h=figure(1);
h.Color=[1 1 1];
h.Position=[50 50 800 600];
colormap parula

%% render movie groupmotion

% loop frames
for n=mset.tspan(1):mset.tspan(2)
    
    % clear figure for new frame
    cla
    
    %%% Select group here G=g which then feeds into S
    
    % Group Set
    G=GroupIndex(2,:)==n;
    
    % Data Segmentation
    S=SegIndex(3,:)==n;
    
    % Trajectory set
    N=TrajIndex(2,:)>=n-floor(traj.tfit/2) & TrajIndex(2,:)<=n+floor(traj.tfit/2) ...
        & ismember(TrajIndex(1,:),SegIndex(2,S));
    
    % Grid node
    gridnode=vecgrid.X; %
    
    % position
    pos=Pos(:,N); % incl scaling dt=1
    
    % Slice
    switch dat
        case 'RotationAxisTrajectories'
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
        otherwise
            switch slc
                case 'Volume'
                    seg=true(1,size(gridnode,2));
                case 'XYsliceZp1'
                    seg=gridnode(3,:)==1;
                case 'XYsliceZ0'
                    seg=gridnode(3,:)==0;
                case 'XYsliceZm1'
                    seg=gridnode(3,:)==-1;
                case 'XYsliceZm2'
                    seg=gridnode(3,:)==-2;
                case 'XYsliceZm3'
                    seg=gridnode(3,:)==-3;
                case 'XYsliceZm4'
                    seg=gridnode(3,:)==-4;
            end
    end
    
    % switch cases
    switch dat
        case 'RotationAxisTrajectories'
            
            % velocity
            vel=Vel(:,N); % incl scaling dt=1
            
            % segment
            pos=pos(:,seg);
            vel=vel(:,seg);
            
            % check data
            if nnz(seg)>0
                
                % color
                col=sqrt(sum(vel.^2,1));
                
                % plot trajectories
                scatter3(pos(1,:),pos(2,:),pos(3,:),[],col,'.')
                
                % get center and axis rotation
                cenrot=CenRot(:,G);
                axrot=3*RotAxis(:,G);
                
                % axis rotation
                hold on
                plot3(cenrot(1,:),cenrot(2,:),cenrot(3,:),'r.','MarkerSize',100)
                quiver3(cenrot(1,:),cenrot(2,:),cenrot(3,:), ...
                    axrot(1,:),axrot(2,:),axrot(3,:),0,'r.','LineWidth',10)
                hold off
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
        case 'FitVelocityField'
            
            % find non-empty grid cells
            for k=find(seg)
                chk=sum((pos-vecgrid.X(:,k)).^2,1)<eulr.sres;
                seg(k)=nnz(chk)>0.5;
            end
            
            % grid data
            xdat=gridnode(:,seg);
            
            % get center and axis rotation
            gravcen=GravCen(:,G);
            
            % get coefficeint data
            coef=reshape(DispField(:,G),7,[]);
            
            % displacement field
            dX=dispeval(coef,[zeros(1,size(xdat,2)) ; xdat-gravcen],[0 0 0 0]);
            
            % Position
            X=xdat+dX;
            
            % displacement field
            U=prop.fps*dispeval(coef,[zeros(1,size(xdat,2)) ; xdat-gravcen],[0 0 0 1]);
            
            % check data
            if nnz(seg)>0
                
                % color
                col=sqrt(sum(U.^2,1));
                
                % plotting
                quiver3c(X(1,:),X(2,:),X(3,:),...
                    U(1,:),U(2,:),U(3,:),0,'LineWidth',1);
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
        case 'GroupVelocity'
            
            % find non-empty grid cells
            for k=find(seg)
                chk=sum((pos-vecgrid.X(:,k)).^2,1)<eulr.sres;
                seg(k)=nnz(chk)>0.5;
            end
            
            % grid data
            xdat=gridnode(:,seg);
            
            % Position
            X=xdat;
            
            % displacement field
            U=repmat(GroupVelocity(:,G),1,nnz(seg));
            
            % check data
            if nnz(seg)>0
                
                % color
                col=sqrt(sum(U.^2,1));
                
                % plotting
                quiver3c(X(1,:),X(2,:),X(3,:),...
                    U(1,:),U(2,:),U(3,:),0,'LineWidth',1);
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
        case 'GroupRotation'
            
            % find non-empty grid cells
            for k=find(seg)
                chk=sum((pos-vecgrid.X(:,k)).^2,1)<eulr.sres;
                seg(k)=nnz(chk)>0.5;
            end
            
            % grid data
            xdat=gridnode(:,seg);
            
            % get center and axis rotation
            cenrot=CenRot(:,G);
            
            % Position
            X=xdat;
            
            % assemble
            W=rotvec2skewmat(GroupRotation(:,G));
            
            % displacement field
            U=(W*(xdat-cenrot));
            
            % check data
            if nnz(seg)>0
                
                % color
                col=sqrt(sum(U.^2,1));
                
                % plotting
                quiver3c(X(1,:),X(2,:),X(3,:),...
                    U(1,:),U(2,:),U(3,:),0,'LineWidth',1);
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
        case 'GroupStrainRate'
            
            % find non-empty grid cells
            for k=find(seg)
                chk=sum((pos-vecgrid.X(:,k)).^2,1)<eulr.sres;
                seg(k)=nnz(chk)>0.5;
            end
            
            % grid data
            xdat=gridnode(:,seg);
            
            % get center and axis rotation
            cenrot=CenRot(:,G);
            
            % Position
            X=xdat;
            
            % assemble
            F=voightvec2symmat(GroupStrainRate(:,G));
            
            % displacement field
            U=(F*(xdat-cenrot));
            
            % check data
            if nnz(seg)>0
                
                % color
                col=sqrt(sum(U.^2,1));
                
                % plotting
                quiver3c(X(1,:),X(2,:),X(3,:),...
                    U(1,:),U(2,:),U(3,:),0,'LineWidth',1);
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
    end
    
    % axis
    view(2)
    axis tight
    axis equal
    
    % Get viewpoint
    switch vpnt
        case 'TopView'
            view([0 90])
        case 'FrontView'
            view([0 0])
        case 'SideView'
            view([90 0])
        case 'CameraMotion'
            view([0 90*(2-erfc( 5*(n-(range(mset.tspan)+1)/2)/((range(mset.tspan)+1)/2)))/2])
    end
    
    % label title and drawing
    xlim(mset.dom(1,:))
    ylim(mset.dom(2,:))
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
            print(h,[folder date rec vsl 'Movie3DTracking' vsl 'GroupMotion' vsl name vsl 'frame_' num2str(n)],mset.rend,mset.qual,mset.ext)
            
    end
    
    % display message
    disp(['Processed ' vpnt ' ' dat ' ' slc ' Frame ' num2str(n)])
    
end

%% close video
close(vidObj)

end

