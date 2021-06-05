function movie_eulerianreference(dat, grd, slc, vpnt)
%movie_eulerianreference Make movie eulerian reference
%
%   Input
%       dat Particular data to plot
%           - Density
%           - BinnedVelocity
%           - Velocity
%           - Accelaration
%           - DisplacementField
%       grd Type of grid position
%           - GridNode
%           - GridPosition
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
global folder date rec prop mset %plotting

%% load grid
vecgrid=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'vecgrid.mat']);

%% get data
% load general data
Index=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Index.dat']);
Time=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Time.dat']);
GridPosition=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Position.dat']);

% load case sensitive data
switch dat
    case 'Density'
        NumberDensity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'NumberDensity.dat']);
    case 'BinnedVelocity'
        BinnedVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'BinnedVelocity.dat']).*1e-3;
    case 'BinnedGraftieaux'
        BinnedVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'BinnedVelocity.dat']);
    case 'BinnedAccelaration'
        BinnedAccelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'BinnedAccelaration.dat']);
    case 'Polarization'
        Polarization=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Polarization.dat']);
    case 'LinearVelocity'
        LinearVelocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'LinearVelocity.dat']);
    case 'Velocity'
        Velocity=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Velocity.dat']);
    case 'Accelaration'
        Accelaration=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'Accelaration.dat']);
    case 'Displacement'
        error('movie_eulerianreference.m: DisplacementField (curved vectors) is under construction.')
        DisplacementField=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'DisplacementField.dat']);
    case 'FitResidual'
        FitResidual=fload([folder date rec vsl 'Postproc3DTracking' vsl 'EulerianReference' vsl 'FitResidual.dat']);
    otherwise
        error('Movie Eulerian Reference: Case not listed')
end

% Curvature=fload([folder date rec vsl 'Postproc' vsl 'EulerianReference' vsl 'Curvature.dat']);

%% create path for processing function storage
if exist([folder date rec vsl 'Movie3DTracking' vsl 'EulerianReference'],'dir')~=7
    mkdir([folder date rec vsl 'Movie3DTracking' vsl 'EulerianReference'])
end

%% start movie

% name
name=[vpnt dat grd slc];
        
% figure or movie
switch mset.ext
    case {'avi' 'mp4'} % movie
        
        % video object
        if strcmp(mset.ext,'avi')
            vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'EulerianReference' vsl name '.' mset.ext],'Motion JPEG AVI');
        elseif strcmp(mset.ext,'mp4')
            vidObj = VideoWriter([folder date rec vsl 'Movie3DTracking' vsl 'EulerianReference' vsl name '.' mset.ext],'MPEG-4');
        end
        vidObj.Quality = mset.qual;
        vidObj.FrameRate = 10;%prop.fps;
        open(vidObj);
        
    case {'-dbmp' '-depsc'} % figure
        
        % folder
        if exist([folder date rec vsl 'Movie3DTracking' vsl 'EulerianReference' vsl name],'dir')~=7
            mkdir([folder date rec vsl 'Movie3DTracking' vsl 'EulerianReference' vsl name])
        end
        
end

%% Initiate figure
h=figure(1);
h.Color=[1 1 1];
h.Position=[50 50 800 600];
colormap parula

%% render movie of 3d tracks

% loop frames
for n=mset.tspan(1):mset.tspan(2)
    
    % Clear figure
    cla
    
    % Time set
    N=Index(2,:)==n ;
    
    % Grid indexing
    I=ismember(vecgrid.IND,Index(1,N));
    
    % Grid node
    gridnode=vecgrid.X(:,I);
    
    % Get grid positions
    switch grd
        case 'GridNode'
            pos=gridnode;
        case 'GridPosition'
            pos=gridnode+GridPosition(1,N);
    end
    
    % Slice
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
    
    % switch cases
    switch dat
        case 'Density'
            
            % get velocity
            den=NumberDensity(:,N);
            
            % check data
            if nnz(seg)>0
                
                % segment
                den=den(:,seg);
                pos=pos(:,seg);
                
                % color
                col=den;
                
                % plotting
                scatter3(pos(1,:),pos(2,:),pos(3,:),500*ones(size(den)),den,'.')
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [\# / ',mset.units,'^3] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
        case 'BinnedGraftieaux'    
            
            % get velocity
            vel = BinnedVelocity(1:3,N);
            % get unique z-positions
            unq_zpos = unique(pos(3,:));
            % get total number of z-planes
            n_zslices = length(unique(pos(3,:)));
            
            % create a new, larger set of grid coordinate (for the vortex data only)
            % step 1: make meshgrid
            gres = unq_zpos(2)-unq_zpos(1); % calculate grid-resolution from data (grid must be isometric)
            [posx_new,posy_new] = meshgrid(-40:gres:40,-40:gres:40); % 40 chosen for vortex test-case only
            % step 2: unwrap meshgrid to column-vector format (use 3rd
            % dimension to bin data over time and slice
            pos_new(1,:) = reshape(posx_new,[1 (size(posx_new,1)*size(posx_new,2))]);
            pos_new(2,:) = reshape(posy_new,[1 (size(posx_new,1)*size(posx_new,2))]);
            % make accumulation matrix
            vec_accum_mat = nan(size(pos_new,1),size(pos_new,2),n_zslices*length(mset.tspan(1):mset.tspan(2)));
            
            
            % estimate the graftieaux function for each slice and
            % reposition
            for iii=1:n_zslices
                iii
                % step 1: find logical index of all vectors that belong to
                % this z-slice
                slc_idx = pos(3,:)== unq_zpos(iii);
                
                % step 2: calculate the direction-vector matrix
                % (0-diagonal)
                dirx_mat = -minus(pos(1,slc_idx),pos(1,slc_idx)');
                diry_mat = -minus(pos(2,slc_idx),pos(2,slc_idx)');
               
                % step 3: reshape the velocity vector data
                velx_slc = vel(1,slc_idx);
                vely_slc = vel(2,slc_idx);
                
                velx_mat = repmat(velx_slc,length(velx_slc),1);
                vely_mat = repmat(vely_slc,length(vely_slc),1);
                
                vel_iii(:,:,1) = velx_mat;
                vel_iii(:,:,2) = vely_mat;
                vel_iii(:,:,3) = zeros(size(vely_mat));
                
                vel_iii = vel_iii./vecnorm(vel_iii,2,3); % normalize
                
                dir_iii(:,:,1) = dirx_mat;
                dir_iii(:,:,2) = diry_mat;
                dir_iii(:,:,3) = zeros(size(dirx_mat));
                
                dir_iii = dir_iii./vecnorm(dir_iii,2,3); % normalize
                
                % step 5.a.: perform cross_-product
                cross_prod = cross(dir_iii,vel_iii,3);
                
                % step 5.b.: calculate vector norm
                sine_theta = vecnorm(cross_prod,2,3);
                
                % step 6: sum over all each row entries
                fcn_graftieaux = nanmean(sine_theta,2);
                
                [~,indMax(iii)] = max(fcn_graftieaux); % this finds the index in only the slice data
                
                % reposition slice with vortex as centre
                slc_pos(1,:) = pos(1,slc_idx);
                slc_pos(2,:) = pos(2,slc_idx);
                vel_pos(1,:) = vel(1,slc_idx);
                vel_pos(2,:) = vel(2,slc_idx);
                
                slc_pos(1,:) = slc_pos(1,:)-slc_pos(1,indMax(iii));
                slc_pos(2,:) = slc_pos(2,:)-slc_pos(2,indMax(iii));
                
%                 figure(1);hold all;scatter(slc_pos(1,:),slc_pos(2,:))
%                 quiver(slc_pos(1,:),slc_pos(2,:),vel(1,slc_idx),vel(2,slc_idx));
%                 drawnow
                % intersect binning positions with shifted vortex positions
                [~,ia,ib] = intersect(slc_pos',pos_new','rows'); % slc_pos(ia) = pos_new(ib)
                    
%                 figure(1)
%                 scatter(slc_pos(1,:),slc_pos(2,:),'r+')
%                 hold on
%                 scatter(pos_new(1,ib),pos_new(2,ib),'bo')
%                 drawnow
%                 pause
                
%                 figure(1);hold all;plot(ia),plot(ib);drawnow;pause
%                 figure(1);hold all;scatter(slc_pos(1,:),slc_pos(2,:),'r+')
%                 scatter(pos_new(1,:),pos_new(2,:),'bo')
%                 drawnow
%                 pause
%                 clf
%                 size(slc_pos)
%                 size(pos_new)
%                 pause
%                 ia_correct = ia<size(vel_pos,2) ; % remove some spurious stuff (check later?)
%                 ind_rem = find(ia<size(vel_pos,2);
%                 size(ia_correct)
%                 size(ib)
%                 size(ia)
%                 pause
                
%                 plot(ia>size(vel_pos,2))
%                 pause
%                 plot(ia_correct>size(vel_pos,2));drawnow;pause

                vec_accum_mat(1,ib,(n-1)*n_zslices + iii) = vel_pos(1,ia);
                vec_accum_mat(2,ib,(n-1)*n_zslices + iii) = vel_pos(2,ia);
%               

%                 vec_accum_mat(1,ib) = vel_pos(1,ia_correct);
%                 vec_accum_mat(2,ib) = vel_pos(2,ia_correct);
                
                
                % clear vars for next slice
                vel_iii=[];
                dir_iii=[];
                slc_pos=[];
                vel_pos=[];
                ia=[];
                ib=[];
            end
                
            
            size(vec_accum_mat)
            
            vel_spc_avg = nanmean(vec_accum_mat,3);

            close all
            figure(10)
            hold all
            scatter(pos_new(1,:),pos_new(2,:))
            quiver(pos_new(1,:),pos_new(2,:),vel_spc_avg(1,:),vel_spc_avg(2,:),0,'LineWidth',1);axis equal
            drawnow
            pause
            
            
            
            
            
        case 'BinnedVelocity'
            
            % get velocity
            vel=BinnedVelocity(1:3,N);
            
            % check data
            if nnz(seg)>0
                
                % segment
                vel=vel(:,seg);
                pos=pos(:,seg);
                
                % color
                col=sqrt(sum(vel.^2,1));
                
                % plotting
                quiver3c(pos(1,:),pos(2,:),pos(3,:),...
                    vel(1,:),vel(2,:),vel(3,:),0,'LineWidth',1); % displacement by tres

%                 quiver3c(pos(1,:),pos(2,:),pos(3,:),...
%                     vel(1,:),vel(2,:),zeros(size(vel(3,:))),0,'LineWidth',1); % displacement by tres


                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s] }$'],'interpreter','latex','FontSize',14)
%                 caxis([0 (nanmean(col)+6*nanstd(col))])
%                 caxis([-(nanmean(abs(col))+6*nanstd(abs(col))) (nanmean(abs(col))+6*nanstd(abs(col)))])
                caxis([0 5])
%                 caxis([-2 2])
%                 colormap(hot)
                
            end
            
        case 'BinnedAccelaration'
            
            % get velocity
            acc=BinnedAccelaration(1:3,N);
            
            % check data
            if nnz(seg)>0
                
                % segment
                acc=acc(:,seg);
                pos=pos(:,seg);
                
                % color
                col=sqrt(sum(acc.^2,1));
                
                % plotting
                quiver3c(pos(1,:),pos(2,:),pos(3,:),...
                    acc(1,:),acc(2,:),acc(3,:),0,'LineWidth',1); % displacement by tres
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s^2] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
        case 'Polarization'
            
            % get velocity
            pol=Polarization(1:3,N);
            
            % check data
            if nnz(seg)>0
                
                % segment
                pol=pol(:,seg);
                pos=pos(:,seg);
                
                % plotting
                quiver3c(pos(1,:),pos(2,:),pos(3,:),...
                    pol(1,:),pol(2,:),pol(3,:),1.5,'LineWidth',1); % scale unit vector grid to resolution
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [-] }$'],'interpreter','latex','FontSize',14)
                caxis([0 1])
                
            end
            
        case 'LinearVelocity'
            
            % get velocity
            vel=LinearVelocity(1,N);
            
            % check data
            if nnz(seg)>0
                
                % segment
                vel=vel(:,seg);
                pos=pos(:,seg);
                
                % color
                col=sqrt(nansum(vel.^2,1));
                
                % plotting
                scatter3(pos(1,:),pos(2,:),pos(3,:),500*ones(size(vel)),vel,'.')
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
        case 'Velocity'
            
            % get velocity
            vel=Velocity(:,N);
            
            % remove nan for color plotting
            seg=seg & ~isnan(sum(vel,1));
            
            % check data
            if nnz(seg)>0
                
                % segment
                vel=vel(:,seg);
                pos=pos(:,seg);
                
                % color
                col=sqrt(nansum(vel.^2,1));
                
                % plotting
                quiver3c(pos(1,:),pos(2,:),pos(3,:),...
                    vel(1,:),vel(2,:),vel(3,:),0,'LineWidth',1); % displacement by tres
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
        case 'Accelaration'
            
            % get velocity
            acc=Accelaration(:,N);
            
            % remove nan for color plotting
            seg=seg & ~isnan(sum(acc,1));
            
            % check data
            if nnz(seg)>0
                
                % segment
                acc=acc(:,seg);
                pos=pos(:,seg);
                
                % color
                col=sqrt(nansum(acc.^2,1));
                
                % plotting
                quiver3c(pos(1,:),pos(2,:),pos(3,:),...
                    acc(1,:),acc(2,:),acc(3,:),0,'LineWidth',1); % displacement by tres
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s^2] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
        case 'Displacement'
            
            % get velocity
            dispfield=DisplacementField(:,N);
            
            % remove nan for color plotting
            seg=seg & ~isnan(sum(dispfield,1));
            
            % check data
            if nnz(seg)>0
                
                % segment
                dispfield=dispfield(:,seg);
                
                % loop grid cells
                for k=1:size(dispfield,2)
                    
                    % reshape
                    coef=reshape(dispfield(:,k),7,[]);
                    
                    % generate positions over time
                    pdat=[0:floor(mset.tlen/2)
                        repmat(pos(:,k),1,floor(mset.tlen/2)+1)];
                    
                    % evaluate displacement
                    dX=dispeval(coef,pdat,[0 0 0 0]);
                    X=pdat(2:4,:)+dX;
                    
                    % evaluate displacement
                    U=dispeval(coef,pdat,[0 0 0 1]);
                    
                    % define color
                    col=sqrt(sum(U.^2,1));
                    
                    % plot like trajectories
                    hold on
                    scatter3(X(1,:),X(2,:),X(3,:),[],col,'.');
                    hold off
                end
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ \rm{ [',mset.units,' / s] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
        case 'FitResidual'
            
            % get velocity
            res=FitResidual(1,N);
            
            % remove nan for color plotting
            seg=seg & ~isnan(sum(res,1));
            
            % check data
            if nnz(seg)>0
                
                % segment
                res=res(:,seg);
                pos=pos(:,seg);
                
                % color
                col=res;
                
                % plotting
                scatter3(pos(1,:),pos(2,:),pos(3,:),500*ones(size(res)),res,'.')
                
                % colorbar
                c=colorbar;
                ylabel(c,[dat ' $ [\rm{ ',mset.units,'] }$'],'interpreter','latex','FontSize',14)
                caxis([0 (nanmean(col)+6*nanstd(col))])
                
            end
            
    end
    
    
  
    
    %%
    
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
        case 'Isometric'
            view([1 1 1])
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
    title([ vpnt ' ' dat ' ' grd ' ' slc ' Frame ' num2str(n) ],'interpreter','latex')
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
            print(h,[folder date rec vsl 'Movie3DTracking' vsl 'EulerianReference' vsl name vsl 'frame_' num2str(n)],mset.rend,mset.qual,mset.ext)
            
    end
    
    % display message
    disp(['Processed ' vpnt ' ' dat ' ' grd ' ' slc ' Frame ' num2str(n)])
    
end

%% testing only
% % calculate avg. over 3rd dimension
% vel_spc_avg = nanmean(vec_accum_mat,3);
% 
% if strcmp(dat,'BinnedGraftieaux')==1
%     close all
%     figure(10)
%     quiver(pos_new(1,:),pos_new(2,:),vel_spc_avg(1,:),vel_spc_avg(2,:),1,'LineWidth',1);axis equal
%     save('Tracking_data.mat','pos_new','vel_spc_avg','vec_accum_mat')
%     drawnow
% end

%% close video
close(vidObj)

end

