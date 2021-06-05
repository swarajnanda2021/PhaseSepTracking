function plot_tracking(rmax,rmin)

%% Code that plots tracks alongside the cavity for a chosen max-radius, r2 (inner), r1(outer) 
% (mm) and trailing which sets the number of timesteps the tracks will
% include

%% get globals
global folder date rec prop mset post %plotting


%% Load tracking data
post=importdata([folder date rec vsl 'Postproc3DTracking' vsl 'post.mat']);
load([folder date rec vsl 'Vortex3DTracking' vsl 'VortexTrajectories' vsl 'Pobj_2_Lmin_15.mat'],'Pobj_2'); % save vecgrid

% Pobj_2 is a matrix with the following entries
% 1) Track-index
% 2) Time-step
% 3-5) position x-z
% 6-7) Centroid of ellipse next to cavity
% 8) Orientation of cross-sectional ellipse
% 9-10) Size of the major and minor axis 
% 11) Distance (euclidean) from centroid of fitted ellipse


% Segment tracks based on radius
if 0
    kk = unique(Pobj_2(1,Pobj_2(end,:) < rmax & Pobj_2(end,:) > rmin));   % Find index
    test = ismember(Pobj_2(1,:),kk);            % Intersect the unique indices with Pobj_2 to find all particles belonging to the unique indices
    segmented_pobj2 = Pobj_2(:,test==1);        % Segment it out of Pobj_2
end

if 1
    
   segmented_pobj2 =Pobj_2(:,Pobj_2(end,:) < rmax & Pobj_2(end,:) > rmin); 
    
end


% Assign a video output variable for the drawn frames
vidfile = VideoWriter(['animation_tracking_shape_rmax' num2str(rmax) '_rmin_' num2str(rmin) ],'MPEG-4');
vidfile.FrameRate = 10; % always write in 10 fps
open(vidfile);



for k=post.tproc(1) : post.tproc(2)
    
    % Load visual hull data for this timestep
    load([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(k) '.mat'],'Vox_int','X_vox', 'Y_vox', 'Z_vox')
    
    
    
    % Begin plotting
    %% First plot the visual hull
    [fo,vo] = isosurface(Y_vox,X_vox,Z_vox,Vox_int);
    [fe,ve,ce] = isocaps(Y_vox,X_vox,Z_vox,Vox_int);
    hfig = figure(1)
    p1 = patch('Faces', fo, 'Vertices', vo);       % draw the outside of the volume
    p1.FaceColor = 'blue';
    p1.EdgeColor = 'none';

    p2 = patch('Faces', fe, 'Vertices', ve, ...    % draw the end caps of the volume
       'FaceVertexCData', ce);
    p2.FaceColor = 'interp';
    p2.EdgeColor = 'none';

    camlight(40,40)                                % create two lights 
    camlight(-20,-20)
    lighting gouraud
    zlim([-10 50])
    axis equal
    axis tight

    hfig.WindowState = 'maximized';
    
    %% Next plot the tracks
    % Segment track indices that are present at the kth timestep
    ia = segmented_pobj2(2,:)==k;
    hold on
    scatter3(segmented_pobj2(3,ia),segmented_pobj2(4,ia),segmented_pobj2(5,ia),[],segmented_pobj2(2,ia),'o','filled')
    hold off
    axis equal
    axis tight
    colorbar
    view(2)
    xlim([-15 15])
    ylim([-15 15])
    
    drawnow
%     pause
    Frm = getframe(gcf);
    clf
    writeVideo(vidfile, Frm);
    
    
%     track_idxs = segmented_pobj2(1,ia==1);
%     ib = ismember(segmented_pobj2(1,:),track_idxs);
%     track_data = segmented_pobj2(:,ib==1);
%    
%     
%     % Find unique track indices
%     [~,ia,~] = unique(track_data(1,:));
%     L = [ia(2:end) - ia(1:end-1) ; length(track_data(1,:)) - ia(end)]'; % Track-length in particles, not timesteps
%     hold on
%     for kk=1:length(ia) %loop over all tracks and plot only (k-trailing_length)th to kth timestep
%         % Take out the track from track_data
%         mini_data = track_data(:,ia(kk):ia(kk)+L(kk)-1); 
%         
%         % segment out part that is relevant to the plotting
%         mini_data_seg = mini_data(:,mini_data(2,:)<=k & mini_data(2,:) >=(k-trailing_length));
%          
% %         % Perform track fitting using third order polynomial for all tracks
% %         % that have a constituent present in the kth timestep
% %         [ coef ,~] = polytraj( mini_data(1:5,:), 3 , 0 ); % estimate residual over all 
% %         % Calculate polynomial coefficient over subpoints
% %         [ tdat ] = polyeval( coef, mini_data(1:5,:), 0 );
% %         % Calculate the euclidean residual between the fitted track and the
% %         % position
% %         plot3(mini_data_seg(3,:),mini_data_seg(4,:),mini_data_seg(5,:),'k')
%         scatter3(mini_data_seg(3,:),mini_data_seg(4,:),mini_data_seg(5,:),[],mini_data_seg(2,:),'.')
%         
%         
%     end
%     
    
%     figure()
    
    
%     pause
    
    
    
    
end

close(vidfile)

end

