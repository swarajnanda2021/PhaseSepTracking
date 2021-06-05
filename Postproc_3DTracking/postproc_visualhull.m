function postproc_visualhull

% Function to reduce the voxel representation into a 3D matrix with time as
% the third dimension, the rows defined by:
% 1) Centroid x
% 2) Centroid y
% 3) Orientation
% 4) Major-axis length
% 5) Minor-axis length
% 6) z-coordinate

% get globals
global folder date rec prop ctrl post traj

%% create path for processing function storage
if exist([folder date rec vsl 'Postproc3DReco' vsl 'ReducedCavityRepresentation'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DReco' vsl 'ReducedCavityRepresentation'])
end

%%

% Loop over time
for n=post.tproc(1):post.tproc(2)
    disp(['Processing frame' num2str(n) ])
    % Load visual hull and estimate axis
    load([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(n) '.mat'],'Vox_int','X_vox', 'Y_vox', 'Z_vox')
    % Reshape 3D matrix to 2D with transverse dimensions squashed    
    reshp_Vox_int = reshape(Vox_int,size(Vox_int,1)*size(Vox_int,2),size(Vox_int,3)); % Intensity
    reshp_X_Vox = reshape(X_vox,size(Vox_int,1)*size(Vox_int,2),size(Vox_int,3)); % X-coordinate
    reshp_Y_Vox = reshape(Y_vox,size(Vox_int,1)*size(Vox_int,2),size(Vox_int,3)); % Y-coordinate

    % Add X and Y coordinate matrices
    Centroid_cav = [ [sum(reshp_X_Vox.*reshp_Vox_int,1) ; sum(reshp_Y_Vox.*reshp_Vox_int,1)]./sum(reshp_Vox_int,1) ; squeeze(Z_vox(1,1,:))'];

    % Recover the scaling of the units (regionprops treats the data as image and hence in pixels units)
    scaling = (max(X_vox(1,:,1)) - min(X_vox(1,:,1)))/length(X_vox(1,:,1));
    
    
    % Loop over cross-sections and fit ellipses to cross-section
    
    for i=1:size(Z_vox,3)
        
        % Calculate centroid, orientation and major/minor axis length of the ellipse
        s = regionprops(Vox_int(:,:,i),{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
        
        % Convert extracted parameters to pixel coordinates
         
        Centx  = interp1(1:length(X_vox(1,:,1)),X_vox(1,:,1),s.Centroid(1));
        Centy  = interp1(1:length(Y_vox(:,1,1)'),Y_vox(:,1,1)',s.Centroid(2));
        Orient = s.Orientation; % orientation in degreees
        MajAx  = scaling*s.MajorAxisLength;
        MinAx  = scaling*s.MinorAxisLength;
        
        % Append data to matrix
        Visual_hull_reduced(1,i,n) = Centx;
        Visual_hull_reduced(2,i,n) = Centy;
        Visual_hull_reduced(3,i,n) = Orient;
        Visual_hull_reduced(4,i,n) = MajAx;
        Visual_hull_reduced(5,i,n) = MinAx;
        Visual_hull_reduced(6,i,n) = Z_vox(1,1,i);
        
    end
  

end

%% Save data
disp('Saving data')
save([folder date rec vsl 'Postproc3DReco' vsl 'ReducedCavityRepresentation' vsl 'Visual_hull_reduced.mat'],'Visual_hull_reduced'); % save vecgrid

%% Rough
if 0 % plot the ellipse fit to the cavity cross-section in 'pixel' coordinates (not real pixels, just how regionprop works)
    figure(1)

    for i=1:601

        % Calculate centroid, orientation and major/minor axis length of the ellipse
        s = regionprops(Vox_int(:,:,i),{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});

        % Calculate the ellipse line
        theta = linspace(0,2*pi);
        col = (s.MajorAxisLength/2)*cos(theta);
        row = (s.MinorAxisLength/2)*sin(theta);
        M = makehgtform('translate',[s.Centroid, 0],'zrotate',deg2rad(-1*s.Orientation));
        D = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];


    %     surf(X_vox(:,:,i),Y_vox(:,:,i),Vox_int(:,:,i));axis equal; view(2)
        imshow(Vox_int(:,:,i))
        hold on
        plot(D(1,:),D(2,:),'r','LineWidth',2)
        hold off
        drawnow

    end


end




end

