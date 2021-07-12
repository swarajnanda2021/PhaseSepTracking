

% This code reconstructs the visual hull of the cavity based on silhouette
% data available in a separate folder in voxel space

close all
clear all
clc

global folder date rec cal ctrl prop plotting ...
    Kmat Rmat tvec Pmat Dpar Dmap Pmap camprop

%% Additional Paths
direc= cd('..\'); % one folder up
addpath(genpath([cd '\Function_Libraries']))
cd(direc) %  back to original folder
clear direc % remove variable

%% Locations
folder='C:\Work\PostGraduation\Scripts'; % Desktop Location
date='\ParticleTracking'; % Date Folder
rec='\Measurement38'; % Recording
cal='\Calibration_200514_141513'; % Calibration files

%% create path for processing function storage
if exist([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull'])
end

if exist([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Segmented_VH'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Segmented_VH'])
end
%% Load/Set variables and files

% Load camera calibration data
[Dmap,Pmap,Dpar,Kmat,Rmat,tvec,camprop]=load_Cal;

%% Controls
plotting='off'; % switch of all figures for cluster! ( 'off' | 'on')
if 1
    [prop,ctrl] = proc_set;
else
    load([folder date rec vsl 'Preproc3DReco' vsl 'prop.mat'])
    load([folder date rec vsl 'Preproc3DReco' vsl 'ctrl.mat'])
end
%%

% 1. Create a voxel space

% 1-D space
x_space = prop.range_x(1):1:prop.range_x(2);
y_space = prop.range_y(1):1:prop.range_y(2);
z_space = prop.range_z(1):1:prop.range_z(2);
% Voxel space meshgrid
[X_vox, Y_vox, Z_vox] = meshgrid(x_space,y_space,z_space);

% 2. Reshape to homogeneous 3D coordinate system
Vox_points = [reshape(X_vox,[1 size(X_vox,1)*size(X_vox,2)*size(X_vox,3)])   ;     reshape(Y_vox,[1 size(Y_vox,1)*size(Y_vox,2)*size(Y_vox,3)])     ;    reshape(Z_vox,[1 size(Z_vox,1)*size(Z_vox,2)*size(Z_vox,3)])    ;     ones(size(reshape(X_vox,[1 size(X_vox,1)*size(X_vox,2)*size(X_vox,3)])))  ];
Vox_int    = ones(1,size(Vox_points,2));


%% Begin loop


% Load deep learning model
load([folder vsl date vsl 'Deep_Learning_Model_Database/' prop.modelname '.mat'])

%%
% Pre-allocate for storing segmented images

% Segmented_img_stack = zeros(size(import_frames({folder date rec},prop.ext,1,1),1),size(import_frames({folder date rec},prop.ext,1,1),2),length(ctrl.lfoc),ctrl.tproc(2)+1-ctrl.tproc(1)    ); % Very expensive, in hindsight
% Segmented_img_stack = zeros(size(import_frames({folder date rec},prop.ext,1,1),1),size(import_frames({folder date rec},prop.ext,1,1),2),length(ctrl.lfoc),1   ); % fixing it
% Occlusion_ratio = zeros(length(ctrl.lfoc),ctrl.tproc(2)+1-ctrl.tproc(1));   

for tstep=ctrl.tproc(1):ctrl.tproc(2)
    tstep
    Vox_int    = ones(1,size(Vox_points,2));
    for iter=1:2
        iter
        for frame=ctrl.lfoc
            frame
            tic

            % Load image 
            disp('Loading distortion param and image')
            % Load distortion calibration
            distortn = load([folder date cal vsl 'wrp_',num2str(frame),'.mat'],'wrp');

            toc
            if 0 % When using silhouette sketches 
                Im1 = double(rgb2gray(imread([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'SilhouetteSketches' vsl 'Frm_' num2str(tstep) '_Cam_' num2str(frame) '.png' ])));
                % Binarize silhouette
                Im1 = double(imcomplement(imbinarize(Im1,200)));
            end

            tic
            disp('Performing Inference')
            if 1 % Using the CNN model

                % Load reference image for normalization (consider replacing this?)
                Img_stack_ref = imread([folder vsl date vsl '\Measurement37' '\Preproc3DReco\ImgSeg\Training\trainingImages\image_1.tiff']);
                % Load image of this timestep and frame 
                Img_stack = import_frames({folder date rec},prop.ext,tstep,frame);
                % Preprocess the images to get it to the range you trained the
                % model at

                if frame==1
                    Img_stack(Img_stack>250) = 0;
                    Img_stack(Img_stack<50) = 0;

                elseif frame==2
                    Img_stack(Img_stack>80)=0;
                    Img_stack(Img_stack<10)=0;

                elseif frame==3
                    Img_stack(Img_stack>100)=0;
                    Img_stack(Img_stack<20)=0;

                elseif frame==4
                    Img_stack(Img_stack>250)=0;
                    Img_stack(Img_stack<40)=0;

                elseif frame==5
                    Img_stack(Img_stack>100)=0;
                    Img_stack(Img_stack<20)=0;

                end

                % Construct block processor function
                fun = @(img_stac) labeloverlay(img_stac.data.*mean2(Img_stack_ref)/mean2(img_stac.data),semanticseg((abs((img_stac.data).*mean2(Img_stack_ref)/mean2(img_stac.data))), net));
                % Perform block inference over 256*n:256*n part of the image
                % (with hopefully good intersection over all frames)
                I2 = blockproc(Img_stack(1:1024,1:1024),[256 256],fun,'TrimBorder', false);

                % Shift and infer again so that two predictions can be
                % multiplied
                shifter = 10;
                I2_shifted = blockproc(Img_stack(1+shifter:1024+shifter,1+shifter:1024+shifter),[256 256],fun,'TrimBorder', false);


                % Binarize and resize the image
                % Part 1
                Segmented_img = zeros(size(Img_stack));
                Segmented_img(1:1024,1:1024) = rgb2gray(I2);
                Segmented_img(Segmented_img>15) = 0;
                Segmented_img(Segmented_img==15) = 1;
                % Part 2
                Segmented_img2 = zeros(size(Img_stack));
                Segmented_img2(1+shifter:1024+shifter,1+shifter:1024+shifter) = rgb2gray(I2_shifted);
                Segmented_img2(Segmented_img2>15) = 0;
                Segmented_img2(Segmented_img2==15) = 1;

                % Store segmented image in Im1
                Im1 = Segmented_img2;%.*Segmented_img;
    %             imshowpair(edge(Im1),Img_stack)
                % Save Segmented image into a 4-D matrix
    %             Segmented_img_stack(:,:,frame,1) = Im1;

                % Calculate occlusion metric
    %             Img_segmented = Im1(max(1,ctrl.croi(frame,1)):min(1024,ctrl.croi(frame,3)),max(1,ctrl.croi(frame,2)):min(1024,ctrl.croi(frame,4)));

    %             Occlusion_ratio(frame,tstep) = length(find(Img_segmented==1))/(size(Img_segmented,1)*size(Img_segmented,2)); %ratio of occluded area to sensor area




                if 0 
                    figure(69)
                    hold all
                    imshowpair(edge(Im1),Img_stack,'ColorChannel','red-cyan')

                    axis equal
                    axis tight
                    drawnow
                    pause
                end


            end
            toc


            tic
            disp('Inverting the distortion function')
            % make the input arrays for Dmap
            imgsize = size(distortn.wrp.x);

            pix_x = 1  : imgsize(1);
            pix_y = 1  : imgsize(2);
            [mgrid_pixx, mgrid_pixy] = meshgrid(pix_x,pix_y);

            input_array_pixx = reshape(mgrid_pixx,1,length(pix_x)*length(pix_y));
            input_array_pixy = reshape(mgrid_pixy,1,length(pix_x)*length(pix_y));

            input_array_pixxy = [input_array_pixx; input_array_pixy];

            % calculate the output arrays for Dmap
            output_array_pixxy  = Dmap(Dpar{1,frame}(1),Dpar{1,frame}(2),Dpar{1,frame}(3),Dpar{1,frame}(4),Dpar{1,frame}(5),Dpar{1,frame}(6),Dpar{1,frame}(7),Dpar{1,frame}(8),Dpar{1,frame}(9),Dpar{1,frame}(10),Dpar{1,frame}(11),Dpar{1,frame}(12),Dpar{1,frame}(13),Dpar{1,frame}(14),Dpar{1,frame}(15),Dpar{1,frame}(16),Dpar{1,frame}(17),Dpar{1,frame}(18),Dpar{1,frame}(19),Dpar{1,frame}(20),input_array_pixx,input_array_pixy);
            % 
            tform = images.geotrans.PolynomialTransformation2D(input_array_pixxy',output_array_pixxy',3);

            toc



            tic
            % Construct projection matrix
            normPMat = [Rmat{frame} tvec{frame}];

            % Project voxels
            disp('Projecting voxels')
            vox_proj = (normPMat*Vox_points);
            % Normalize
            vox_proj = vox_proj./vox_proj(3,:);
            % Bring to camera frame
            vox_proj_dew = Kmat{frame}*vox_proj;
            % Normalize to get dewarped pixel representation
            vox_proj_dew = vox_proj_dew./vox_proj_dew(3,:);
            % Distort to pixel frame
            vox_proj_war = transformPointsInverse(tform,vox_proj_dew(1:2,:)');
            % Find index of points outside of the camera frame
            ind_outside = vox_proj_war(:,1)>imgsize(1) | vox_proj_war(:,2)>imgsize(2);

            toc

    %         % Plot check
    %         figure(2); hold all; imagesc(Im1);colormap(gray);scatter(vox_proj_war(:,2),vox_proj_war(:,1),'b.');scatter(vox_proj_war(ind_outside,2),vox_proj_war(ind_outside,1),'r.')
    %         drawnow
    %         pause

            tic
            % Query intensity
            disp('Querying intensity')
            % Convert subscript to indices 
            vox_proj_cam_ind = sub2ind(size(Im1),ceil(vox_proj_war(:,1)),ceil(vox_proj_war(:,2)));
            vox_proj_int = Im1(vox_proj_cam_ind);

            % Query index of voxel that is value 1
            Vox_int = Vox_int.*(vox_proj_int'==1);

            toc
        end



        Vox_int = reshape(Vox_int,size(X_vox,1),size(X_vox,2),size(X_vox,3));



        tic 
        disp('Estimating axis')
        % Axis estimate by averaging the transverse coordinates of the shape

        % Reshape 3D matrix to 2D with transverse dimensions squashed    
        reshp_Vox_int = reshape(Vox_int,size(Vox_int,1)*size(Vox_int,2),size(Vox_int,3)); % Intensity
        reshp_X_Vox = reshape(X_vox,size(Vox_int,1)*size(Vox_int,2),size(Vox_int,3)); % X-coordinate
        reshp_Y_Vox = reshape(Y_vox,size(Vox_int,1)*size(Vox_int,2),size(Vox_int,3)); % Y-coordinate

        % Add X and Y coordinate matrices
        Centroid_cav = [ [sum(reshp_X_Vox.*reshp_Vox_int,1) ; sum(reshp_Y_Vox.*reshp_Vox_int,1)]./sum(reshp_Vox_int,1) ; squeeze(Z_vox(1,1,:))'];

        % Estimate surface distance from 
        surf_dist  = ((X_vox- reshape(repmat(Centroid_cav(1,:),size(X_vox,1),1,size(X_vox,1)),[size(X_vox,1) ,size(X_vox,2), size(X_vox,3)]   )).^2 +...
                     (Y_vox- reshape(repmat(Centroid_cav(2,:),size(Y_vox,1),1,size(Y_vox,1)),[size(Y_vox,1) ,size(Y_vox,2), size(Y_vox,3)]   )).^2).^0.5;

        if 0 %plotcheck
            plot3(Centroid_cav(1,:),Centroid_cav(2,:),squeeze(Z_vox(1,1,:)),'r*'); axis equal
        end

        toc

        if iter==2
            
            
            %% Save visual hull data
            % If you want to save the visual hull, toggle this. But it is
            % limited by hard drive memory. 
            if 0

                tic
                disp('Saving VH data')
                save([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(tstep) '.mat'],'Vox_int','X_vox', 'Y_vox', 'Z_vox','Centroid_cav','surf_dist')
                toc 
                
            end
            
            %% Post-process it on the fly and save reduced visual-hull data
            
            if tstep==ctrl.tproc(1)
                Visual_hull_reduced = zeros(8,size(Z_vox,3),0);
            end
            % Recover the scaling of the units (regionprops treats the data as image and hence in pixels units)
            scaling = (max(X_vox(1,:,1)) - min(X_vox(1,:,1)))/length(X_vox(1,:,1));


            % Loop over cross-sections and fit ellipses to cross-section
            Centx = squeeze(sum(Y_vox.*Vox_int,[1 2]))./length(Vox_int(find(Vox_int==1))) ; % units object space
            Centy = squeeze(sum(X_vox.*Vox_int,[1 2]))./length(Vox_int(find(Vox_int==1))) ;

            for i=1:size(Z_vox,3) 


                % Handle the necessary transformations
                edges = edge(squeeze(permute(Vox_int(:,:,i),[2 1 3])   ));
        %         edges = edge(squeeze(permute(Vox_int(:,:,i),[1 2 3])   ));
                [ind2,ind1] = find(edges);
                if length(find(edges)) > 30 % some threshold for minimum number of points needed for the ellipse fit
                    % find x and y coordinates
                    edgpts = [squeeze(X_vox(1,ind2,i));squeeze(Y_vox(ind1,1,i))'];

                    % REQUIRES conic and quadric fitting toolbox
                    fobj=ellipticalFit(edgpts); %Perform the fit

                    if 0
                        figure(1)
                        hold on
                        imagesc(squeeze(X_vox(1,:,i)),squeeze(Y_vox(:,1,i)),edge(squeeze(Vox_int(:,:,i))));
                        scatter(edgpts(1,:),edgpts(2,:))
                        [hFit,hData]=plot(fobj,{'LineWidth',2},{'MarkerFaceColor','c','MarkerEdgeColor','k'}); %Visualize the fit
                        drawnow
                        pause
                    end

                    MajAx = fobj.a;
                    MinAx = fobj.b;
                    Orient= fobj.angle; % in degrees
                    Area  = sum(sum(squeeze(Vox_int(:,:,i)))).*scaling^2; % units^2 (in this work, mm^2)


                    % Append data to matrix
                    Visual_hull_reduced(1,i,tstep) = tstep;
                    Visual_hull_reduced(2,i,tstep) = Centx(i);
                    Visual_hull_reduced(3,i,tstep) = Centy(i);
                    Visual_hull_reduced(4,i,tstep) = Orient;
                    Visual_hull_reduced(5,i,tstep) = MajAx;
                    Visual_hull_reduced(6,i,tstep) = MinAx;
                    Visual_hull_reduced(7,i,tstep) = Area;
                    Visual_hull_reduced(8,i,tstep) = Z_vox(1,1,i);
                else

                    Visual_hull_reduced(1,i,tstep) = tstep;
                    Visual_hull_reduced(2,i,tstep) = nan;
                    Visual_hull_reduced(3,i,tstep) = nan;
                    Visual_hull_reduced(4,i,tstep) = nan;
                    Visual_hull_reduced(5,i,tstep) = nan;
                    Visual_hull_reduced(6,i,tstep) = nan;
                    Visual_hull_reduced(7,i,tstep) = nan;
                    Visual_hull_reduced(8,i,tstep) = Z_vox(1,1,i);






                end

            end
            
            
            
            
            
            
        end



        if 0
            %     Plotting
    %             load([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(tstep) '.mat'],'Vox_int','X_vox', 'Y_vox', 'Z_vox')
                [fo,vo] = isosurface(X_vox,Y_vox,Z_vox,Vox_int);
                [fe,ve,ce] = isocaps(X_vox,Y_vox,Z_vox,Vox_int);
                figure(1)
                hold all
                p1 = patch('Faces', fo, 'Vertices', vo);       % draw the outside of the volume
                isocolors(X_vox,Y_vox,Z_vox,surf_dist,p1)
                p1.FaceColor = 'interp';
                p1.EdgeColor = 'none';

                p2 = patch('Faces', fe, 'Vertices', ve, ...    % draw the end caps of the volume
                   'FaceVertexCData', ce);
                p2.FaceColor = 'interp';
                p2.EdgeColor = 'none';

    %             camlight(40,40)                                % create two lights 
    %             camlight(-20,-10)
    %             lighting gouraud
                axis equal
                colormap gray

                plot3(Centroid_cav(1,:),Centroid_cav(2,:),squeeze(Z_vox(1,1,:)),'r*'); axis equal

                pause
                clf
        end
        
        
        % Calculate cavity location and update voxel space
        cavity_x = nanmean(Centroid_cav(1,:));
        cavity_y = nanmean(Centroid_cav(2,:));
        
        
        % Update voxel space to fine
        
        
        % 1. Initialize Voxel space

        % 1-D space
        x_space = (cavity_x-4):prop.del_vox:(cavity_x+4);
        y_space = (cavity_y-4):prop.del_vox:(cavity_y+4);
        z_space = prop.range_z(1):prop.del_vox:prop.range_z(2);
        % Voxel space meshgrid
        [X_vox, Y_vox, Z_vox] = meshgrid(x_space,y_space,z_space);

        % 2. Reshape to homogeneous 3D coordinate system
        Vox_points = [reshape(X_vox,[1 size(X_vox,1)*size(X_vox,2)*size(X_vox,3)])   ;     reshape(Y_vox,[1 size(Y_vox,1)*size(Y_vox,2)*size(Y_vox,3)])     ;    reshape(Z_vox,[1 size(Z_vox,1)*size(Z_vox,2)*size(Z_vox,3)])    ;     ones(size(reshape(X_vox,[1 size(X_vox,1)*size(X_vox,2)*size(X_vox,3)])))  ];
        Vox_int    = ones(1,size(Vox_points,2));


        
    end
    
    % Back to coarse
    % 1-D space
    x_space = prop.range_x(1):1:prop.range_x(2);
    y_space = prop.range_y(1):1:prop.range_y(2);
    z_space = prop.range_z(1):1:prop.range_z(2);
    % Voxel space meshgrid
    [X_vox, Y_vox, Z_vox] = meshgrid(x_space,y_space,z_space);
    
    % 2. Reshape to homogeneous 3D coordinate system
    Vox_points = [reshape(X_vox,[1 size(X_vox,1)*size(X_vox,2)*size(X_vox,3)])   ;     reshape(Y_vox,[1 size(Y_vox,1)*size(Y_vox,2)*size(Y_vox,3)])     ;    reshape(Z_vox,[1 size(Z_vox,1)*size(Z_vox,2)*size(Z_vox,3)])    ;     ones(size(reshape(X_vox,[1 size(X_vox,1)*size(X_vox,2)*size(X_vox,3)])))  ];
    Vox_int    = ones(1,size(Vox_points,2));
    
    
    
    if tstep==1000
        if exist([folder date rec vsl 'Postproc3DReco'],'dir')~=7
            mkdir([folder date rec vsl 'Postproc3DReco'])
        end
        save([folder date rec vsl 'Postproc3DReco' vsl 'VHReduced.mat'],'Visual_hull_reduced');
    end
    
end



%% save data

% create path for processing function storage

if exist([folder date rec vsl 'Postproc3DReco'],'dir')~=7
    mkdir([folder date rec vsl 'Postproc3DReco'])
end


if 1
    save([folder date rec vsl 'Postproc3DReco' vsl 'VHReduced.mat'],'Visual_hull_reduced');
end


% Save occlusion_ratio
% save([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Occlusion_ratio.mat'],'Occlusion_ratio')

% Save occlusion_ratio
% save([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Occlusion_ratio.mat'],'Occlusion_ratio')

%%
% 

if 0
    
    vidfile = VideoWriter('visual_hull_CNN.mp4','MPEG-4');
    vidfile.FrameRate = 20; % always write in 10 fps
    % vidfile.LosslessCompression = 'true';
    open(vidfile);
    % pause

    for tstep=ctrl.tproc(1):ctrl.tproc(2)

        load([folder date rec vsl 'Preproc3DReco' vsl 'VisualHull' vsl 'tstep_' num2str(tstep) '.mat'],'Vox_int','X_vox', 'Y_vox', 'Z_vox','surf_dist')
        [fo,vo] = isosurface(X_vox,Y_vox,Z_vox,Vox_int);
        [fe,ve,ce] = isocaps(X_vox,Y_vox,Z_vox,Vox_int);
        hfig = figure(1)
        p1 = patch('Faces', fo, 'Vertices', vo);       % draw the outside of the volume
        isocolors(X_vox,Y_vox,Z_vox,surf_dist,p1)
        p1.FaceColor = 'interp';
        p1.EdgeColor = 'none';
        
        p2 = patch('Faces', fe, 'Vertices', ve, ...    % draw the end caps of the volume
           'FaceVertexCData', ce);
        p2.FaceColor = 'interp';
        p2.EdgeColor = 'none';

        camlight(40,40)                                % create two lights 
        camlight(-20,-20)
        lighting gouraud
        xlim([-5 9])
        xlim([-5 9])
        zlim([-10 70])
        axis equal
        axis tight
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        hfig.WindowState = 'maximized';
        view([1 1 1])

        Frm = getframe(gcf);

        writeVideo(vidfile, Frm);
        clf

    end

    close(vidfile)
end
















