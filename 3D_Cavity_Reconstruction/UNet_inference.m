clear all
close all
clc

%% Code that prepares the training image set, background images and the cavity silhouette

direc= cd('..\'); % one folder up
addpath(genpath([cd '\Function_Libraries']))
cd(direc) %  back to original folder
clear direc % remove variable


folder='C:\Work\PostGraduation\Scripts\ParticleTracking\Koen_testdata'; % Desktop Location
date='\Test_Swaraj'; % Date Folder
rec='\Experiment_Raw_28'; % Recording
cal='\Calibration_200514_141513'; % Calibration files


% Load the DL network
% load('Training_27012021_batchsize96_epochs200_ResNet18_titanx.mat')
% load('Training_28012021_batchsize40_epochs200_ResNet50_titanx_optionset_2.mat')
load('Training_29012021_batchsize48_epochs50_MobileNet50_titanx_optionset_2.mat')
%% Perform inference over a loop

% v = VideoWriter('ResNet18_6.avi','uncompressed AVI');
% v.FrameRate = 10;
framenum=4;
save_imgmean2 = [];
% v.LosslessCompression = 'true';
% open(v)
counter = 1;
for tstep = 1:150
    % Load the frame 
    Img_stack_ref = imread(['C:\Work\PostGraduation\Scripts\ParticleTracking\Koen_testdata\Test_Swaraj\Experiment_Raw_28\Preproc3DReco\ImgSeg\Training\trainingImages\image_1.tiff']);
    Img_stack = import_frames({folder date rec},'im7',tstep,framenum);
    
%     % Test ROI selection
%     ROI = [131 500];
%     
%     Img_stack = (Img_stack(ROI(1):ROI(1)+256,ROI(2):ROI(2)+256));
%     
%     save_imgmean2(counter) = mean2(Img_stack);
    counter = counter + 1;
    
    if framenum==1
        Img_stack(Img_stack>250) = 0;
        Img_stack(Img_stack<50) = 0;

    elseif framenum==2
        Img_stack(Img_stack>80)=0;
        Img_stack(Img_stack<10)=0;

    elseif framenum==3
        Img_stack(Img_stack>100)=0;
        Img_stack(Img_stack<20)=0;

    elseif framenum==4
        Img_stack(Img_stack>250)=0;
        Img_stack(Img_stack<40)=0;

    elseif framenum==5
        Img_stack(Img_stack>100)=0;
        Img_stack(Img_stack<20)=0;

    end
    
    
    
    % Normalize image
%     Img_stack = (Img_stack'.*mean2(Img_stack_ref)/mean2(Img_stack));

    % Perform Inference
    C = semanticseg((abs((Img_stack)')), net);
    B = labeloverlay(Img_stack',C);
    
    
    
    fun = @(img_stac) labeloverlay(img_stac.data.*mean2(Img_stack_ref)/mean2(img_stac.data),semanticseg((abs((img_stac.data).*mean2(Img_stack_ref)/mean2(img_stac.data))), net));
    
    

    I2 = blockproc(Img_stack(1:1024,1:1024),[256 256],fun,'TrimBorder', false);

%     imagesc(I2)
    
    % Binarize the image
    Segmented_img = zeros(size(Img_stack));
    
    Segmented_img(1:1024,1:1024) = rgb2gray(I2);
    
    Segmented_img(Segmented_img>15) = 0;
    
    Segmented_img(Segmented_img==15) = 1;
    
    
    
%     
%     
%     
    % Calculate the edge
    Edge_B = edge(rgb2gray(I2),'canny');
    
    % Store edges
    edges = find(Edge_B==1);
    [edges_I,edges_J] = ind2sub(size(Edge_B),edges);
    
    hfig = figure(1)
    hold all
    imagesc((Img_stack))
    plot(edges_J,edges_I,'r.')
    colormap gray
    axis equal
    axis tight
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    hfig.WindowState = 'maximized';
    drawnow
%     pause




%     
%     % Plot
%     hfig = figure(1)
%     subplot(1,2,1)
%     imagesc(flipud(B))
%     
%     colormap gray
%     axis equal
%     axis tight
%     set(gca,'xtick',[])
%     set(gca,'xticklabel',[])
%     set(gca,'ytick',[])
%     set(gca,'yticklabel',[])
%     
%     subplot(1,2,2)
%     imagesc(flipud(Img_stack(1:768,:)'))
% 
%     colormap gray
%     axis equal
%     axis tight
%     set(gca,'xtick',[])
%     set(gca,'xticklabel',[])
%     set(gca,'ytick',[])
%     set(gca,'yticklabel',[])
%     
%     hfig.WindowState = 'maximized';
%     
%     drawnow
    
%     
%     frame = getframe(gcf);
%     writeVideo(v,frame)
    
end

% close(v)





