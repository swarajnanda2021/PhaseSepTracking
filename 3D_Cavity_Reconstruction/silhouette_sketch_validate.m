%% Plots the silhouette sketch edges to the real image for visual validation of the training set

clear all
close all
clc


%% Set up the correct directories
direc= cd('..\'); % one folder up
addpath(genpath([cd '\Function_Libraries']))
cd(direc) %  back to original folder
clear direc % remove variable


folder='C:\Work\PostGraduation\Scripts\ParticleTracking\Koen_testdata'; % Desktop Location
date='\Test_Swaraj'; % Date Folder
rec='\Experiment_Raw_28'; % Recording
cal='\Calibration_200514_141513'; % Calibration files


%% Begin loop for inspection

tsteprange = 1:20;
frmrange = 1:5;

invalid_cases = [];

for tstep=tsteprange
    % Load particle-free images
    Img_stack = import_frames( {folder date rec} , 'im7', tstep, 1:5);
    
    for frame=1:frmrange(end)
        
        
        
        
        ImSketch = imbinarize(double(rgb2gray(imread([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'SilhouetteSketches' vsl 'Frm_' num2str(tstep) '_Cam_' num2str(frame) '.png' ]))),230);
    
        % Calculate the edge
        Edge_B = edge((ImSketch),'canny');

        % Store edges
        edges = find(Edge_B==1);
        [edges_I,edges_J] = ind2sub(size(Edge_B),edges);


        
        
        hfig = figure(1)
        hold all
        imagesc(imadjust(Img_stack(:,:,1,frame)))
        plot(edges_J,edges_I,'.')
        colormap gray
        axis equal
        axis tight
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        hfig.WindowState = 'maximized';
        title(['tstep' num2str(tstep) ', frame' num2str(frame)])
        drawnow
        
        % Query validity from user
        validity = questdlg('Is this segmentation satisfactory?','Yes','No');
        
        switch validity
            case 'Yes'
                 invalid_cases = invalid_cases;
            case 'No'
                 invalid_cases = [invalid_cases ; [tstep frame]] ;
        end
        
        
    
    end
    

end







%% Pull out bad images and re-export them in a separate folder 

% 
% % Create path for storing bad images in ipad importable format and quality
% if exist([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'ReSketchSet_images'],'dir')~=7
%     mkdir([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'ReSketchSet_images'])
% end
% 
% 
% for set=1:size(invalid_cases,1)
%     
%     set
%     
%     Img_stack = import_frames( {folder date rec} , 'im7', invalid_cases(set,1), invalid_cases(set,2));
%     
%     
%     imwrite(imadjust(Img_stack),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'ReSketchSet_images' vsl 'Frm_' num2str(invalid_cases(set,1)) '_Cam_' num2str(invalid_cases(set,2)) '.tiff']);
% 
% 
% 
% end
% 
% 


