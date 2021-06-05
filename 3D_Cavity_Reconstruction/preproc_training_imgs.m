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
rec='\Experiment_Raw_23'; % Recording
cal='\Calibration_200514_141513'; % Calibration files

%% create path for processing function storage


if exist([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg'])
end

% Folder for storing particle removed images 
if exist([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'TrainingImages'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'TrainingImages'])
end

% Folder for storing binarized background images
if exist([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'LabelledImages'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'LabelledImages'])
end

finalsize=256;

%% 
tsteprange = 1:20;
frmrange = 1:5;
counter = 1;
for tstep=1%:tsteprange(end)
    % Load particle-free images
    Img_stack = import_frames( {folder date rec} , 'im7', tstep, 1:5);
    
    
    for frame=1%1:frmrange(end)
        % Load sketched silhouettes
        
        
        %% Part 1: Convert training images to scaled format of size 1024X1024 and good intensity range
        ImTraining = double(squeeze(Img_stack(:,:,1,frame)));
        if frame==1
            ImTraining(ImTraining>250)=0;
            ImTraining(ImTraining<50)=0;

        elseif frame==2
            ImTraining(ImTraining>80)=0;
            ImTraining(ImTraining<10)=0;

        elseif frame==3
            ImTraining(ImTraining>100)=0;
            ImTraining(ImTraining<20)=0;

        elseif frame==4
            ImTraining(ImTraining>250)=0;
            ImTraining(ImTraining<40)=0;

        elseif frame==5
            ImTraining(ImTraining>100)=0;
            ImTraining(ImTraining<20)=0;

        end
        
        % apply a basic median filter to remove particles from the image
        % (salt and pepper noise)
%         ImTraining = wiener2(medfilt2(ImTraining));
%         
%         ImTraining = abs(ImTraining./max(ImTraining(:)));
%         ImTraining = (imgaussfilt(imgaussfilt(ImTraining,[2 2]) - imgaussfilt(ImTraining)));
%         ImTraining = ImTraining + abs(min(ImTraining(:)));
%         
% %         
%         ImTraining = ImTraining + abs(min(ImTraining(:)));
%         
%         pause
        
        %% Part 2: Convert sketched background images to scaled format of size 1024X1024
        % load and binarize
        ImSketch = imbinarize(double(rgb2gray(imread([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'SilhouetteSketches' vsl 'Frm_' num2str(tstep) '_Cam_' num2str(frame) '.png' ]))),230);
                
%         pause
        
%         ImSketch = double(imcomplement(imbinarize(ImSketch,200))); % use imcomplement if you want this line to give you non-zero boundary and zero cavity intensity
        
        ImLabel = ImSketch;
        
        %% Break-up ImLabel and ImTraining into images of 64X64 pixels
        
%         imagesc(ImTraining(64*n + 1:64*(n+1),64*n + 1:64*(n+1))); axis equal; colormap gray
        % image information
        % frame:    row     col
        % im 1 :   1:768    1:1100
        % im 2 :   1:1104   333:1100
        % im 3 :   1:766    1:1100
        % im 4 :   1:1152   333:1100
        % im 5 :   1:1104   349:1100
        
        if frame==1
            
            for n =2:floor(768/64)-5
                for m = 2:floor(1100/64)-2
                    % Segment
                    
                    ImTraining64x64 = ImTraining(64*n + 1:64*(n+1),64*m + 1:64*(m+1));
                    ImLabel64x64    = ImLabel(64*n + 1:64*(n+1),64*m + 1:64*(m+1));
                    
                    
%                     ImTraining64x64 = ImTraining64x64 - imgaussfilt(ImTraining64x64);
% %                     
%                     ImTraining64x64 = ImTraining64x64 + abs(min(ImTraining64x64(:)));
                    
                    % Write file
                    imwrite(uint8(ImTraining64x64),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
                    imwrite(ImLabel64x64,[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);
                    
                    % Increment imgfilename counter
                    counter = counter+1;
                    
                    % Plotting
                    platting=0;
                    if platting == 1
                        % plot-check
                        figure(1)
                        subplot(1,2,1)
                        imagesc(ImTraining(64*n + 1:64*(n+1),64*m + 1:64*(m+1))); axis equal; colormap gray

                        subplot(1,2,2)
                        imagesc(ImLabel(64*n + 1:64*(n+1),64*m + 1:64*(m+1))); axis equal; colormap gray

                        pause
                    end
                end
            end
            
        elseif frame==2
             for n = 2:floor(1104/64)-2
                for m = floor(334/64)+2:floor(1100/64)-2
                    % Segment
                    
                    ImTraining64x64 = ImTraining(64*n + 1:64*(n+1),64*m + 1:64*(m+1));
                    ImLabel64x64    = ImLabel(64*n + 1:64*(n+1),64*m + 1:64*(m+1));
                    
                    % Write file
                    imwrite(uint8(ImTraining64x64),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
                    imwrite(uint8(ImLabel64x64),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);
                    
                    % Increment imgfilename counter
                    counter = counter+1;
                    
                    % Plotting
                    platting=0;
                    if platting == 1
                        % plot-check
                        figure(1)
                        subplot(1,2,1)
                        imagesc(ImTraining(64*n + 1:64*(n+1),64*m + 1:64*(m+1))); axis equal; colormap gray

                        subplot(1,2,2)
                        imagesc(ImLabel(64*n + 1:64*(n+1),64*m + 1:64*(m+1))); axis equal; colormap gray

                        pause
                    end
                end
            end
           
        elseif frame==3
             for n = 2:floor(766/64)-2
                for m = 2:floor(1100/64)-2
                    % Segment
                    
                    ImTraining64x64 = ImTraining(64*n + 1:64*(n+1),64*m + 1:64*(m+1));
                    ImLabel64x64    = ImLabel(64*n + 1:64*(n+1),64*m + 1:64*(m+1));
                    
                    % Write file
                    imwrite(uint8(ImTraining64x64),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
                    imwrite(uint8(ImLabel64x64),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);
                    
                    % Increment imgfilename counter
                    counter = counter+1;
                    
                    % Plotting
                    platting=0;
                    if platting == 1
                        % plot-check
                        figure(1)
                        subplot(1,2,1)
                        imagesc(ImTraining(64*n + 1:64*(n+1),64*m + 1:64*(m+1))); axis equal; colormap gray

                        subplot(1,2,2)
                        imagesc(ImLabel(64*n + 1:64*(n+1),64*m + 1:64*(m+1))); axis equal; colormap gray

                        pause
                    end
                end
            end
            
        elseif frame==4
            for n = 2:floor(1152/64)-2
                for m = floor(333/64)+2:floor(1100/64)-2
                    % Segment
                    
                    ImTraining64x64 = ImTraining(64*n + 1:64*(n+1),64*m + 1:64*(m+1));
                    ImLabel64x64    = ImLabel(64*n + 1:64*(n+1),64*m + 1:64*(m+1));
                    
                    % Write file
                    imwrite(uint8(ImTraining64x64),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
                    imwrite(uint8(ImLabel64x64),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);
                    
                    % Increment imgfilename counter
                    counter = counter+1;
                    
                    % Plotting
                    platting=0;
                    if platting == 1
                        % plot-check
                        figure(1)
                        subplot(1,2,1)
                        imagesc(ImTraining(64*n + 1:64*(n+1),64*m + 1:64*(m+1))); axis equal; colormap gray

                        subplot(1,2,2)
                        imagesc(ImLabel(64*n + 1:64*(n+1),64*m + 1:64*(m+1))); axis equal; colormap gray

                        pause
                    end
                end
            end
           
        elseif frame==5
            for n = 2:floor(1104/64)-2
                for m = floor(349/64)+2:floor(1100/64)-2
                    % Segment
                    
                    ImTraining64x64 = ImTraining(64*n + 1:64*(n+1),64*m + 1:64*(m+1));
                    ImLabel64x64    = ImLabel(64*n + 1:64*(n+1),64*m + 1:64*(m+1));
                    
                    % Write file
                    imwrite(uint8(ImTraining64x64),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
                    imwrite(uint8(ImLabel64x64),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);
                    
                    % Increment imgfilename counter
                    counter = counter+1;
                    
                    % Plotting
                    platting=0;
                    if platting == 1
                        % plot-check
                        figure(1)
                        subplot(1,2,1)
                        imagesc(ImTraining(64*n + 1:64*(n+1),64*m + 1:64*(m+1))); axis equal; colormap gray

                        subplot(1,2,2)
                        imagesc(ImLabel(64*n + 1:64*(n+1),64*m + 1:64*(m+1))); axis equal; colormap gray

                        pause
                    end
                end
            end
           
        end
        
        
        %% Write images to disk, update counter and clear images
        % export
%         imwrite(ImTraining,[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.png']);
        
%         imwrite(uint8(ImLabel+ImSketch),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.png']);
% 
%         counter = counter+1;
%         pause
        
    end
    
    
    
    
end











