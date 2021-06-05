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



%% 
segimsize = 256;
tsteprange = 1:10;%randi(100,1,10);
frmrange = 1:5;

% pre-allocate the trainingImages and trainingLabels matrices
% trainingImages = zeros(segimsize,segimsize,max(tsteprange),max(frmrange));
% trainingLabels = zeros(segimsize,segimsize,max(tsteprange),max(frmrange));

counter = 1;
for tstep=tsteprange
    % Load particle-free images
    Img_stack = import_frames( {folder date rec} , 'im7', tstep, 1:5);
    tstep
    
    for frame=1:frmrange(end)
        % Load sketched silhouettes
        frame
        
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
            
            for n =1:floor(768/segimsize)-2
                for m = 0:floor(1100/segimsize)-1
                    
                    % load background image
                    load([folder date rec vsl 'bgr_' num2str(frame) '.mat'], 'bgr')
                    ImTraining_nobgr = ImTraining;%-bgr.avg;
                    
                    % Segment
                    
                    ImTrainingsegmented = ImTraining_nobgr(segimsize*n + 1:segimsize*(n+1),segimsize*m + 1:segimsize*(m+1));
                    ImLabelsegmented    = ImLabel(segimsize*n + 1:segimsize*(n+1),segimsize*m + 1:segimsize*(m+1));
                    
%                     ImTrainingsegmented = abs(wiener2(medfilt2(ImTrainingsegmented)));
%                     ImTrainingsegmented = ImTrainingsegmented + min(ImTrainingsegmented(:));
                    
%                     ImTraining64x64 = ImTraining64x64 - imgaussfilt(ImTraining64x64);
% %                     
%                     ImTraining64x64 = ImTraining64x64 + abs(min(ImTraining64x64(:)));
                    
                    % Write file
%                     imwrite(uint8(ImTrainingsegmented),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
%                     imwrite(ImLabelsegmented,[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);
                    
                    % Store in matrix
                    trainingImages(:,:,counter) = ImTrainingsegmented;
                    trainingLabels(:,:,counter) = ImLabelsegmented;

                    % Increment imgfilename counter
                    counter = counter+1;
                    
                    % Plotting
                    platting=0;
                    if platting == 1
                        % plot-check
                        figure(1)
                        subplot(1,3,1)
                        imagesc(ImTrainingsegmented); axis equal; colormap gray

                        subplot(1,3,2)
                        imagesc(ImLabelsegmented); axis equal; colormap gray
                        
                        subplot(1,3,3)
                        imagesc(wiener2(medfilt2(ImTrainingsegmented))); axis equal; colormap gray
                        
%                         caxis([0 50])
                        
                        pause
                    end
                end
            end
            
        elseif frame==2
             for n = 1:floor(1104/segimsize)-2
                for m = 2%floor(334/segimsize):floor(1100/segimsize)-2
                    
                    % load background image
                    load([folder date rec vsl 'bgr_' num2str(frame) '.mat'], 'bgr')
                    ImTraining_nobgr = ImTraining;%-bgr.avg;
                    
                    % Segment
                    
                    ImTrainingsegmented = ImTraining_nobgr(segimsize*n + 1:segimsize*(n+1),segimsize*m + 1:segimsize*(m+1));
                    ImLabelsegmented    = ImLabel(segimsize*n + 1:segimsize*(n+1),segimsize*m + 1:segimsize*(m+1));
                    
%                     ImTrainingsegmented = abs(wiener2(medfilt2(ImTrainingsegmented)));
%                     ImTraining64x64 = ImTraining64x64 - imgaussfilt(ImTraining64x64);
% %                     
%                     ImTraining64x64 = ImTraining64x64 + abs(min(ImTraining64x64(:)));
                    
                    % Write file
%                     imwrite(uint8(ImTrainingsegmented),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
%                     imwrite(ImLabelsegmented,[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);

                    % Store in matrix
                    trainingImages(:,:,counter) = ImTrainingsegmented;
                    trainingLabels(:,:,counter) = ImLabelsegmented;

                    % Increment imgfilename counter
                    counter = counter+1;
                    
                    % Plotting
                    platting=0;
                    if platting == 1
                        % plot-check
                        figure(1)
                        subplot(1,2,1)
                        imagesc(ImTrainingsegmented); axis equal; colormap gray

                        subplot(1,2,2)
                        imagesc(ImLabelsegmented); axis equal; colormap gray

                        pause
                    end
                    platting=0;
                end
            end
           
        elseif frame==3
             for n = 1:2% 1:floor(766/segimsize)
                for m = 0:floor(1100/segimsize)-1
                    
%                     if counter==12
%                         pause
%                         
%                     end
                    
                    % load background image
                    load([folder date rec vsl 'bgr_' num2str(frame) '.mat'], 'bgr')
                    ImTraining_nobgr = ImTraining;%-bgr.avg;
                    
                    % Segment
                    
                    ImTrainingsegmented = ImTraining_nobgr(segimsize*n + 1:segimsize*(n+1),segimsize*m + 1:segimsize*(m+1));
                    ImLabelsegmented    = ImLabel(segimsize*n + 1:segimsize*(n+1),segimsize*m + 1:segimsize*(m+1));
                    
%                     ImTrainingsegmented = abs(wiener2(medfilt2(ImTrainingsegmented)));
                    
%                     ImTraining64x64 = ImTraining64x64 - imgaussfilt(ImTraining64x64);
% %                     
%                     ImTraining64x64 = ImTraining64x64 + abs(min(ImTraining64x64(:)));
                    
                    % Write file
%                     imwrite(uint8(ImTrainingsegmented),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
%                     imwrite(ImLabelsegmented,[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);

                    % Store in matrix
                    trainingImages(:,:,counter) = ImTrainingsegmented;
                    trainingLabels(:,:,counter) = ImLabelsegmented;

                    % Increment imgfilename counter
                    counter = counter+1;
                    
                    % Plotting
                    platting=0;
                    if platting == 1
                        % plot-check
                        figure(1)
                        subplot(1,2,1)
                        imagesc(ImTrainingsegmented); axis equal; colormap gray

                        subplot(1,2,2)
                        imagesc(ImLabelsegmented); axis equal; colormap gray

                        pause
                    end
                    platting=0;
                end
            end
            
        elseif frame==4
            for n = 1:floor(1152/segimsize)-1
                for m = 2%floor(333/segimsize):floor(1100/segimsize)-2
                    
                    % load background image
                    load([folder date rec vsl 'bgr_' num2str(frame) '.mat'], 'bgr')
                    ImTraining_nobgr = ImTraining;%-bgr.avg;
                    
                    % Segment
                    
                    ImTrainingsegmented = ImTraining_nobgr(segimsize*n + 1:segimsize*(n+1),segimsize*m + 1:segimsize*(m+1));
                    ImLabelsegmented    = ImLabel(segimsize*n + 1:segimsize*(n+1),segimsize*m + 1:segimsize*(m+1));
                    
%                     ImTrainingsegmented = abs(wiener2(medfilt2(ImTrainingsegmented)));
                    
%                     ImTraining64x64 = ImTraining64x64 - imgaussfilt(ImTraining64x64);
% %                     
%                     ImTraining64x64 = ImTraining64x64 + abs(min(ImTraining64x64(:)));
                    
                    % Write file
%                     imwrite(uint8(ImTrainingsegmented),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
%                     imwrite(ImLabelsegmented,[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);

                    % Store in matrix
                    trainingImages(:,:,counter) = ImTrainingsegmented;
                    trainingLabels(:,:,counter) = ImLabelsegmented;

                    % Increment imgfilename counter
                    counter = counter+1;
                    
                    % Plotting
                    platting=0;
                    if platting == 1
                        % plot-check
                        figure(1)
                        subplot(1,2,1)
                        imagesc(ImTrainingsegmented); axis equal; colormap gray

                        subplot(1,2,2)
                        imagesc(ImLabelsegmented); axis equal; colormap gray

                        pause
                    end
                    platting=0;
                end
            end
           
        elseif frame==5
            for n = 1:floor(1104/segimsize)-1
                for m = 2%floor(349/segimsize):floor(1100/segimsize)-1
                    
                    % load background image
                    load([folder date rec vsl 'bgr_' num2str(frame) '.mat'], 'bgr')
                    ImTraining_nobgr = ImTraining;%-bgr.avg;
                    
                    % Segment
                    
                    ImTrainingsegmented = ImTraining_nobgr(segimsize*n + 1:segimsize*(n+1),segimsize*m + 1:segimsize*(m+1));
                    ImLabelsegmented    = ImLabel(segimsize*n + 1:segimsize*(n+1),segimsize*m + 1:segimsize*(m+1));
                    
%                     ImTrainingsegmented = abs(wiener2(medfilt2(ImTrainingsegmented)));
                    
%                     ImTraining64x64 = ImTraining64x64 - imgaussfilt(ImTraining64x64);
% %                     
%                     ImTraining64x64 = ImTraining64x64 + abs(min(ImTraining64x64(:)));
                    
                    % Write file
%                     imwrite(uint8(ImTrainingsegmented),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
%                     imwrite(ImLabelsegmented,[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);

                    % Store in matrix
                    trainingImages(:,:,counter) = ImTrainingsegmented;
                    trainingLabels(:,:,counter) = ImLabelsegmented;

                    % Increment imgfilename counter
                    counter = counter+1;
                    
                    % Plotting
                    platting=0;
                    if platting == 1
                        % plot-check
                        figure(1)
                        subplot(1,2,1)
                        imagesc(ImTrainingsegmented); axis equal; colormap gray

                        subplot(1,2,2)
                        imagesc(ImLabelsegmented); axis equal; colormap gray

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



%% Image-intensity normalization across all frames, all timesteps, all image-segments




% counter=1;
for i=1:size(trainingImages,3)
    i
    % Compute image mean
    mns(i) = mean2(trainingImages(:,:,i));
    % Normalize to intensity of first frame
    intnorm_trIms(:,:,i) = trainingImages(:,:,i).*(mns(1)/mns(i));
    
    % Rotate by 90
    trIms_rot90(:,:,i) = imrotate(intnorm_trIms(:,:,i),90);
    trLbs_rot90(:,:,i) = imrotate(trainingLabels(:,:,i),90);
    
    % Rotate by 180
    trIms_rot180(:,:,i) = imrotate(intnorm_trIms(:,:,i),90*2);
    trLbs_rot180(:,:,i) = imrotate(trainingLabels(:,:,i),90*2);

    % Rotate by 270
    trIms_rot270(:,:,i) = imrotate(intnorm_trIms(:,:,i),90*3);
    trLbs_rot270(:,:,i) = imrotate(trainingLabels(:,:,i),90*3);
    
    

end




%% Write images to disk
writing=1;
if writing==1
    counter=1;

    for i=1:size(trIms_rot90,3)
        counter
        % write normal orientation to disk

        imwrite(uint8(intnorm_trIms(:,:,i)),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
        imwrite(trainingLabels(:,:,i),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);


        counter=counter+1;


        % write rot90 orientation to disk

        imwrite(uint8(trIms_rot90(:,:,i)),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
        imwrite(trLbs_rot90(:,:,i),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);


        counter=counter+1;


        % write rot180 orientation to disk
        imwrite(uint8(trIms_rot180(:,:,i)),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
        imwrite(trLbs_rot180(:,:,i),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);

        counter=counter+1;


        % write rot270 orientation to disk
        imwrite(uint8(trIms_rot270(:,:,i)),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingImages' vsl 'image_' num2str(counter) '.tiff']);
        imwrite(trLbs_rot270(:,:,i),[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Training' vsl 'trainingLabels' vsl 'labelled_image_' num2str(counter) '.tiff']);

        counter=counter+1;



    end

end