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
rec='\Measurement37'; % Recording
cal='\Calibration_200514_141513'; % Calibration files

dataSetDir = [folder date rec '\Preproc3DReco\BubSeg'];
imageDir = [dataSetDir '\Images'];
labelDir = [dataSetDir '\Labels'];



%% Network model


data = load('fasterRCNNVehicleTrainingData.mat');
lgraph = layerGraph(data.detector.Network);

% vehicleDetector = load('yolov2VehicleDetector.mat');
% lgraph = vehicleDetector.lgraph

% analyzeNetwork(lgraph)




load([labelDir '\test2.mat'])
% 
% imds = imageDatastore(gTruth.DataSource.Source);
% 
% blds = boxLabelDatastore(gTruth.LabelData);

[imds,blds] = objectDetectorTrainingData(gTruth); 

cds = combine(imds,blds);


options = trainingOptions('sgdm', ...
       'InitialLearnRate', 0.001, ...
       'Verbose',true, ...
       'MiniBatchSize',2, ...
       'MaxEpochs',30,'Plots', 'training-progress', ...
    'ExecutionEnvironment', 'gpu'); 
   
[detector,info] = trainFasterRCNNObjectDetector(cds,lgraph,options);
    
    
%% test

 
%     vidfile = VideoWriter('bubble_detector_FasterRCNN.mp4','MPEG-4');
%     vidfile.FrameRate = 5; % always write in 10 fps
%     % vidfile.LosslessCompression = 'true';
%     open(vidfile);

for i=1:100
    I = imread([imageDir '\Image_' num2str(i) '.jpg']);


    [bboxes,scores,labels] = detect(detector,I)


    detectedI = insertObjectAnnotation(I,'Rectangle',bboxes,cellstr(labels));
    figure(1)
    imshow(detectedI)
    drawnow
    
%     Frm = getframet(gcf);

%     writeVideo(vidfile, Frm);
    pause
    clf
    
end

% close(vidfile)





    

