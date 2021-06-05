clear all
close all
clc



folder='C:\Work\PostGraduation\Scripts'; % Desktop Location
date='\ParticleTracking'; % Date Folder
rec='\Measurement37'; % Recording

dataSetDir = [folder date rec '\Preproc3DReco\ImgSeg\Training'];
imageDir = [dataSetDir '\trainingImages'];
labelDir = [dataSetDir '\trainingLabels'];

testImagesDir = [dataSetDir '\testImages'];
testLabelsDir = [dataSetDir '\testLabels'];



classNames = ["cavity","background"];
labelIDs   = [0 1];

finalsize = 256;
% Create Unet 32by32 is the example
imageSize = [finalsize finalsize 3];
numClasses = 2;
encoderDepth = 5;

% Model number 1: U-Net
% lgraph = unetLayers(imageSize, numClasses, 'EncoderDepth',encoderDepth);
% Model number 2: ResNet-18
lgraph = deeplabv3plusLayers(imageSize, numClasses, "resnet18");
% Model number 3: ResNet-50
% lgraph = deeplabv3plusLayers(imageSize, numClasses, "resnet50");
% Model number 4: MobileNet-50
% lgraph = deeplabv3plusLayers(imageSize, numClasses, "mobilenetv2");

% Change loss function of the network
layer = dicePixelClassificationLayer('Name','pixclass');

lgraph = replaceLayer(lgraph,'classification',layer);



% Create a datastore for training the network.
imds = imageDatastore(imageDir);
pxds = pixelLabelDatastore(labelDir,classNames,labelIDs);


% Augment images around the rotation with some added scaling
imageAugmenter = imageDataAugmenter( ...
    'RandRotation',[-45,45],'RandScale', [0.7 1.3],'RandXReflection',1,...
    'RandYReflection',1, 'RandXShear',[-10 10], 'RandYShear',[-10 10]);







augds = pixelLabelImageDatastore(imds,pxds,'DataAugmentation',imageAugmenter,'OutputSize',imageSize,...
    'ColorPreprocessing','gray2rgb');


ds = pixelLabelImageDatastore(imds,pxds,'OutputSize',imageSize,...
    'ColorPreprocessing','gray2rgb');




% %Set training options 1

% options = trainingOptions('sgdm', ...
%     'InitialLearnRate',1e-3, ...
%     'MaxEpochs',50, ...
%     'VerboseFrequency',1, ...
%     'MiniBatchSize',96, ...
%     'Plots', 'training-progress', ...
%     'ExecutionEnvironment', 'gpu');

% Set training options 2

% options = trainingOptions('sgdm', ...
%     'LearnRateSchedule','piecewise', ...
%     'LearnRateDropFactor',0.2, ...
%     'LearnRateDropPeriod',5, ...
%     'InitialLearnRate',0.02, ...
%     'MaxEpochs',50, ...
%     'VerboseFrequency',1, ...
%     'MiniBatchSize',48, ...
%     'Plots', 'training-progress', ...
%     'ExecutionEnvironment', 'gpu');



options = trainingOptions('sgdm', ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.94, ...
    'LearnRateDropPeriod',2, ...
    'InitialLearnRate',0.05, ...
    'L2Regularization', 5e-4 , ...
    'MaxEpochs',50, ...
    'Shuffle', 'every-epoch',...
    'VerboseFrequency',1, ...
    'MiniBatchSize',96, ...
    'Plots', 'training-progress', ...
    'ExecutionEnvironment', 'gpu');


%% Training
% Train the network
[net,info] = trainNetwork(augds,lgraph,options)



if 1
    
    save([folder vsl date vsl 'Deep_Learning_Model_Database' vsl modelname '.mat'])

end


%% Evaluate semantic segmentation network 


counter = 1;
for i=1001:1600
    
    testim = imread([testImagesDir '\image_' num2str(i) '.tiff']);
    C = semanticseg(abs(testim), net);
    B = labeloverlay(testim,C);
    
    
    labelim = imread([testLabelsDir '\labelled_image_' num2str(i) '.tiff']);
    
    

    
    
    similarity_bfscore(counter) = bfscore(double(C)-1, double(labelim));
    similarity_jaccardindex(counter) = jaccard(double(C)-1, double(labelim));
    
    
    counter = counter+1;
    
    figure(1)
    
%     subplot(1,3,1)
    imagesc(B)

    colormap gray
    axis equal
    axis tight
    title(num2str((i-1)/4 + 1))
%     subplot(1,3,2)
% %     figure(2)
%     imagesc(testim)
% 
%     colormap gray
%     axis equal
%     axis tight
%     
%     subplot(1,3,3)
% %     figure(3)
%     imagesc(labelim)
% 
%     colormap gray
%     axis equal
%     axis tight
%     title(['Image:' num2str(i)])
%     % 
    drawnow
%     pause
end






%%
% 
% layer = 1;
% channels = 1;
% 
% I = deepDreamImage(net,layer,channels);
% 
% figure
% for i = 1:25
%     subplot(5,5,i)
%     imshow(I(:,:,:,i))
% end
% 



