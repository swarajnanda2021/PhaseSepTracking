function convert_jpeg(frms)

% Converts and saves images to jpeg format after rescaling the intensity
% range of the image


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

% load calibration and define problem
[Dmap,Pmap,Dpar,Kmat,Rmat,tvec,camprop]=load_Cal;

% compute projection matrix
Pmat=cell(size(Rmat));
for c=1:size(Pmat,2)
    Pmat{c}=krtm2pmat(eye(3),Rmat{c},tvec{c});
end



%% create path for processing function storage
if exist([folder date rec vsl 'Preproc3DBubble' vsl 'BubSeg'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'BubSeg'])
end

if exist([folder date rec vsl 'Preproc3D3DBubble' vsl 'BubSeg' vsl 'Images'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'BubSeg' vsl 'Images'])
end

if exist([folder date rec vsl 'Preproc3D3DBubble' vsl 'BubSeg' vsl 'Images' vsl 'Test'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'BubSeg' vsl 'Images' vsl 'Test'])
end

if exist([folder date rec vsl 'Preproc3D3DBubble' vsl 'BubSeg' vsl 'Images' vsl 'Train'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'BubSeg' vsl 'Images' vsl 'Train'])
end

counter=1;
trainingImages =[];
for frame=1:10
    frame
    % Input frames to memory
    Img_stack = import_frames( {folder date rec} , 'im7', randi(1000,1,1), 1:5);
    
    for cam=1:5
        
%         Imfin = imadjust(uint8(Img_stack(:,:,1,cam))); % convert too uint8 for jpeg writing
        % image information
        % frame:    row     col
        % im 1 :   1:768    1:1100
        % im 2 :   1:1104   333:1100
        % im 3 :   1:766    1:1100
        % im 4 :   1:1152   333:1100
        % im 5 :   1:1104   349:1100
        
        ImTraining = double(squeeze(Img_stack(:,:,1,cam)));
        if cam==1
            ImTraining_crop = ImTraining(1:700,1:1050); 

        elseif cam==2
            ImTraining_crop = ImTraining(1:1050,333:333+699)'; 

        elseif cam==3
            ImTraining_crop = ImTraining(1:700,1:1050); 

        elseif cam==4
            ImTraining_crop = ImTraining(1:1050,333:333+699)'; 

        elseif cam==5
            ImTraining_crop = ImTraining(1:1050,349:349+699)'; 

        end


        
        
        trainingImages(:,:,counter) = ImTraining_crop;
        counter = counter+1;

%         imwrite(Imfin,[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Images' vsl 'Frm_' num2str(i) '_Cam_' num2str(cam) '.jpg']);
    end
    
     
    

end

%% Normalize and segment images
segimsize = 256;
counter=1;
intnorm_trIms=[];

for i=1:size(trainingImages,3)
    i
    % Compute image mean
    mns(i) = mean2(trainingImages(:,:,i));
    % Normalize to intensity of first frame
    intnorm_trIms(:,:,i) = trainingImages(:,:,i).*(mns(1)/mns(i));
    
    
    for n =0:floor(700/segimsize)-1
        for m = 0:floor(1050/segimsize)-1

            % load background image
            ImTraining_nobgr = intnorm_trIms(:,:,i);
            
            % Segment

            ImTrainingsegmented = ImTraining_nobgr(segimsize*n + 1:segimsize*(n+1),segimsize*m + 1:segimsize*(m+1));
            


            % Store in matrix
            trainingImages_segmented(:,:,counter) = ImTrainingsegmented;
            

            % Increment imgfilename counter
            counter = counter+1;

           
        end
    end
   
    
    
    
end


%% save
for i = 1:ceil(0.7*size(trainingImages_segmented,3))
    imwrite(uint8(trainingImages_segmented(:,:,i)),[folder date rec vsl 'Preproc3DReco' vsl 'BubSeg' vsl 'Images' vsl '\Train\Image_' num2str(i) '.jpg'])
end

for i = ceil(0.7*size(trainingImages_segmented,3))+1:size(trainingImages_segmented,3)
    imwrite(uint8(trainingImages_segmented(:,:,i)),[folder date rec vsl 'Preproc3DReco' vsl 'BubSeg' vsl 'Images' vsl '\Test\Image_' num2str(i) '.jpg'])
end


%% Label bounding boxes


imds = imageDatastore([folder date rec vsl 'Preproc3DReco' vsl 'BubSeg' vsl 'Images' vsl 'Train']);

imageLabeler(imds)




end

