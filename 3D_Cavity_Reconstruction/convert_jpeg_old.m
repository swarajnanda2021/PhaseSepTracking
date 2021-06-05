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
% folder='C:\Work\PostGraduation\Scripts'; % Desktop Location
% date='\ParticleTracking'; % Date Folder
% rec='\Measurement37'; % Recording
% cal='\Calibration_200514_141513'; % Calibration files

% load calibration and define problem
[Dmap,Pmap,Dpar,Kmat,Rmat,tvec,camprop]=load_Cal;

% compute projection matrix
Pmat=cell(size(Rmat));
for c=1:size(Pmat,2)
    Pmat{c}=krtm2pmat(eye(3),Rmat{c},tvec{c});
end



%% create path for processing function storage
if exist([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg'])
end

if exist([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Images'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Images'])
end

if exist([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'SilhouetteData'],'dir')~=7
    mkdir([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'SilhouetteData'])
end

for i=1:max(frms)

    % Input frames to memory
    Img_stack = import_frames( {folder date rec} , 'im7', i, 1:5);
    
    for cam=1:5
        
        Imfin = imadjust(uint8(Img_stack(:,:,1,cam))); % convert too uint8 for jpeg writing
        
        
%         imshow(imadjust(Img_stack(:,:,1,cam)))
%         drawnow
        
%         pause
        
        
%         bgr = load([folder date rec vsl 'bgr_',num2str(cam),'.mat'],'bgr');
%         Imsub = Im-bgr.bgr.avg;

%         if cam==1
% %             Im(Im>200) = 0;
% %             Im(Im<100) = 0;
%             Imfin =  imminmax(Imsub,strel('disk',10,0),0.03); %ctrl.fwin(end)
%         elseif cam==2
% %             Imfin(Imfin>120) = 0;
% %             Imfin(Imfin<10) = 0;
%             Imfin =  imminmax(Imsub,strel('disk',10,0),0.06); %ctrl.fwin(end)
%         elseif cam==3
% %             Imfin(Imfin>120) = 0;
% %             Imfin(Imfin<10) = 0;
%             Imfin =  imminmax(Imsub,strel('disk',12,0),0.005); %ctrl.fwin(end)    
%         elseif cam==4
% %             Imfin(Imfin>250) = 0;
% %             Imfin(Imfin<20) = 0;
%             Imfin =  imminmax(Imsub,strel('disk',10,0),0.03); %ctrl.fwin(end)               
%             
%         elseif cam==5
% %             Imfin(Imfin>150) = 0;
% %             Imfin(Imfin<10) = 0;
%             Imfin =  imminmax(Imsub,strel('disk',5,0),0.001); %ctrl.fwin(end)               
%             
%         end
        imwrite(Imfin,[folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'Images' vsl 'Frm_' num2str(i) '_Cam_' num2str(cam) '.jpg']);
    end
    
     
    

end





end

