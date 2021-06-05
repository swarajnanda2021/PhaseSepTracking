function [prop,ctrl] = proc_set
%proc_set Processing settings
%
%   Output:
%       prop - Properties
%       ctrl - Control
%
%   Full commenting first data set in switch case
%   Play around with settings to optimize, tip: keep it simple.
%
%   Note: future releases can make parameters vary per camera, or let the
%       cameras find settings itself.
                

%% global
global folder date rec camprop %plotting

%% Specific


prop.modelname = 'noEmpty_thick_optionset3_augmented_rand_resnet18_07042021';
prop.range_x = [-10  10];
prop.range_y = [-10  10];
prop.range_z = [-10  50];
prop.del_vox = 0.05 ; % delta space in millimeters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% properties from the image data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prop.ext='im7'; % extension file format
prop.fps=5000; % frame per sec
[vid,prop.res]=import_frames({folder date rec},prop.ext,1,1); % image plane - number frames - number cameras resolution
dum=class(vid);
prop.bit=str2double(dum(5:end)); % dynamic range in bits
prop.roi=[1 1 prop.res(1) prop.res(2)]; % region of interest w.r.t. the calibration images
prop.ecal=cat(2,camprop.reprjavg)+cat(2,camprop.reprjstd); % estimated reprojection error calibration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% controls for processing the data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% camera setup
ctrl.lfoc=[1 2 3 4 5]; % choose the cameras you want to use for the reconstruction
% time windowing [long time memory]
ctrl.tproc=[1 3000];%prop.res(3)]; % span of frames to process (no window processing needed)


% spatial region of interest (for occlusion metric)
ctrl.croi= [1           1         756        1100
           1         350        1089        1100
           1           1         753        1100
           1         349        1152        1100
           1         368        1089        1100];


%% save proc set
disp('Save properties and controls')
save([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'prop'],'prop'); % save properties
save([folder date rec vsl 'Preproc3DReco' vsl 'ImgSeg' vsl 'ctrl'],'ctrl'); % save controls

end

% % Code snippet to select croi by clicking and figure inspection
% if strcmp(plotting,'on') % to select ctrl.croi
%     for c=1:prop.res(4) % loop views
%         
%         % make figure
%         figure(c)
%         
%         for n=ctrl.tproc % loop first and last frame
%             
%             % load image
%             img=double(import_frames({folder date rec},prop.ext,n,c));
%             
%             % plot
%             plt=log(img+1);
%             surf(0*plt',plt')
%             view(2)
%             shading flat
%             axis tight
%             axis equal
%             colormap gray
%             hold on
%             caxis([mean(plt(:))-3*std(plt(:)) mean(plt(:))+6*std(plt(:))])
%             
%             % draw existing roi
%             roi=ctrl.croi(c,:);
%             roi(3:4)=roi(3:4)-roi(1:2);
%             h = imrect(gca,roi);
%             
%             % pause
%             title(['Camera ',num2str(c),' Frame ',num2str(n),' Modify rectagle to ROI, and press spacebar'])
%             pause
%             
%             % get adjusted new roi
%             title(['Camera ',num2str(c),' Frame ',num2str(n),' Selected ROI'])
%             roi=round(getPosition(h));
%             roi(3:4)=roi(3:4)+roi(1:2);
%             ctrl.croi(c,:)=roi;
%             
%         end % n
%         
%     end % c
% end