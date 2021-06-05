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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% properties from the image data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prop.ext='im7'; % extension file format
prop.fps=5000; % frame per sec
[vid,prop.res]=import_frames({folder date rec},prop.ext,1,1); % image plane - number frames - number cameras resolution
dum=class(vid);
prop.bit=str2double(dum(5:end)); % dynamic range in bits
prop.roi=[1 1 prop.res(1) prop.res(2)]; % region of interest w.r.t. the calibration images
prop.ecal=cat(2,camprop.reprjavg);%+cat(2,camprop.reprjstd); % estimated reprojection error calibration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% controls for processing the data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% camera setup
ctrl.lfoc=[3 5]; % limit minimum and maximum camera focality
ctrl.ldof=[0 inf]; % limit to the depth of field [only speeds up things] [can also be done per camera] [algea can have domain -inf inf with negative focal length]
% time windowing [long time memory]
ctrl.tproc=[1 1000];%prop.res(3)]; % span of frames to process
ctrl.twin=50; % span of the window for processing
ctrl.tstep=25; % step-size for shifting the window
ctrl.rv = 10;
% spatial region of interest
% ctrl.croi=[184           1         452        1100
%            1         504        1089         845
%          248           1         537        1100
%            1         450        1152         808
%            1         549        1089         839]; % in each camera within prop.roi [speed up, not neccesssary improved qual, does have an effect on correct image segmentation glb length scl]
ctrl.croi= [1           1         756        1100
           1         350        1089        1100
           1           1         753        1100
           1         349        1152        1100
           1         368        1089        1100];

% resolution of the data
ctrl.tres=1; % time resolution in frames w.r.t. average object displacement
ctrl.rbox=[1 1]; % spatial resolution w.r.t. smallest pixel scale sampling objects (also allows increase resolution by 1/n if signal to noise allows)

% multiple scale image processing [fwin and dres set scale table with average growth, e.g. x 1.2]
ctrl.nscl=1; % number of scales in image
ctrl.fwin=5; % number of window sizes to filters % for diffraction limited
% ctrl.fwin=7; % number of window sizes to filters % for bubbles
ctrl.dres=[1 1]; % box filter to downgrade resolution image-plane
ctrl.brightsel = 'true'; % caps the total particles considered by brightness, to reduce RAM consumption
ctrl.nparts = [9000 9000 9000 9000 9000]; % Bright particle count upper bound % for diffraction limited only
% ctrl.nparts = [500 500 500 500 500];
ctrl.nKernel = 11; % Local correlation box size

% tracking and tracking matching
ctrl.dspl=2; % displacement by bodylength in image-plane
ctrl.dspr=0.5; % camera disparity in body-length
ctrl.ovlp=0.5; % percentage of overlap bounding ellipse

% trajectories shape filter and prediction [short time memory, keep it small - Kalman filter-ish]
ctrl.tfit=3; % stencil size trajectory fit [filter quality by # frames, e.g. 5, 10, 15]
ctrl.ord=1; % polynomial order [filter degrees of freedom, i.e. 0-position, 1-velocity, 2-accelaration, 3-minimal curvature]
ctrl.prd=2; % allowed timesteps to predict path [skip frames, e.g. 1, 2, 3]

%% save proc set
disp('Save properties and controls')
save([folder date rec vsl 'Preproc3DTracking' vsl 'prop'],'prop'); % save properties
save([folder date rec vsl 'Preproc3DTracking' vsl 'ctrl'],'ctrl'); % save controls

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